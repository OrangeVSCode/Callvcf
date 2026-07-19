#!/usr/bin/env python3
"""Lead-variant LD, annotation and regional association helpers."""

import csv
import gzip
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

from vcf_service import VCFError, parse_loci


def genotype_dosage(gt):
    value = str(gt or ".").split(":", 1)[0].replace("|", "/")
    if "." in value:
        return None
    alleles = value.split("/")
    if len(alleles) != 2 or any(x not in {"0", "1"} for x in alleles):
        return None
    return sum(int(x) for x in alleles)


def pairwise_r2(left, right, min_samples=20):
    pairs = [(a, b) for a, b in zip(left, right) if a is not None and b is not None]
    n = len(pairs)
    if n < min_samples:
        return None, n
    mean_a = sum(a for a, _ in pairs) / n
    mean_b = sum(b for _, b in pairs) / n
    var_a = sum((a - mean_a) ** 2 for a, _ in pairs)
    var_b = sum((b - mean_b) ** 2 for _, b in pairs)
    if var_a <= 0 or var_b <= 0:
        return None, n
    cov = sum((a - mean_a) * (b - mean_b) for a, b in pairs)
    return max(0.0, min(1.0, (cov * cov) / (var_a * var_b))), n


def _open_text(path):
    with Path(path).open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    if compressed:
        return gzip.open(str(path), "rt", encoding="utf-8", errors="replace")
    return Path(path).open("r", encoding="utf-8", errors="replace")


def _existing_file(path_text, label):
    path = Path(os.path.expandvars(os.path.expanduser(str(path_text or "")))).resolve()
    if not path.is_file():
        raise VCFError("{}不存在：{}".format(label, path))
    return path


def extract_info_annotations(record):
    info = record.get("info") or {}
    output = []
    for raw in str(info.get("ANN", "")).split(","):
        if not raw:
            continue
        values = raw.split("|")
        output.append({
            "source": "ANN/SnpEff", "allele": values[0] if values else "",
            "effect": values[1] if len(values) > 1 else "",
            "impact": values[2] if len(values) > 2 else "",
            "gene": values[3] if len(values) > 3 else "",
            "gene_id": values[4] if len(values) > 4 else "",
            "feature": values[6] if len(values) > 6 else "",
            "hgvs_c": values[9] if len(values) > 9 else "",
            "hgvs_p": values[10] if len(values) > 10 else "",
            "raw": raw,
        })
    for raw in str(info.get("CSQ", "")).split(","):
        if not raw:
            continue
        values = raw.split("|")
        output.append({
            "source": "CSQ/VEP", "allele": values[0] if values else "",
            "effect": values[1] if len(values) > 1 else "",
            "impact": values[2] if len(values) > 2 else "",
            "gene": values[3] if len(values) > 3 else "",
            "gene_id": values[4] if len(values) > 4 else "",
            "feature": values[6] if len(values) > 6 else "",
            "hgvs_c": "", "hgvs_p": "", "raw": raw,
        })
    for raw in str(info.get("BCSQ", "")).split(","):
        if raw:
            values = raw.split("|")
            output.append({
                "source": "BCSQ/bcftools", "allele": "",
                "effect": values[0] if values else "", "impact": "",
                "gene": values[1] if len(values) > 1 else "",
                "gene_id": "", "feature": values[2] if len(values) > 2 else "",
                "hgvs_c": values[3] if len(values) > 3 else "",
                "hgvs_p": values[4] if len(values) > 4 else "", "raw": raw,
            })
    return output


def _norm_header(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def annotation_tsv_matches(path_text, chrom, pos):
    if not path_text:
        return []
    path = _existing_file(path_text, "功能注释表")
    with _open_text(path) as handle:
        first = handle.readline()
        if not first:
            return []
        delimiter = "\t" if first.count("\t") >= first.count(",") else ","
        header = [x.strip() for x in first.rstrip("\r\n").split(delimiter)]
        normalized = [_norm_header(x) for x in header]
        chr_i = next((i for i, x in enumerate(normalized) if x in {"chr", "chrom", "chromosome"}), None)
        pos_i = next((i for i, x in enumerate(normalized) if x in {"pos", "position", "bp", "start"}), None)
        if chr_i is None or pos_i is None:
            raise VCFError("功能注释表必须包含 chr/chrom 与 pos/position 列")
        results = []
        for line in handle:
            values = line.rstrip("\r\n").split(delimiter)
            if len(values) <= max(chr_i, pos_i):
                continue
            try:
                matches = values[chr_i].strip() == str(chrom) and int(float(values[pos_i])) == int(pos)
            except ValueError:
                continue
            if matches:
                results.append({header[i]: values[i] if i < len(values) else "" for i in range(len(header))})
                if len(results) >= 100:
                    break
        return results


def _parse_attributes(text):
    attrs = {}
    for token in str(text).strip().strip(";").split(";"):
        token = token.strip()
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
        elif " " in token:
            key, value = token.split(None, 1)
        else:
            continue
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def locate_gene(path_text, chrom, pos):
    if not path_text:
        return {"status": "not_configured", "matches": []}
    path = _existing_file(path_text, "GFF3/GTF")
    genes, overlapping_features = [], []
    with _open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 9 or parts[0] != str(chrom):
                continue
            try:
                start, end = int(parts[3]), int(parts[4])
            except ValueError:
                continue
            feature_type, strand = parts[2], parts[6]
            attrs = _parse_attributes(parts[8])
            name = (attrs.get("gene_name") or attrs.get("Name") or attrs.get("gene") or
                    attrs.get("gene_id") or attrs.get("ID") or attrs.get("Parent") or "unknown")
            item = {"type": feature_type, "start": start, "end": end, "strand": strand,
                    "name": name, "id": attrs.get("gene_id") or attrs.get("ID") or name,
                    "parent": attrs.get("Parent") or attrs.get("transcript_id") or ""}
            if feature_type.lower() in {"gene", "pseudogene"}:
                genes.append(item)
            if start <= pos <= end and feature_type.lower() in {
                "gene", "pseudogene", "mrna", "transcript", "exon", "cds",
                "five_prime_utr", "three_prime_utr", "utr"
            }:
                overlapping_features.append(item)

    containing = [g for g in genes if g["start"] <= pos <= g["end"]]
    priority = {"cds": 0, "five_prime_utr": 1, "three_prime_utr": 1, "utr": 1, "exon": 2,
                "mrna": 3, "transcript": 3, "gene": 4, "pseudogene": 4}
    feature = min(overlapping_features, key=lambda x: priority.get(x["type"].lower(), 9), default=None)
    matches = []
    if containing:
        for gene in containing[:20]:
            relation = "gene_body/intron"
            if feature:
                kind = feature["type"].lower()
                relation = {"cds": "CDS", "exon": "exon", "five_prime_utr": "5'UTR",
                            "three_prime_utr": "3'UTR", "utr": "UTR"}.get(kind, relation)
            matches.append({**gene, "relation": relation, "distance": 0,
                            "overlapping_feature": feature})
    else:
        ranked = sorted(genes, key=lambda g: min(abs(pos - g["start"]), abs(pos - g["end"])))[:5]
        for gene in ranked:
            if pos < gene["start"]:
                distance = gene["start"] - pos
                relation = "upstream" if gene["strand"] != "-" else "downstream"
            else:
                distance = pos - gene["end"]
                relation = "downstream" if gene["strand"] != "-" else "upstream"
            matches.append({**gene, "relation": relation, "distance": distance,
                            "overlapping_feature": None})
    return {"status": "ok", "path": str(path), "matches": matches}


def domain_matches(path_text, identifiers):
    if not path_text:
        return []
    path = _existing_file(path_text, "结构域注释表")
    wanted = {str(x).strip().lower() for x in identifiers if str(x).strip()}
    if not wanted:
        return []
    with _open_text(path) as handle:
        first = handle.readline()
        if not first:
            return []
        delimiter = "\t" if first.count("\t") >= first.count(",") else ","
        header = [x.strip() for x in first.rstrip("\r\n").split(delimiter)]
        normalized = [_norm_header(x) for x in header]
        key_indices = [i for i, x in enumerate(normalized) if x in {
            "gene", "geneid", "genename", "protein", "proteinid", "transcript", "transcriptid"
        }]
        if not key_indices:
            key_indices = [0]
        output = []
        for line in handle:
            values = line.rstrip("\r\n").split(delimiter)
            keys = {values[i].strip().lower() for i in key_indices if i < len(values)}
            if wanted.intersection(keys):
                output.append({header[i]: values[i] if i < len(values) else "" for i in range(len(header))})
                if len(output) >= 500:
                    break
        return output


def _ps_record(parts, header):
    norm = [_norm_header(x) for x in header] if header else []
    chr_i = next((i for i, x in enumerate(norm) if x in {"chr", "chrom", "chromosome"}), None)
    pos_i = next((i for i, x in enumerate(norm) if x in {"pos", "position", "bp"}), None)
    p_i = next((i for i, x in enumerate(norm) if x in {"p", "pvalue", "pval", "pvaluewald"}), None)
    marker_i = next((i for i, x in enumerate(norm) if x in {"snp", "marker", "markerid", "id", "rs"}), None)
    chrom = pos = marker = None
    if chr_i is not None and pos_i is not None and max(chr_i, pos_i) < len(parts):
        chrom, pos = parts[chr_i], parts[pos_i]
        marker = parts[marker_i] if marker_i is not None and marker_i < len(parts) else "{}:{}".format(chrom, pos)
    else:
        for token in parts[:4]:
            match = re.match(r"^([^:\s]+):(\d+)(?::.*)?$", token)
            if match:
                chrom, pos, marker = match.group(1), match.group(2), token
                break
        if chrom is None and len(parts) >= 3 and re.fullmatch(r"\d+", parts[1] or ""):
            chrom, pos, marker = parts[0], parts[1], "{}:{}".format(parts[0], parts[1])
    try:
        pos = int(float(pos))
    except (TypeError, ValueError):
        return None
    pvalue = None
    if p_i is not None and p_i < len(parts):
        try:
            pvalue = float(parts[p_i])
        except ValueError:
            pass
    if pvalue is None:
        for token in reversed(parts):
            try:
                value = float(token)
            except ValueError:
                continue
            if 0 <= value <= 1:
                pvalue = value
                break
    if pvalue is None:
        return None
    return {"chrom": str(chrom), "pos": pos, "marker": marker, "pvalue": pvalue}


def phenotype_region(path_text, chrom, start, end, max_files=100, max_records=10000):
    if not path_text:
        return {"status": "not_configured", "files": [], "records": [], "traits": []}
    root = Path(os.path.expandvars(os.path.expanduser(str(path_text)))).resolve()
    if root.is_file():
        paths = [root]
    elif root.is_dir():
        paths = sorted(p for p in root.rglob("*.ps") if p.is_file())[:max_files]
    else:
        raise VCFError("表型 .ps 路径不存在：{}".format(root))
    records, traits = [], []
    for path in paths:
        trait_records = []
        with _open_text(path) as handle:
            header = None
            for line_no, line in enumerate(handle):
                if not line.strip() or line.startswith("#"):
                    continue
                parts = re.split(r"\s+", line.strip())
                if header is None:
                    normalized = {_norm_header(x) for x in parts}
                    if normalized.intersection({"pvalue", "pval", "chrom", "chr", "position", "bp", "marker"}):
                        header = parts
                        continue
                    header = []
                record = _ps_record(parts, header)
                if record and record["chrom"] == str(chrom) and start <= record["pos"] <= end:
                    record["trait"] = path.stem
                    record["source_file"] = str(path)
                    trait_records.append(record)
                    records.append(record)
                    if len(records) >= max_records:
                        break
            if trait_records:
                best = min(trait_records, key=lambda x: x["pvalue"])
                traits.append({"trait": path.stem, "record_count": len(trait_records),
                               "min_p": best["pvalue"], "top_marker": best["marker"],
                               "top_pos": best["pos"], "source_file": str(path)})
        if len(records) >= max_records:
            break
    records.sort(key=lambda x: (x["trait"], x["pos"]))
    return {"status": "ok", "path": str(root), "files": [str(x) for x in paths],
            "records": records, "traits": traits, "truncated": len(records) >= max_records}


class AdvancedAnalyzer:
    def __init__(self, service):
        self.service = service

    def calculate_ld(self, path, chrom, lead_pos, window_bp, threshold, min_samples):
        start, end = max(1, lead_pos - window_bp), lead_pos + window_bp
        metadata, samples, records = self.service.query_region(path, chrom, start, end)
        usable = []
        for record in records:
            if "," in record["alt"]:
                continue
            dosages = [genotype_dosage(record["genotypes"].get(sample)) for sample in samples]
            if sum(x is not None for x in dosages) >= min_samples:
                usable.append((record, dosages))
        lead_candidates = [(r, d) for r, d in usable if r["pos"] == lead_pos]
        if not lead_candidates:
            raise VCFError("Lead 位点不存在，或没有足够的二等位基因型用于 LD 计算")
        lead_record, lead_dosage = lead_candidates[0]
        points = []
        for record, dosage in usable:
            if record is lead_record:
                r2, n = 1.0, sum(x is not None for x in dosage)
            else:
                r2, n = pairwise_r2(lead_dosage, dosage, min_samples)
            if r2 is None:
                continue
            points.append({
                "key": record["key"], "chrom": record["chrom"], "pos": record["pos"],
                "id": record.get("id"), "ref": record["ref"], "alt": record["alt"],
                "variant_type": record["variant_type"], "r2": round(r2, 6), "n": n,
                "distance": record["pos"] - lead_pos,
            })
        points.sort(key=lambda x: x["pos"])
        linked = [x for x in points if x["r2"] >= threshold]
        block_start = min((x["pos"] for x in linked), default=lead_pos)
        block_end = max((x["pos"] for x in linked), default=lead_pos)
        return {
            "lead_record": {k: lead_record.get(k) for k in (
                "key", "chrom", "pos", "id", "ref", "alt", "variant_type", "svtype", "svlen"
            )},
            "search_region": {"chrom": chrom, "start": start, "end": end, "window_bp": window_bp},
            "linkage_region": {"chrom": chrom, "start": block_start, "end": block_end,
                               "span_bp": block_end - block_start + 1},
            "threshold": threshold, "min_samples": min_samples,
            "sample_count": len(samples), "tested_variant_count": len(points),
            "linked_variant_count": len(linked),
            "linked_snp_count": sum(x["variant_type"] == "SNP" for x in linked),
            "linked_variants": linked, "points": points,
        }

    def run_ldblockshow(self, payload, ld_result, phenotype):
        executable_text = str(payload.get("ldblockshow_path") or "").strip()
        use_wsl = bool(payload.get("ldblockshow_wsl"))
        executable = None
        wsl_executable = None
        if use_wsl:
            wsl_executable = shutil.which("wsl.exe") or shutil.which("wsl")
            if not wsl_executable:
                raise VCFError("未找到 WSL；请先启用 Windows Subsystem for Linux，或取消 WSL 模式")
            executable = executable_text or "LDBlockShow"
        else:
            executable = executable_text if executable_text and Path(executable_text).is_file() else shutil.which(executable_text or "LDBlockShow")
            if not executable:
                raise VCFError("未找到 LDBlockShow；请填写可执行文件路径，Windows 可勾选 WSL 模式")
        output_text = str(payload.get("output_dir") or "").strip()
        if not output_text:
            output_text = str(Path(payload["path"]).resolve().parent / "CallVCF_LDBlockShow")
        output_dir = Path(os.path.expandvars(os.path.expanduser(output_text))).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        lead = ld_result["lead_record"]
        prefix_name = re.sub(r"[^A-Za-z0-9_.-]", "_", "lead_{}_{}".format(lead["chrom"], lead["pos"]))
        prefix = output_dir / prefix_name
        region = ld_result["linkage_region"]
        input_vcf = str(Path(payload["path"]).resolve())
        prefix_arg = str(prefix)

        def wsl_path(value):
            proc = subprocess.run([str(wsl_executable), "wslpath", "-a", str(value)],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=30)
            if proc.returncode != 0 or not proc.stdout.strip():
                raise VCFError("WSL 路径转换失败：{}".format(value))
            return proc.stdout.strip()

        if use_wsl:
            input_vcf, prefix_arg = wsl_path(input_vcf), wsl_path(prefix_arg)
        block_cut = "{}:0.90".format(payload.get("r2_threshold", 0.6))
        tool_args = ["-InVCF", input_vcf, "-OutPut", prefix_arg,
                     "-Region", "{}:{}:{}".format(region["chrom"], region["start"], region["end"]),
                     "-SeleVar", "2", "-BlockType", "3", "-BlockCut", block_cut,
                     "-TopSite", "{}:{}".format(lead["chrom"], lead["pos"]), "-OutPng"]
        gff_path = str(payload.get("gff_path") or "").strip()
        if gff_path:
            gff_value = str(_existing_file(gff_path, "GFF3/GTF"))
            tool_args.extend(["-InGFF", wsl_path(gff_value) if use_wsl else gff_value])
        if phenotype and phenotype.get("records"):
            gwas_path = output_dir / (prefix_name + ".gwas.tsv")
            with gwas_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                for item in phenotype["records"]:
                    writer.writerow([item["chrom"], item["pos"], item["pvalue"]])
            tool_args.extend(["-InGWAS", wsl_path(gwas_path) if use_wsl else str(gwas_path)])
        command = ([str(wsl_executable), str(executable)] if use_wsl else [str(executable)]) + tool_args
        try:
            proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  text=True, timeout=int(payload.get("ldblockshow_timeout", 1800)))
        except subprocess.TimeoutExpired:
            raise VCFError("LDBlockShow 运行超过设定时间，已停止")
        outputs = sorted(str(p) for p in output_dir.glob(prefix_name + "*"))
        return {"ok": proc.returncode == 0, "return_code": proc.returncode,
                "command": command, "stdout": proc.stdout[-5000:], "stderr": proc.stderr[-5000:],
                "output_prefix": str(prefix), "output_files": outputs}

    def analyze(self, payload):
        path = payload.get("path")
        loci = parse_loci(payload.get("lead_locus"), limit=1)
        if len(loci) != 1:
            raise VCFError("Lead 分析一次只接受一个位点")
        chrom, lead_pos = loci[0]
        window_bp = int(float(payload.get("window_kb", 500)) * 1000)
        if not 1000 <= window_bp <= 50_000_000:
            raise VCFError("LD 窗口必须在 1 kb 到 50 Mb 之间")
        threshold = float(payload.get("r2_threshold", 0.6))
        if not 0 <= threshold <= 1:
            raise VCFError("r² 阈值必须在 0 到 1 之间")
        min_samples = int(payload.get("min_samples", 20))
        if min_samples < 3:
            raise VCFError("LD 最少有效样本数不能低于 3")
        options = payload.get("options") or {}
        need_ld = bool(options.get("ld") or options.get("phenotype") or options.get("ldblockshow"))
        result = {"lead_locus": "{}:{}".format(chrom, lead_pos), "options": options}

        _, _, _, lead_records = self.service.query_records(path, "{}:{}".format(chrom, lead_pos))
        if not lead_records:
            raise VCFError("VCF 中不存在 Lead 位点 {}:{}".format(chrom, lead_pos))
        lead_record = lead_records[0]
        result["lead_record"] = {k: lead_record.get(k) for k in (
            "key", "chrom", "pos", "end", "id", "ref", "alt", "variant_type", "svtype", "svlen"
        )}

        ld_result = None
        if need_ld:
            ld_result = self.calculate_ld(path, chrom, lead_pos, window_bp, threshold, min_samples)
            result["ld"] = ld_result
        region = ld_result["linkage_region"] if ld_result else {"chrom": chrom, "start": lead_pos, "end": lead_pos}

        gene_result = None
        if options.get("gene") or options.get("domain"):
            gene_result = locate_gene(payload.get("gff_path"), chrom, lead_pos)
            result["gene"] = gene_result
        annotations = []
        if options.get("function") or options.get("domain"):
            annotations = extract_info_annotations(lead_record)
            tsv = annotation_tsv_matches(payload.get("annotation_path"), chrom, lead_pos)
            result["function"] = {"vcf_info": annotations, "external_table": tsv,
                                  "available_info_keys": sorted((lead_record.get("info") or {}).keys())}
        if options.get("domain"):
            identifiers = set()
            for item in annotations:
                identifiers.update([item.get("gene", ""), item.get("gene_id", ""), item.get("feature", "")])
            if gene_result:
                for item in gene_result.get("matches", []):
                    identifiers.update([item.get("name", ""), item.get("id", "")])
            result["domain"] = {"matches": domain_matches(payload.get("domain_path"), identifiers),
                                "identifiers": sorted(x for x in identifiers if x)}

        phenotype = None
        if options.get("phenotype") or options.get("ldblockshow"):
            phenotype = phenotype_region(payload.get("phenotype_path"), chrom, region["start"], region["end"])
            result["phenotype"] = phenotype
        if options.get("ldblockshow"):
            try:
                result["ldblockshow"] = self.run_ldblockshow(payload, ld_result, phenotype)
            except VCFError as exc:
                result["ldblockshow"] = {"ok": False, "error": str(exc), "output_files": []}
        return result
