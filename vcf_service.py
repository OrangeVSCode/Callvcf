#!/usr/bin/env python3
import os
import gzip
import re
import shutil
import subprocess
import tempfile
import threading
import time
from collections import Counter, defaultdict
from pathlib import Path


class VCFError(RuntimeError):
    pass


GENOTYPE_LABELS = {
    "HOM_REF": "0/0（纯合参考）",
    "HET": "0/1（杂合）",
    "HOM_ALT": "1/1（纯合变异）",
    "MISSING": "缺失",
    "OTHER": "其他/多等位",
}


def normalize_genotype(gt):
    gt = (gt or ".").split(":", 1)[0].replace("|", "/")
    if gt in {".", "./."} or "." in gt:
        return "MISSING"
    alleles = gt.split("/")
    if len(alleles) != 2:
        return "OTHER"
    if alleles == ["0", "0"]:
        return "HOM_REF"
    if set(alleles) == {"0", "1"}:
        return "HET"
    if alleles == ["1", "1"]:
        return "HOM_ALT"
    return "OTHER"


def genotype_alleles(gt, ref, alt):
    gt = (gt or ".").split(":", 1)[0]
    if "." in gt:
        return "."
    sep = "|" if "|" in gt else "/"
    parts = gt.replace("|", "/").split("/")
    alleles = [ref] + alt.split(",")
    try:
        return sep.join(alleles[int(x)] for x in parts)
    except (ValueError, IndexError):
        return gt


def parse_info(info_text):
    result = {}
    for token in str(info_text or "").split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            result[key] = value
    return result


def classify_variant(ref, alt, bcftools_type="", info_text=""):
    alt_values = alt.split(",")
    type_upper = (bcftools_type or "").upper()
    info = parse_info(info_text)
    if info.get("SVTYPE") not in {None, "", "."}:
        return "SV"
    try:
        if any(abs(int(float(x))) >= 50 for x in info.get("SVLEN", "").split(",") if x not in {"", "."}):
            return "SV"
    except ValueError:
        pass
    if any(a.startswith("<") or "[" in a or "]" in a or a == "*" for a in alt_values):
        return "SV"
    if any(abs(len(a) - len(ref)) >= 50 for a in alt_values):
        return "SV"
    if "BND" in type_upper or "OTHER" in type_upper:
        return "SV" if any("[" in a or "]" in a for a in alt_values) else "OTHER"
    if "INDEL" in type_upper:
        return "INDEL"
    if "SNP" in type_upper:
        return "SNP"
    if all(len(ref) == 1 and len(a) == 1 for a in alt_values):
        return "SNP"
    if any(len(ref) != len(a) for a in alt_values):
        return "INDEL"
    return "OTHER"


def parse_loci(raw, limit=500):
    if isinstance(raw, list):
        text = "\n".join(str(x) for x in raw)
    else:
        text = str(raw or "")
    found = []
    seen = set()
    pattern = re.compile(r"([A-Za-z0-9_.-]+)\s*(?::|\s)\s*([0-9]+)")
    for chrom, pos_text in pattern.findall(text.replace(",", "\n").replace(";", "\n")):
        item = (chrom, int(pos_text))
        if item not in seen:
            seen.add(item)
            found.append(item)
        if len(found) > limit:
            raise VCFError("一次最多查询 {} 个位点".format(limit))
    if not found:
        raise VCFError("未识别到位点；请使用 chr:pos 或 chr pos 格式")
    return found


def detect_compression(path):
    """Return plain, gzip, or bgzf based on file bytes rather than the suffix."""
    path = Path(path)
    with path.open("rb") as handle:
        header = handle.read(12)
        if len(header) < 2 or header[:2] != b"\x1f\x8b":
            return "plain"
        if len(header) < 12 or not (header[3] & 0x04):
            return "gzip"
        extra_length = int.from_bytes(header[10:12], "little")
        extra = handle.read(extra_length)

    offset = 0
    while offset + 4 <= len(extra):
        subfield_id = extra[offset:offset + 2]
        subfield_length = int.from_bytes(extra[offset + 2:offset + 4], "little")
        offset += 4
        if offset + subfield_length > len(extra):
            break
        if subfield_id == b"BC" and subfield_length == 2:
            return "bgzf"
        offset += subfield_length
    return "gzip"


class VCFService:
    def __init__(self, bcftools=None):
        self.bcftools = bcftools or os.environ.get("BCFTOOLS") or shutil.which("bcftools")
        if not self.bcftools:
            raise VCFError("未找到 bcftools")
        self._cache = {}
        self._cache_lock = threading.Lock()

    def _run(self, args, timeout=180, check=True):
        try:
            proc = subprocess.run(
                args,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            raise VCFError("操作超时；请检查 VCF 是否过大、损坏或缺少索引")
        if check and proc.returncode != 0:
            detail = (proc.stderr or proc.stdout or "bcftools 执行失败").strip()
            raise VCFError(detail[-3000:])
        return proc

    def _run_pipeline(self, producer_args, consumer_args, timeout=900):
        """Run bcftools-to-bcftools in memory without writing a large subset file."""
        producer = subprocess.Popen(
            producer_args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        producer_errors = []

        def collect_producer_stderr():
            if producer.stderr:
                producer_errors.append(producer.stderr.read())

        stderr_thread = threading.Thread(target=collect_producer_stderr, daemon=True)
        stderr_thread.start()
        try:
            consumer = subprocess.Popen(
                consumer_args,
                stdin=producer.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if producer.stdout:
                producer.stdout.close()
            try:
                stdout, consumer_stderr = consumer.communicate(timeout=timeout)
            except subprocess.TimeoutExpired:
                consumer.kill()
                producer.kill()
                consumer.communicate()
                raise VCFError("操作超时；请为大型 VCF.GZ 建立 .tbi/.csi 索引后重试")
            producer_code = producer.wait(timeout=30)
            stderr_thread.join(timeout=5)
            producer_stderr = b"".join(producer_errors).decode("utf-8", errors="replace")
            if producer_code != 0 or consumer.returncode != 0:
                detail = (consumer_stderr or producer_stderr or "bcftools 管道执行失败").strip()
                raise VCFError(detail[-3000:])
            return stdout
        finally:
            if producer.poll() is None:
                producer.terminate()

    @staticmethod
    def _validate_file(path_text):
        if not path_text:
            raise VCFError("请输入 VCF 路径")
        path = Path(os.path.expandvars(os.path.expanduser(str(path_text)))).resolve()
        if not path.is_file():
            raise VCFError("文件不存在：{}".format(path))
        return path

    @staticmethod
    def _index_path(path):
        for suffix in (".tbi", ".csi"):
            candidate = Path(str(path) + suffix)
            if candidate.is_file():
                return candidate
        return None

    def inspect(self, path_text, force=False):
        path = self._validate_file(path_text)
        stat = path.stat()
        cache_key = (str(path), stat.st_mtime_ns, stat.st_size)
        if not force:
            with self._cache_lock:
                if cache_key in self._cache:
                    return self._cache[cache_key]

        header = self._run([self.bcftools, "view", "-h", str(path)], timeout=180).stdout
        samples_out = self._run([self.bcftools, "query", "-l", str(path)], timeout=180).stdout
        samples = [x for x in samples_out.splitlines() if x]
        index_path = self._index_path(path)
        record_count = None
        if index_path:
            count_proc = self._run([self.bcftools, "index", "-n", str(path)], timeout=60, check=False)
            if count_proc.returncode == 0 and count_proc.stdout.strip().isdigit():
                record_count = int(count_proc.stdout.strip())

        format_match = re.search(r"##fileformat=([^\r\n]+)", header)
        contigs = re.findall(r"##contig=<ID=([^,>]+)", header)
        observed = self._observe_types(path, limit=500)
        observed_types = Counter(x["variant_type"] for x in observed)

        name_lower = path.name.lower()
        compression = "bcf" if name_lower.endswith(".bcf") else detect_compression(path)
        if name_lower.endswith(".bcf"):
            storage = "BCF"
        elif compression == "bgzf":
            storage = "BGZF VCF（压缩直读）"
        elif compression == "gzip":
            storage = "gzip VCF（压缩直读）"
        elif name_lower.endswith(".vcf") or compression == "plain":
            storage = "plain VCF"
        else:
            storage = "VCF-compatible"

        result = {
            "path": str(path),
            "name": path.name,
            "file_size": stat.st_size,
            "storage": storage,
            "compressed": compression in {"gzip", "bgzf", "bcf"},
            "compression": compression,
            "vcf_version": format_match.group(1) if format_match else "unknown",
            "indexed": bool(index_path),
            "index_usable": bool(index_path),
            "index_path": str(index_path) if index_path else None,
            "query_mode": "indexed random access" if index_path else "streaming target filter",
            "space_mode": "direct source read; no decompressed VCF copy",
            "record_count": record_count,
            "sample_count": len(samples),
            "samples": samples,
            "contig_count": len(contigs),
            "contigs": contigs,
            "observed_types": dict(observed_types),
            "observed_records": len(observed),
            "backend": "bcftools",
        }
        with self._cache_lock:
            self._cache = {cache_key: result}
        return result

    def _observe_types(self, path, limit=500):
        fmt = "%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%INFO\n"
        proc = subprocess.Popen(
            [self.bcftools, "query", "-f", fmt, str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        rows = []
        try:
            for line in proc.stdout:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 6:
                    rows.append({
                        "chrom": parts[0],
                        "pos": int(parts[1]),
                        "ref": parts[2],
                        "alt": parts[3],
                        "variant_type": classify_variant(parts[2], parts[3], parts[4], parts[5]),
                    })
                if len(rows) >= limit:
                    break
        finally:
            if proc.stdout:
                proc.stdout.close()
            proc.terminate()
            try:
                proc.wait(timeout=3)
            except subprocess.TimeoutExpired:
                proc.kill()
        return rows

    def discover(self, root_text, max_depth=5, max_files=500):
        root = Path(os.path.expandvars(os.path.expanduser(str(root_text or ".")))).resolve()
        if not root.is_dir():
            raise VCFError("目录不存在：{}".format(root))
        results = []

        def walk(directory, depth):
            if depth > max_depth or len(results) >= max_files:
                return
            try:
                entries = sorted(os.scandir(str(directory)), key=lambda x: (not x.is_dir(), x.name.lower()))
            except PermissionError:
                return
            for entry in entries:
                if entry.name.startswith("."):
                    continue
                if entry.is_dir(follow_symlinks=False):
                    walk(Path(entry.path), depth + 1)
                elif entry.is_file(follow_symlinks=False):
                    lower = entry.name.lower()
                    if lower.endswith((".vcf", ".vcf.gz", ".vcf.bgz", ".bcf")):
                        stat = entry.stat()
                        results.append({"path": entry.path, "name": entry.name, "size": stat.st_size})
                        if len(results) >= max_files:
                            return

        walk(root, 0)
        return {"root": str(root), "files": results, "truncated": len(results) >= max_files}

    def _validate_samples(self, metadata, selected):
        selected = [str(x).strip() for x in (selected or []) if str(x).strip()]
        seen = set()
        selected = [x for x in selected if not (x in seen or seen.add(x))]
        available = set(metadata["samples"])
        missing = [x for x in selected if x not in available]
        if missing:
            raise VCFError("VCF 中不存在这些样本：{}".format(", ".join(missing[:20])))
        return selected

    def query_records(self, path_text, raw_loci, selected_samples=None):
        metadata = self.inspect(path_text)
        path = Path(metadata["path"])
        loci = parse_loci(raw_loci)
        samples = metadata["samples"] if selected_samples is None else self._validate_samples(metadata, selected_samples)
        include_genotypes = selected_samples is None or bool(samples)
        fmt = "%CHROM\t%POS\t%END\t%ID\t%REF\t%ALT\t%TYPE\t%INFO"
        if include_genotypes:
            fmt += "[\t%GT]"
        fmt += "\n"
        sample_args = ["-s", ",".join(samples)] if samples else []

        if metadata["index_usable"]:
            regions = ",".join("{}:{}".format(chrom, pos) for chrom, pos in loci)
            cmd = [self.bcftools, "query", "-r", regions] + sample_args + ["-f", fmt, str(path)]
            output = self._run(cmd, timeout=300).stdout
        else:
            with tempfile.TemporaryDirectory(prefix="vcf_query_") as tmpdir:
                targets = Path(tmpdir) / "targets.tsv"
                with targets.open("w", encoding="utf-8") as handle:
                    for chrom, pos in loci:
                        handle.write("{}\t{}\t{}\n".format(chrom, pos, pos))
                output = self._run_pipeline(
                    [self.bcftools, "view", "-T", str(targets), "-Ou", str(path)],
                    [self.bcftools, "query"] + sample_args + ["-f", fmt, "-"],
                    timeout=900,
                )

        records = []
        for line in output.splitlines():
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom, pos_text, end_text, rec_id, ref, alt, raw_type, info_text = parts[:8]
            genotypes = parts[8:]
            info = parse_info(info_text)
            record = {
                "chrom": chrom,
                "pos": int(pos_text),
                "end": int(end_text) if end_text.isdigit() else int(pos_text),
                "id": None if rec_id == "." else rec_id,
                "ref": ref,
                "alt": alt,
                "raw_type": raw_type,
                "variant_type": classify_variant(ref, alt, raw_type, info_text),
                "svtype": info.get("SVTYPE"),
                "svlen": info.get("SVLEN"),
                "info": info,
                "info_text": info_text,
                "key": "{}:{}:{}>{}".format(chrom, pos_text, ref, alt),
                "genotypes": dict(zip(samples, genotypes)),
            }
            records.append(record)
        return metadata, loci, samples, records

    def query_region(self, path_text, chrom, start, end, selected_samples=None, max_records=20000):
        """Return records in an interval, including GT for the requested samples."""
        metadata = self.inspect(path_text)
        path = Path(metadata["path"])
        chrom = str(chrom).strip()
        start, end = max(1, int(start)), int(end)
        if not chrom or end < start:
            raise VCFError("无效的查询区间")
        samples = metadata["samples"] if selected_samples is None else self._validate_samples(metadata, selected_samples)
        fmt = "%CHROM\t%POS\t%END\t%ID\t%REF\t%ALT\t%TYPE\t%INFO"
        if samples:
            fmt += "[\t%GT]"
        fmt += "\n"
        sample_args = ["-s", ",".join(samples)] if samples else []
        region = "{}:{}-{}".format(chrom, start, end)
        if metadata["index_usable"]:
            cmd = [self.bcftools, "query", "-r", region] + sample_args + ["-f", fmt, str(path)]
            output = self._run(cmd, timeout=900).stdout
        else:
            output = self._run_pipeline(
                [self.bcftools, "view", "-t", region, "-Ou", str(path)],
                [self.bcftools, "query"] + sample_args + ["-f", fmt, "-"],
                timeout=900,
            )
        records = []
        for line in output.splitlines():
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom_value, pos_text, end_text, rec_id, ref, alt, raw_type, info_text = parts[:8]
            info = parse_info(info_text)
            pos = int(pos_text)
            records.append({
                "chrom": chrom_value,
                "pos": pos,
                "end": int(end_text) if end_text.isdigit() else pos,
                "id": None if rec_id == "." else rec_id,
                "ref": ref,
                "alt": alt,
                "raw_type": raw_type,
                "variant_type": classify_variant(ref, alt, raw_type, info_text),
                "svtype": info.get("SVTYPE"),
                "svlen": info.get("SVLEN"),
                "info": info,
                "info_text": info_text,
                "key": "{}:{}:{}>{}".format(chrom_value, pos, ref, alt),
                "genotypes": dict(zip(samples, parts[8:])),
            })
            if len(records) > int(max_records):
                raise VCFError("区间内变异超过 {} 条；请缩小窗口或提高过滤强度".format(max_records))
        return metadata, samples, records

    def check_loci(self, path_text, raw_loci):
        metadata, loci, _, records = self.query_records(path_text, raw_loci, selected_samples=[])
        by_site = defaultdict(list)
        for record in records:
            by_site[(record["chrom"], record["pos"])].append({k: record[k] for k in (
                "key", "chrom", "pos", "end", "id", "ref", "alt", "variant_type", "raw_type", "svtype", "svlen"
            )})
        results = []
        for chrom, pos in loci:
            matches = by_site.get((chrom, pos), [])
            results.append({"query": "{}:{}".format(chrom, pos), "exists": bool(matches), "records": matches})
        return {"vcf": metadata["path"], "results": results}

    def distributions(self, path_text, raw_loci):
        metadata, loci, samples, records = self.query_records(path_text, raw_loci)
        output = []
        for record in records:
            counts = Counter()
            names = defaultdict(list)
            for sample in samples:
                category = normalize_genotype(record["genotypes"].get(sample, "."))
                counts[category] += 1
                names[category].append(sample)
            output.append({
                "record": {k: record[k] for k in (
                    "key", "chrom", "pos", "end", "id", "ref", "alt", "variant_type", "raw_type", "svtype", "svlen"
                )},
                "total_samples": len(samples),
                "counts": {key: counts.get(key, 0) for key in GENOTYPE_LABELS},
                "samples": {key: names.get(key, []) for key in GENOTYPE_LABELS},
                "labels": GENOTYPE_LABELS,
            })
        found_sites = {(r["chrom"], r["pos"]) for r in records}
        missing_loci = ["{}:{}".format(c, p) for c, p in loci if (c, p) not in found_sites]
        return {"vcf": metadata["path"], "records": output, "missing_loci": missing_loci}

    def sample_locus_matrix(self, path_text, raw_loci, selected_samples):
        metadata = self.inspect(path_text)
        samples = self._validate_samples(metadata, selected_samples)
        if not samples:
            raise VCFError("请至少输入一个样本")
        _, loci, _, records = self.query_records(path_text, raw_loci, selected_samples=samples)
        columns = [{k: r[k] for k in ("key", "chrom", "pos", "ref", "alt", "variant_type", "svtype", "svlen")} for r in records]
        rows = []
        for sample in samples:
            values = []
            for record in records:
                gt = record["genotypes"].get(sample, ".")
                category = normalize_genotype(gt)
                values.append({
                    "record_key": record["key"],
                    "gt": gt,
                    "alleles": genotype_alleles(gt, record["ref"], record["alt"]),
                    "category": category,
                    "label": GENOTYPE_LABELS[category],
                })
            rows.append({"sample": sample, "values": values})
        found_sites = {(r["chrom"], r["pos"]) for r in records}
        missing_loci = ["{}:{}".format(c, p) for c, p in loci if (c, p) not in found_sites]
        return {"vcf": metadata["path"], "columns": columns, "rows": rows, "missing_loci": missing_loci}

    def sample_stats(self, path_text, selected_samples, raw_loci=None):
        metadata = self.inspect(path_text)
        samples = self._validate_samples(metadata, selected_samples)
        if not samples:
            raise VCFError("请至少输入一个样本")
        if len(samples) > 50:
            raise VCFError("全 VCF 统计一次最多选择 50 个样本")

        stats = {sample: {
            "sample": sample,
            "total_records": 0,
            "counts": Counter(),
            "by_variant_type": defaultdict(Counter),
        } for sample in samples}

        if raw_loci:
            _, _, _, records = self.query_records(path_text, raw_loci, selected_samples=samples)
            iterator = ((record["variant_type"], [record["genotypes"].get(s, ".") for s in samples]) for record in records)
        else:
            iterator = self._iter_all_genotypes(metadata["path"], samples)

        processed = 0
        for variant_type, genotypes in iterator:
            processed += 1
            for sample, gt in zip(samples, genotypes):
                category = normalize_genotype(gt)
                item = stats[sample]
                item["total_records"] += 1
                item["counts"][category] += 1
                item["by_variant_type"][variant_type][category] += 1

        results = []
        for sample in samples:
            item = stats[sample]
            results.append({
                "sample": sample,
                "total_records": item["total_records"],
                "counts": {key: item["counts"].get(key, 0) for key in GENOTYPE_LABELS},
                "by_variant_type": {
                    vtype: {key: counter.get(key, 0) for key in GENOTYPE_LABELS}
                    for vtype, counter in sorted(item["by_variant_type"].items())
                },
            })
        return {"vcf": metadata["path"], "processed_records": processed, "labels": GENOTYPE_LABELS, "results": results}

    def _iter_all_genotypes(self, path, samples):
        fmt = "%REF\t%ALT\t%TYPE\t%INFO[\t%GT]\n"
        cmd = [self.bcftools, "query", "-s", ",".join(samples), "-f", fmt, str(path)]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        started = time.time()
        try:
            for line in proc.stdout:
                if time.time() - started > 1800:
                    raise VCFError("全 VCF 样本统计超过 30 分钟，已停止")
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                ref, alt, raw_type, info_text = parts[:4]
                yield classify_variant(ref, alt, raw_type, info_text), parts[4:]
            return_code = proc.wait()
            stderr = proc.stderr.read() if proc.stderr else ""
            if return_code != 0:
                raise VCFError((stderr or "bcftools 查询失败")[-3000:])
        finally:
            if proc.poll() is None:
                proc.terminate()


class PurePythonVCFService(VCFService):
    """Sequential VCF/VCF.GZ reader used when bcftools is unavailable."""

    def __init__(self):
        self.bcftools = None
        self._cache = {}
        self._cache_lock = threading.Lock()

    @staticmethod
    def _open_vcf(path):
        with Path(path).open("rb") as handle:
            magic = handle.read(2)
        if magic == b"\x1f\x8b":
            return gzip.open(str(path), "rt", encoding="utf-8", errors="replace")
        return Path(path).open("r", encoding="utf-8", errors="replace")

    def inspect(self, path_text, force=False):
        path = self._validate_file(path_text)
        if path.name.lower().endswith(".bcf"):
            raise VCFError("本机未找到 bcftools；BCF 需要安装 bcftools，VCF/VCF.GZ 可直接读取")
        stat = path.stat()
        cache_key = (str(path), stat.st_mtime_ns, stat.st_size, "pure")
        if not force:
            with self._cache_lock:
                if cache_key in self._cache:
                    return self._cache[cache_key]

        samples = []
        contigs = []
        version = "unknown"
        observed_types = Counter()
        observed_records = 0
        with self._open_vcf(path) as handle:
            for line in handle:
                if line.startswith("##fileformat="):
                    version = line.strip().split("=", 1)[1]
                elif line.startswith("##contig=<ID="):
                    contigs.append(line.split("##contig=<ID=", 1)[1].split(",", 1)[0].split(">", 1)[0])
                elif line.startswith("#CHROM"):
                    samples = line.rstrip("\r\n").split("\t")[9:]
                elif line.startswith("#"):
                    continue
                else:
                    parts = line.rstrip("\r\n").split("\t")
                    if len(parts) >= 8:
                        observed_types[classify_variant(parts[3], parts[4], "", parts[7])] += 1
                        observed_records += 1
                    if observed_records >= 500:
                        break

        compression = detect_compression(path)
        if compression == "bgzf":
            storage = "BGZF VCF（压缩直读）"
        elif compression == "gzip":
            storage = "gzip VCF（压缩直读）"
        else:
            storage = "plain VCF"
        index_path = self._index_path(path)
        result = {
            "path": str(path),
            "name": path.name,
            "file_size": stat.st_size,
            "storage": storage,
            "compressed": compression in {"gzip", "bgzf"},
            "compression": compression,
            "vcf_version": version,
            "indexed": bool(index_path),
            "index_usable": False,
            "index_path": str(index_path) if index_path else None,
            "query_mode": (
                "direct compressed stream; index detected but bcftools is required to use it"
                if index_path else "direct sequential stream; no decompressed copy"
            ),
            "space_mode": "direct source read; no decompressed VCF copy",
            "record_count": None,
            "sample_count": len(samples),
            "samples": samples,
            "contig_count": len(contigs),
            "contigs": contigs,
            "observed_types": dict(observed_types),
            "observed_records": observed_records,
            "backend": "Pure Python",
        }
        with self._cache_lock:
            self._cache = {cache_key: result}
        return result

    def _iter_vcf_records(self, path, selected_samples, target_set=None, region=None):
        metadata = self.inspect(str(path))
        all_samples = metadata["samples"]
        selected_indices = [all_samples.index(s) for s in selected_samples]
        region_chrom, region_start, region_end = region if region else (None, None, None)
        seen_region_chrom = False
        target_stop = None
        contig_rank = {chrom: index for index, chrom in enumerate(metadata.get("contigs") or [])}
        if target_set and all(chrom in contig_rank for chrom, _ in target_set):
            last_rank = max(contig_rank[chrom] for chrom, _ in target_set)
            last_pos = max(pos for chrom, pos in target_set if contig_rank[chrom] == last_rank)
            target_stop = (last_rank, last_pos)
        with self._open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 8:
                    continue
                chrom = parts[0]
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue
                if target_stop is not None and chrom in contig_rank:
                    rank = contig_rank[chrom]
                    if rank > target_stop[0] or (rank == target_stop[0] and pos > target_stop[1]):
                        break
                if region is not None:
                    if chrom != region_chrom:
                        if seen_region_chrom:
                            break
                        continue
                    seen_region_chrom = True
                    if pos < region_start:
                        continue
                    if pos > region_end:
                        break
                if target_set is not None and (chrom, pos) not in target_set:
                    continue
                ref, alt, info_text = parts[3], parts[4], parts[7]
                info = parse_info(info_text)
                format_keys = parts[8].split(":") if len(parts) > 8 else []
                gt_index = format_keys.index("GT") if "GT" in format_keys else None
                genotypes = {}
                for sample, index in zip(selected_samples, selected_indices):
                    sample_field = parts[9 + index] if 9 + index < len(parts) else "."
                    values = sample_field.split(":")
                    gt = values[gt_index] if gt_index is not None and gt_index < len(values) else "."
                    genotypes[sample] = gt
                end_value = info.get("END", "")
                end = int(end_value) if str(end_value).isdigit() else pos + max(len(ref), 1) - 1
                raw_type = ""
                yield {
                    "chrom": chrom,
                    "pos": pos,
                    "end": end,
                    "id": None if parts[2] == "." else parts[2],
                    "ref": ref,
                    "alt": alt,
                    "raw_type": raw_type,
                    "variant_type": classify_variant(ref, alt, raw_type, info_text),
                    "svtype": info.get("SVTYPE"),
                    "svlen": info.get("SVLEN"),
                    "info": info,
                    "info_text": info_text,
                    "key": "{}:{}:{}>{}".format(chrom, pos, ref, alt),
                    "genotypes": genotypes,
                }

    def query_records(self, path_text, raw_loci, selected_samples=None):
        metadata = self.inspect(path_text)
        path = Path(metadata["path"])
        loci = parse_loci(raw_loci)
        samples = metadata["samples"] if selected_samples is None else self._validate_samples(metadata, selected_samples)
        records = list(self._iter_vcf_records(path, samples, set(loci)))
        return metadata, loci, samples, records

    def query_region(self, path_text, chrom, start, end, selected_samples=None, max_records=20000):
        metadata = self.inspect(path_text)
        path = Path(metadata["path"])
        chrom = str(chrom).strip()
        start, end = max(1, int(start)), int(end)
        if not chrom or end < start:
            raise VCFError("无效的查询区间")
        samples = metadata["samples"] if selected_samples is None else self._validate_samples(metadata, selected_samples)
        records = []
        for record in self._iter_vcf_records(path, samples, None, (chrom, start, end)):
            records.append(record)
            if len(records) > int(max_records):
                raise VCFError("区间内变异超过 {} 条；请缩小窗口或安装 bcftools".format(max_records))
        return metadata, samples, records

    def _iter_all_genotypes(self, path, samples):
        for record in self._iter_vcf_records(Path(path), samples, None):
            yield record["variant_type"], [record["genotypes"].get(sample, ".") for sample in samples]


def create_service(bcftools=None):
    executable = bcftools or os.environ.get("BCFTOOLS") or shutil.which("bcftools")
    if executable:
        service = VCFService(executable)
        service.backend = "bcftools"
        return service
    service = PurePythonVCFService()
    service.backend = "Pure Python"
    return service
