#!/usr/bin/env python3
"""Multi-VCF IBS, IBD and KING-robust evidence integration."""

import bisect
import csv
import html
import json
import math
import os
import shutil
import statistics
import subprocess
import threading
import time
import uuid
import zipfile
from collections import defaultdict
from pathlib import Path

from tool_manager import resolve_plink, resolve_plink2
from vcf_service import VCFError


KINSHIP_THRESHOLDS = [
    (0.354, "疑似重复/同一材料"),
    (0.177, "KING一级亲缘参照区间"),
    (0.0884, "KING二级亲缘参照区间"),
    (0.0442, "KING三级亲缘参照区间"),
]
SOURCE_TYPES = ("snp", "indel", "sv")


def _number(value):
    try:
        result=float(value)
        return result if math.isfinite(result) else None
    except (TypeError,ValueError):
        return None


def _integer(value):
    try: return int(float(value))
    except (TypeError,ValueError): return None


def _pair(sample_1,sample_2):
    left=str(sample_1 or "").strip(); right=str(sample_2 or "").strip()
    if not left or not right or left==right: return None
    return (left,right) if left<right else (right,left)


def _read_space_table(path):
    path=Path(path)
    if not path.is_file() or not path.stat().st_size: return []
    with path.open("r",encoding="utf-8-sig",errors="replace") as handle:
        header=handle.readline().strip().lstrip("#").split()
        return [dict(zip(header,line.split())) for line in handle if line.strip()]


def _parse_kin0(path,source="external_kin0",marker_type="external"):
    output=[]
    for row in _read_space_table(path):
        sample_1=row.get("IID1") or row.get("ID1") or row.get("FID1")
        sample_2=row.get("IID2") or row.get("ID2") or row.get("FID2")
        key=_pair(sample_1,sample_2)
        if not key: continue
        nsnp=_integer(row.get("NSNP")); hethet=_number(row.get("HETHET")); ibs0=_number(row.get("IBS0"))
        het1=_number(row.get("HET1_HOM2")); het2=_number(row.get("HET2_HOM1")); ibs_distance=_number(row.get("IBS"))
        original=(str(sample_1).strip(),str(sample_2).strip())
        if original!=key: het1,het2=het2,het1
        counts=source=="vcf_plink2" or bool(nsnp and any(value is not None and value>1 for value in (hethet,ibs0,het1,het2)))
        def count(value):
            if value is None: return None
            return value if counts or not nsnp else value*nsnp
        ibs0_rate=(ibs0/nsnp if counts and nsnp else ibs0)
        output.append({
            "sample_1":key[0],"sample_2":key[1],"marker_type":marker_type,"source":source,"nsnp":nsnp,
            "hethet_count":count(hethet),"ibs0_count":count(ibs0),"het1_hom2_count":count(het1),"het2_hom1_count":count(het2),
            "ibs0_rate":ibs0_rate,"ibs_distance":ibs_distance/nsnp if counts and nsnp and ibs_distance is not None else ibs_distance,
            "ibs_similarity":None if ibs_distance is None else 1-(ibs_distance/nsnp if counts and nsnp else ibs_distance),
            "king_kinship":_number(row.get("KINSHIP")),"pi_hat":None,"z0":None,"z1":None,"z2":None,
            "raw_counts":counts,
        })
    return output


def _parse_genome(path,marker_type):
    output=[]
    for row in _read_space_table(path):
        key=_pair(row.get("IID1"),row.get("IID2"))
        if not key: continue
        ibs0=_integer(row.get("IBS0")); ibs1=_integer(row.get("IBS1")); ibs2=_integer(row.get("IBS2")); compared=sum(x or 0 for x in (ibs0,ibs1,ibs2))
        similarity=_number(row.get("DST"))
        if similarity is None and compared: similarity=((ibs2 or 0)+.5*(ibs1 or 0))/compared
        output.append({
            "sample_1":key[0],"sample_2":key[1],"marker_type":marker_type,"source":"vcf_plink19","nsnp":compared or None,
            "ibs0_count":ibs0,"ibs1_count":ibs1,"ibs2_count":ibs2,"ibs0_rate":ibs0/compared if compared else None,
            "ibs_similarity":similarity,"genotype_discordance":((ibs0 or 0)+(ibs1 or 0))/compared if compared else None,
            "pi_hat":_number(row.get("PI_HAT")),"z0":_number(row.get("Z0")),"z1":_number(row.get("Z1")),"z2":_number(row.get("Z2")),
            "king_kinship":None,
        })
    return output


def _relationship(kinship,ibs0_rate=None,pi_hat=None):
    if kinship is not None:
        for threshold,label in KINSHIP_THRESHOLDS:
            if kinship>=threshold:
                if threshold==.177 and ibs0_rate is not None:
                    return ("一级亲缘参照：亲子型" if ibs0_rate<.005 else "一级亲缘参照：全同胞型",threshold)
                return label,threshold
        return "未达到KING三级亲缘参照线",0
    if pi_hat is not None:
        if pi_hat>=.90: return "PLINK IBD疑似重复",.90
        if pi_hat>=.375: return "PLINK IBD高亲缘提示",.375
        if pi_hat>=.177: return "PLINK IBD中等亲缘提示",.177
    return "无明确近亲提示",0


def _write_tsv(path,rows,fields):
    rows=list(rows)
    with Path(path).open("w",encoding="utf-8-sig",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",extrasaction="ignore"); writer.writeheader(); writer.writerows(rows)
    return Path(path)


def _write_matrix(path,samples,pairs,value_key,diagonal):
    lookup={_pair(row["sample_1"],row["sample_2"]):row.get(value_key) for row in pairs}
    fields=["sample_id"]+samples
    with Path(path).open("w",encoding="utf-8-sig",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t"); writer.writeheader()
        for left in samples:
            row={"sample_id":left}
            for right in samples: row[right]=diagonal if left==right else lookup.get(_pair(left,right),"")
            writer.writerow(row)
    return Path(path)


def _svg_heatmap(samples,pairs,value_key,title,minimum,maximum):
    selected=samples[:160]; n=len(selected)
    if n<2: return "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 900 220'><text x='30' y='60'>至少需要两个共同样本。</text></svg>"
    lookup={_pair(row["sample_1"],row["sample_2"]):row.get(value_key) for row in pairs}; size=760; cell=size/n; left=95; top=70; width=900; height=860
    parts=["<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {} {}'><rect width='100%' height='100%' fill='white'/><text x='30' y='30' font-size='21' font-weight='700' fill='#173b31'>{}</text>".format(width,height,html.escape(title))]
    for i,left_sample in enumerate(selected):
        if n<=55: parts.append("<text x='{:.1f}' y='62' transform='rotate(-55 {:.1f} 62)' font-size='7' fill='#64776f'>{}</text>".format(left+(i+.5)*cell,left+(i+.5)*cell,html.escape(left_sample[:12])))
        for j,right_sample in enumerate(selected):
            value=(maximum if i==j and value_key=="ibs_similarity" else (.5 if i==j and value_key=="king_kinship" else lookup.get(_pair(left_sample,right_sample))))
            if value is None: color="#eef2f0"
            else:
                q=min(1,max(0,(value-minimum)/max(1e-12,maximum-minimum)))
                if value_key=="king_kinship": color="rgb({},{},{})".format(round(58+190*q),round(121-55*q),round(139-85*q))
                else: color="rgb({},{},{})".format(round(235-185*q),round(244-112*q),round(239-137*q))
            parts.append("<rect x='{:.2f}' y='{:.2f}' width='{:.2f}' height='{:.2f}' fill='{}'><title>{} × {}: {}</title></rect>".format(left+j*cell,top+i*cell,cell+.15,cell+.15,color,html.escape(left_sample),html.escape(right_sample),"NA" if value is None else "{:.5f}".format(value)))
        if n<=55: parts.append("<text x='90' y='{:.1f}' text-anchor='end' font-size='7' fill='#64776f'>{}</text>".format(top+(i+.5)*cell+2,html.escape(left_sample[:12])))
    note="显示全部{}个样本".format(n) if len(samples)<=160 else "样本共{}个；热图展示排序后的前160个，完整矩阵见TSV".format(len(samples))
    parts.append("<text x='30' y='845' font-size='11' fill='#687b73'>{}</text></svg>".format(note)); return "".join(parts)


def _svg_scatter(pairs):
    points=[row for row in pairs if row.get("king_kinship") is not None and row.get("ibs0_rate") is not None]
    if not points: return "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 900 260'><text x='35' y='65'>缺少同时含KING系数与IBS0的样本对。</text></svg>"
    total_points=len(points)
    if total_points>12000:
        priority=[row for row in points if row.get("king_kinship",-9)>=.0442]
        remaining=[row for row in points if row.get("king_kinship",-9)<.0442]
        slots=max(0,12000-len(priority)); step=max(1,math.ceil(len(remaining)/max(1,slots)))
        points=(priority+remaining[::step][:slots])[:max(12000,len(priority))]
    width=900;height=520;left=75;top=45;pw=760;ph=390;xmax=max(.05,max(row["ibs0_rate"] for row in points));ymin=min(-.25,min(row["king_kinship"] for row in points));ymax=max(.52,max(row["king_kinship"] for row in points))
    parts=["<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 900 520'><rect width='100%' height='100%' fill='white'/><text x='28' y='28' font-size='21' font-weight='700' fill='#173b31'>KING-robust kinship × IBS0</text>"]
    for threshold,_ in KINSHIP_THRESHOLDS:
        y=top+ph*(1-(threshold-ymin)/(ymax-ymin)); parts.append("<line x1='{}' y1='{:.1f}' x2='{}' y2='{:.1f}' stroke='#d7a36b' stroke-dasharray='5 4'/><text x='{}' y='{:.1f}' font-size='9' fill='#8a5f33'>{}</text>".format(left,y,left+pw,y,left+pw-2,y-3,threshold))
    for row in points:
        x=left+pw*row["ibs0_rate"]/xmax; y=top+ph*(1-(row["king_kinship"]-ymin)/(ymax-ymin)); color="#b55245" if row["kinship_reference_threshold"]>=.177 else "#2e7a67"
        parts.append("<circle cx='{:.2f}' cy='{:.2f}' r='2.6' fill='{}' opacity='.7'><title>{} × {} · KING {:.4f} · IBS0 {:.4f}</title></circle>".format(x,y,color,html.escape(row["sample_1"]),html.escape(row["sample_2"]),row["king_kinship"],row["ibs0_rate"]))
    note="全部{}个样本对".format(total_points) if total_points==len(points) else "显示{} / {:,}个样本对；保留全部KING≥0.0442候选并对其余点等距抽样，完整数值见TSV".format(len(points),total_points)
    parts.append("<line x1='{0}' y1='{1}' x2='{2}' y2='{1}' stroke='#8ca097'/><line x1='{0}' y1='{3}' x2='{0}' y2='{1}' stroke='#8ca097'/><text x='455' y='502' text-anchor='middle' font-size='12'>IBS0比例</text><text transform='translate(18 250) rotate(-90)' text-anchor='middle' font-size='12'>KING-robust kinship coefficient</text><text x='75' y='518' font-size='10' fill='#687b73'>{4}</text></svg>".format(left,top+ph,left+pw,top,html.escape(note))); return "".join(parts)


class SimilarityJobManager:
    def __init__(self):
        self.plink=resolve_plink(); self.plink2=resolve_plink2(); self._jobs={}; self._lock=threading.Lock()

    @staticmethod
    def default_report_root():
        base=os.environ.get("LOCALAPPDATA") or str(Path.home()/".local"/"share")
        return Path(base)/"GPA-Accelerator"/"similarity-reports"

    def catalog(self):
        self.plink=resolve_plink(); self.plink2=resolve_plink2()
        return {"plink":self.plink,"plink2":self.plink2,"default_report_root":str(self.default_report_root()),"source_types":list(SOURCE_TYPES),"kinship_thresholds":[{"minimum":x,"label":y} for x,y in KINSHIP_THRESHOLDS]}

    @staticmethod
    def _validate_file(value,extensions=None):
        path=Path(os.path.expandvars(os.path.expanduser(str(value or "")))).resolve()
        if not path.is_file(): raise VCFError("输入文件不存在：{}".format(path))
        if extensions and not any(str(path).lower().endswith(ext) for ext in extensions): raise VCFError("文件格式不受支持：{}".format(path.name))
        return path

    def start(self,payload):
        raw_sources=payload.get("sources") if isinstance(payload.get("sources"),dict) else {}
        sources={}
        for marker_type in SOURCE_TYPES:
            value=str(raw_sources.get(marker_type) or "").strip()
            if value: sources[marker_type]=self._validate_file(value,(".vcf",".vcf.gz",".vcf.bgz",".bcf"))
        kin0_text=str(payload.get("kin0_path") or "").strip(); kin0=self._validate_file(kin0_text,(".kin0",".txt",".tsv")) if kin0_text else None
        if not sources and not kin0: raise VCFError("请至少选择一种VCF，或输入一个KING .kin0文件")
        if sources and not self.plink: raise VCFError("根据VCF计算IBS需要PLINK 1.9；请先点一键安装")
        options=payload.get("options") if isinstance(payload.get("options"),dict) else {}
        config={
            "sources":{key:str(value) for key,value in sources.items()},"kin0_path":str(kin0) if kin0 else "",
            "site_missing":min(1,max(0,float(options.get("site_missing") if options.get("site_missing") is not None else .10))),
            "maf":min(.5,max(0,float(options.get("maf") if options.get("maf") is not None else .01))),
            "ld_prune":options.get("ld_prune",True) is not False,"prune_window":max(2,int(options.get("prune_window") or 50)),
            "prune_step":max(1,int(options.get("prune_step") or 5)),"prune_r2":min(1,max(0,float(options.get("prune_r2") or .20))),
            "chromosome_count":max(1,min(200,int(options.get("chromosome_count") or 26))),"nearest_neighbors":max(1,min(50,int(options.get("nearest_neighbors") or 10))),
        }
        root_text=str(payload.get("output_dir") or "").strip(); root=Path(os.path.expandvars(os.path.expanduser(root_text))).resolve() if root_text else self.default_report_root().resolve()
        run_id=time.strftime("%Y%m%d-%H%M%S-")+uuid.uuid4().hex[:8]; run_dir=root/run_id
        job={"id":run_id,"status":"queued","message":"等待运行","progress":0.0,"created_at":time.strftime("%Y-%m-%d %H:%M:%S"),"run_dir":str(run_dir),"artifacts":[],"result":None,"error":None,"cancel":threading.Event()}
        with self._lock: self._jobs[run_id]=job
        thread=threading.Thread(target=self._run,args=(job,config),daemon=True); job["thread"]=thread; thread.start(); return self._public(job)

    def _command(self,args,job,label,log_path):
        command=[str(x) for x in args]; flags=getattr(subprocess,"CREATE_NO_WINDOW",0) if os.name=="nt" else 0
        with log_path.open("a",encoding="utf-8",errors="replace") as log:
            log.write("\n=== {} ===\n$ {}\n".format(label," ".join(command))); log.flush()
            proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,creationflags=flags)
            while proc.poll() is None:
                if job["cancel"].is_set(): proc.kill(); proc.wait(); raise VCFError("样本相似性分析已取消")
                time.sleep(.2)
            log.write("exit={}\n".format(proc.returncode))
        if proc.returncode: raise VCFError("{}失败（退出码{}；详见similarity_analysis.log）".format(label,proc.returncode))

    def _run_source(self,marker_type,vcf,work,run_dir,config,job,log_path,index,total):
        prefix=work/marker_type; panel=prefix/"panel"; prefix.mkdir(parents=True,exist_ok=True)
        input_flag="--bcf" if str(vcf).lower().endswith(".bcf") else "--vcf"
        base=[self.plink,input_flag,vcf,"--chr-set",config["chromosome_count"],"no-x","no-y","no-xy","no-mt","--allow-extra-chr","--double-id","--biallelic-only","strict","--set-missing-var-ids","@:#:$1:$2"]
        if config["site_missing"]<1: base += ["--geno",config["site_missing"]]
        if config["maf"]>0: base += ["--maf",config["maf"]]
        base += ["--make-bed","--out",panel]
        job["message"]="正在建立{}标记面板".format(marker_type.upper()); job["progress"]=round(5+55*(index-.8)/total,2)
        self._command(base,job,"建立{}面板".format(marker_type.upper()),log_path)
        fam=Path(str(panel)+".fam"); bim=Path(str(panel)+".bim")
        samples=sum(1 for _ in fam.open("r",encoding="utf-8",errors="replace")); variants=sum(1 for _ in bim.open("r",encoding="utf-8",errors="replace"))
        if samples<2 or variants<1: raise VCFError("{}过滤后样本或位点不足".format(marker_type.upper()))
        extract=[]; pruned=variants
        if config["ld_prune"] and variants>1:
            prune=prefix/"prune"
            self._command([self.plink,"--bfile",panel,"--chr-set",config["chromosome_count"],"no-x","no-y","no-xy","no-mt","--allow-extra-chr","--indep-pairwise",config["prune_window"],config["prune_step"],config["prune_r2"],"--out",prune],job,"{} LD剪枝".format(marker_type.upper()),log_path)
            prune_file=Path(str(prune)+".prune.in"); pruned=sum(1 for line in prune_file.open("r",encoding="utf-8",errors="replace") if line.strip())
            if pruned: extract=["--extract",prune_file]
        genome=prefix/"genome"
        job["message"]="正在计算{} IBS/IBD".format(marker_type.upper()); job["progress"]=round(5+55*index/total,2)
        self._command([self.plink,"--bfile",panel,"--chr-set",config["chromosome_count"],"no-x","no-y","no-xy","no-mt","--allow-extra-chr",*extract,"--genome","full","--out",genome],job,"{} IBS/IBD".format(marker_type.upper()),log_path)
        ibs_rows=_parse_genome(Path(str(genome)+".genome"),marker_type); king_rows=[]
        if self.plink2:
            king=prefix/"king"
            command=[self.plink2,"--bfile",panel,"--chr-set",config["chromosome_count"],"no-x","no-y","no-xy","no-mt","--allow-extra-chr",*extract,"--make-king-table","counts","cols=fid,id,nsnp,hethet,ibs0,ibs1,ibs,kinship","--out",king]
            try:
                self._command(command,job,"{} KING-robust".format(marker_type.upper()),log_path)
            except VCFError:
                self._command([self.plink2,"--bfile",panel,"--chr-set",config["chromosome_count"],"no-x","no-y","no-xy","no-mt","--allow-extra-chr",*extract,"--make-king-table","counts","--out",king],job,"{} KING-robust（默认列重试）".format(marker_type.upper()),log_path)
            king_rows=_parse_kin0(Path(str(king)+".kin0"),"vcf_plink2",marker_type)
        king_lookup={_pair(row["sample_1"],row["sample_2"]):row for row in king_rows}
        for row in ibs_rows:
            king=king_lookup.get(_pair(row["sample_1"],row["sample_2"]))
            if king:
                for key in ("king_kinship","hethet_count","het1_hom2_count","het2_hom1_count","ibs_distance"):
                    if king.get(key) is not None: row[key]=king[key]
        return ibs_rows,{"marker_type":marker_type,"path":str(vcf),"samples":samples,"variants_after_qc":variants,"variants_used":pruned,"pairs":len(ibs_rows),"ibs_backend":"PLINK 1.9","king_backend":"PLINK 2 KING-robust" if self.plink2 else "not_available"}

    @staticmethod
    def _integrate(source_rows,external_rows):
        grouped=defaultdict(dict)
        for row in source_rows: grouped[_pair(row["sample_1"],row["sample_2"])][row["marker_type"]]=row
        external={_pair(row["sample_1"],row["sample_2"]):row for row in external_rows}
        pairs=[]
        for key in sorted(set(grouped)|set(external)):
            sources=grouped.get(key,{})
            weights=[(row,row.get("nsnp") or 0) for row in sources.values() if row.get("nsnp")]
            total=sum(weight for _,weight in weights)
            def weighted(field):
                current=[(row.get(field),weight) for row,weight in weights if row.get(field) is not None]
                den=sum(weight for _,weight in current); return sum(value*weight for value,weight in current)/den if den else None
            ibs0_count=sum((row.get("ibs0_count") or 0) for row in sources.values()); ibs1_count=sum((row.get("ibs1_count") or 0) for row in sources.values()); ibs2_count=sum((row.get("ibs2_count") or 0) for row in sources.values())
            count_total=ibs0_count+ibs1_count+ibs2_count
            ibs_similarity=((ibs2_count+.5*ibs1_count)/count_total) if count_total else weighted("ibs_similarity")
            ibs0_rate=ibs0_count/count_total if count_total else weighted("ibs0_rate")
            vcf_king=weighted("king_kinship"); ext=external.get(key); king=ext.get("king_kinship") if ext and ext.get("king_kinship") is not None else vcf_king
            if ext and ext.get("ibs0_rate") is not None: external_ibs0=ext["ibs0_rate"]
            else: external_ibs0=None
            if ibs0_rate is None: ibs0_rate=external_ibs0
            pi_hat=weighted("pi_hat"); relation,threshold=_relationship(king,ibs0_rate,pi_hat)
            evidence="external_kin0" if ext and ext.get("king_kinship") is not None else ("vcf_plink2" if vcf_king is not None else "plink19_ibs_ibd")
            concordance="not_comparable"
            if ext and ext.get("king_kinship") is not None and vcf_king is not None: concordance="concordant" if abs(ext["king_kinship"]-vcf_king)<=.05 else "discordant_review"
            item={"sample_1":key[0],"sample_2":key[1],"source_count":len(sources),"marker_types":"+".join(sorted(sources)),"nsnp":total or (ext.get("nsnp") if ext else None),"ibs_similarity":ibs_similarity,"ibs0_rate":ibs0_rate,"pi_hat":pi_hat,"king_kinship":king,"king_evidence":evidence,"vcf_weighted_king":vcf_king,"external_king":ext.get("king_kinship") if ext else None,"evidence_concordance":concordance,"relationship_hint":relation,"kinship_reference_threshold":threshold}
            for marker_type in SOURCE_TYPES:
                row=sources.get(marker_type) or {}; item[marker_type+"_nsnp"]=row.get("nsnp"); item[marker_type+"_ibs_similarity"]=row.get("ibs_similarity"); item[marker_type+"_ibs0_rate"]=row.get("ibs0_rate"); item[marker_type+"_pi_hat"]=row.get("pi_hat"); item[marker_type+"_king"]=row.get("king_kinship")
            pairs.append(item)
        differences=sorted(1-row["ibs_similarity"] for row in pairs if row.get("ibs_similarity") is not None)
        for row in pairs:
            if row.get("ibs_similarity") is None or not differences: row["difference_percentile"]=None
            else:
                difference=1-row["ibs_similarity"]; lo=bisect.bisect_left(differences,difference); hi=bisect.bisect_right(differences,difference); row["difference_percentile"]=round(100*(lo+.5*(hi-lo))/len(differences),3)
        pairs.sort(key=lambda row:(row.get("king_kinship") is not None,row.get("king_kinship") or -9,row.get("ibs_similarity") or -9),reverse=True)
        return pairs

    def _run(self,job,config):
        try:
            self.plink=resolve_plink(); self.plink2=resolve_plink2(); run_dir=Path(job["run_dir"]); run_dir.mkdir(parents=True,exist_ok=False); work=run_dir/".work"; work.mkdir()
            log_path=run_dir/"similarity_analysis.log"; log_path.write_text("GPA-Accelerator sample similarity audit log\n",encoding="utf-8")
            job.update({"status":"running","message":"正在准备样本相似性分析","progress":2.0})
            source_rows=[]; source_summaries=[]; items=list(config["sources"].items())
            for index,(marker_type,path) in enumerate(items,1):
                rows,summary=self._run_source(marker_type,path,work,run_dir,config,job,log_path,index,max(1,len(items))); source_rows.extend(rows); source_summaries.append(summary)
            external_rows=[]
            if config["kin0_path"]:
                job["message"]="正在读取外部KING .kin0"; job["progress"]=65
                external_rows=_parse_kin0(config["kin0_path"])
                source_summaries.append({"marker_type":"external_kin0","path":config["kin0_path"],"samples":len({x["sample_1"] for x in external_rows}|{x["sample_2"] for x in external_rows}),"variants_after_qc":None,"variants_used":max((x.get("nsnp") or 0 for x in external_rows),default=None),"pairs":len(external_rows),"ibs_backend":"IBS0 imported","king_backend":"KING-robust imported"})
            job["message"]="正在整合SNP/INDEL/SV与KING证据"; job["progress"]=72
            integrated=self._integrate(source_rows,external_rows)
            if not integrated: raise VCFError("没有得到可比较的样本对；请检查共同样本和过滤阈值")
            samples=sorted({row["sample_1"] for row in integrated}|{row["sample_2"] for row in integrated})
            adjacency=defaultdict(list)
            for row in integrated: adjacency[row["sample_1"]].append((row["sample_2"],row)); adjacency[row["sample_2"]].append((row["sample_1"],row))
            nearest=[]; sample_summary=[]
            for sample in samples:
                pairs=[row for _,row in adjacency[sample]]; ibs=[row["ibs_similarity"] for row in pairs if row.get("ibs_similarity") is not None]; king=[row["king_kinship"] for row in pairs if row.get("king_kinship") is not None]
                sample_summary.append({"sample_id":sample,"pair_count":len(pairs),"median_ibs_similarity":statistics.median(ibs) if ibs else None,"max_ibs_similarity":max(ibs) if ibs else None,"median_king_kinship":statistics.median(king) if king else None,"max_king_kinship":max(king) if king else None,"king_related_neighbors":sum((row.get("king_kinship") or -9)>=.0442 for row in pairs),"duplicate_candidates":sum((row.get("king_kinship") or -9)>=.354 or (row.get("ibs_similarity") or 0)>=.99 for row in pairs)})
                ranked=sorted(adjacency[sample],key=lambda item:(item[1].get("king_kinship") is not None,item[1].get("king_kinship") or -9,item[1].get("ibs_similarity") or -9),reverse=True)
                for rank,(neighbor,row) in enumerate(ranked[:config["nearest_neighbors"]],1): nearest.append({"sample_id":sample,"rank":rank,"neighbor":neighbor,"marker_types":row["marker_types"],"ibs_similarity":row["ibs_similarity"],"ibs0_rate":row["ibs0_rate"],"pi_hat":row["pi_hat"],"king_kinship":row["king_kinship"],"relationship_hint":row["relationship_hint"],"evidence_concordance":row["evidence_concordance"]})
            suspicious=[row for row in integrated if (row.get("king_kinship") or -9)>=.0442 or (row.get("pi_hat") or -9)>=.177 or (row.get("ibs_similarity") or 0)>=.98]
            pair_fields=["sample_1","sample_2","source_count","marker_types","nsnp","ibs_similarity","ibs0_rate","pi_hat","king_kinship","king_evidence","vcf_weighted_king","external_king","evidence_concordance","relationship_hint","kinship_reference_threshold","difference_percentile"]
            for marker_type in SOURCE_TYPES: pair_fields += [marker_type+"_nsnp",marker_type+"_ibs_similarity",marker_type+"_ibs0_rate",marker_type+"_pi_hat",marker_type+"_king"]
            source_fields=["sample_1","sample_2","marker_type","nsnp","ibs_similarity","ibs0_rate","genotype_discordance","pi_hat","z0","z1","z2","king_kinship","source"]
            source_path=_write_tsv(run_dir/"pairwise_by_marker_type.tsv",source_rows+external_rows,source_fields)
            integrated_path=_write_tsv(run_dir/"pairwise_integrated.tsv",integrated,pair_fields); suspicious_path=_write_tsv(run_dir/"pairs_for_review.tsv",suspicious,pair_fields)
            nearest_path=_write_tsv(run_dir/"sample_nearest_neighbors.tsv",nearest,["sample_id","rank","neighbor","marker_types","ibs_similarity","ibs0_rate","pi_hat","king_kinship","relationship_hint","evidence_concordance"])
            sample_path=_write_tsv(run_dir/"sample_similarity_summary.tsv",sample_summary,["sample_id","pair_count","median_ibs_similarity","max_ibs_similarity","median_king_kinship","max_king_kinship","king_related_neighbors","duplicate_candidates"])
            source_summary_path=_write_tsv(run_dir/"source_summary.tsv",source_summaries,["marker_type","path","samples","variants_after_qc","variants_used","pairs","ibs_backend","king_backend"])
            ibs_matrix=_write_matrix(run_dir/"integrated_ibs_similarity_matrix.tsv",samples,integrated,"ibs_similarity",1)
            king_matrix=_write_matrix(run_dir/"king_robust_kinship_matrix.tsv",samples,integrated,"king_kinship",.5)
            ibs_svg=run_dir/"integrated_ibs_heatmap.svg"; ibs_svg.write_text(_svg_heatmap(samples,integrated,"ibs_similarity","SNP / INDEL / SV 综合 IBS 相似性",.5,1),encoding="utf-8")
            king_svg=run_dir/"king_robust_heatmap.svg"; king_svg.write_text(_svg_heatmap(samples,integrated,"king_kinship","KING-robust kinship coefficient",-.20,.50),encoding="utf-8")
            scatter_svg=run_dir/"king_ibs0_scatter.svg"; scatter_svg.write_text(_svg_scatter(integrated),encoding="utf-8")
            summary={"samples":len(samples),"pairs":len(integrated),"marker_sources":len(items),"external_kin0":bool(external_rows),"pairs_with_ibs":sum(row.get("ibs_similarity") is not None for row in integrated),"pairs_with_king":sum(row.get("king_kinship") is not None for row in integrated),"pairs_at_king_3rd_reference":sum((row.get("king_kinship") or -9)>=.0442 for row in integrated),"duplicate_candidates":sum((row.get("king_kinship") or -9)>=.354 or (row.get("ibs_similarity") or 0)>=.99 for row in integrated),"discordant_external_vcf_king":sum(row["evidence_concordance"]=="discordant_review" for row in integrated)}
            result={"schema_version":"gpa-sample-similarity-1.0","run_id":job["id"],"generated_at":time.strftime("%Y-%m-%d %H:%M:%S"),"configuration":config,"backends":{"plink19":self.plink,"plink2":self.plink2},"summary":summary,"sources":source_summaries,"sample_summary":sample_summary,"preview":integrated[:100],"method_notes":["IBS/PI_HAT按每类VCF单独计算；组合IBS由各类可比较基因型计数精确加权。","有外部.kin0时，外部KING-robust系数作为综合KING主证据；否则使用PLINK 2按各类VCF计算并按NSNP汇总。","KING常用亲缘区间源于人类二倍体研究；植物、自交系和多倍体只作为相对相似性复核线，不直接赋予谱系称谓。","PLINK 1.9 PI_HAT依赖样本等位基因频率假设，不等同于KING-robust。","结果不会自动删除任何材料。"]}
            json_path=run_dir/"similarity_summary.json"; json_path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
            report_path=run_dir/"sample_similarity_report.html"; report_path.write_text(self._render_report(result,integrated,run_dir),encoding="utf-8")
            core=[report_path,json_path,source_summary_path,source_path,integrated_path,suspicious_path,nearest_path,sample_path,ibs_matrix,king_matrix,ibs_svg,king_svg,scatter_svg,log_path]
            zip_path=run_dir/"GPA_Accelerator_sample_similarity.zip"
            with zipfile.ZipFile(zip_path,"w",compression=zipfile.ZIP_DEFLATED) as archive:
                for item in core: archive.write(item,item.name)
            artifacts=[]
            for item in core+[zip_path]:
                base="/api/similarity/artifact?run_id={}&name={}".format(job["id"],item.name); artifacts.append({"name":item.name,"size":item.stat().st_size,"url":base,"view_url":base+"&view=1"})
            job.update({"result":result,"artifacts":artifacts,"status":"complete","message":"样本相似性与亲缘证据报告已生成","progress":100.0})
            shutil.rmtree(work,ignore_errors=True)
        except Exception as exc:
            job["status"]="cancelled" if job["cancel"].is_set() else "failed"; job["error"]=str(exc); job["message"]=str(exc)

    @staticmethod
    def _render_report(result,pairs,run_dir):
        top=pairs[:100]
        rows="".join("<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(html.escape(row["sample_1"]),html.escape(row["sample_2"]),html.escape(row["marker_types"] or "external"),"—" if row["ibs_similarity"] is None else "{:.5f}".format(row["ibs_similarity"]),"—" if row["ibs0_rate"] is None else "{:.5f}".format(row["ibs0_rate"]),"—" if row["pi_hat"] is None else "{:.5f}".format(row["pi_hat"]),"—" if row["king_kinship"] is None else "{:.5f}".format(row["king_kinship"]),html.escape(row["relationship_hint"])) for row in top)
        source_rows="".join("<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(html.escape(str(row["marker_type"])),html.escape(Path(row["path"]).name),row.get("samples") or "—",row.get("variants_used") or "—",row.get("pairs") or 0,html.escape(str(row.get("king_backend") or ""))) for row in result["sources"])
        template="""<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>GPA-Accelerator 样本相似性与亲缘报告</title><style>
body{margin:0;background:#edf4f0;color:#17372e;font-family:system-ui,'Microsoft YaHei',sans-serif}main{max-width:1200px;margin:auto;padding:28px}.hero,.card{margin-bottom:18px;padding:24px;border:1px solid #d5e2dc;border-radius:18px;background:white}.hero{color:white;background:linear-gradient(135deg,#153f33,#397e6a)}h1,h2{margin-top:0}.kpis{display:grid;grid-template-columns:repeat(6,1fr);gap:10px}.kpi{padding:13px;border-radius:12px;background:#edf6f1}.kpi b{display:block;font-size:22px}.chart svg{width:100%;height:auto;max-height:900px}.table{overflow:auto;max-height:620px}table{width:100%;border-collapse:collapse;font-size:12px}th,td{padding:8px;border-bottom:1px solid #e3ebe7;text-align:left;white-space:nowrap}.note{padding:12px;border-radius:10px;background:#fff5e7;border-left:4px solid #d48035;margin:8px 0}.muted{color:#687b73;font-size:12px}@media(max-width:800px){.kpis{grid-template-columns:repeat(2,1fr)}}@media print{body{background:white}main{max-width:none;padding:0}.card,.hero{break-inside:avoid}}
</style></head><body><main><section class='hero'><p>GPA-Accelerator · SAMPLE SIMILARITY &amp; KINSHIP</p><h1>样本相似性、IBS与KING-robust证据整合</h1><p>__GENERATED__</p></section><section class='card'><div class='kpis'><div class='kpi'>样本<b>__SAMPLES__</b></div><div class='kpi'>样本对<b>__PAIRS__</b></div><div class='kpi'>VCF类型<b>__SOURCES__</b></div><div class='kpi'>含IBS<b>__IBS__</b></div><div class='kpi'>含KING<b>__KING__</b></div><div class='kpi'>复核候选<b>__REVIEW__</b></div></div></section>
<section class='card'><h2>输入与计算后端</h2><div class='table'><table><thead><tr><th>来源</th><th>文件</th><th>样本</th><th>使用位点</th><th>样本对</th><th>KING后端</th></tr></thead><tbody>__SOURCE_ROWS__</tbody></table></div></section>
<section class='card'><h2>综合IBS相似性</h2><div class='chart'>__IBS_SVG__</div></section><section class='card'><h2>KING-robust亲缘系数</h2><div class='chart'>__KING_SVG__</div></section><section class='card'><h2>KING × IBS0关系图</h2><div class='chart'>__SCATTER_SVG__</div></section>
<section class='card'><h2>优先复核样本对（前100）</h2><div class='table'><table><thead><tr><th>样本1</th><th>样本2</th><th>标记类型</th><th>IBS相似性</th><th>IBS0</th><th>PI_HAT</th><th>KING</th><th>解释</th></tr></thead><tbody>__PAIR_ROWS__</tbody></table></div></section>
<section class='card'><h2>解释边界</h2><div class='note'>KING的0.354、0.177、0.0884、0.0442阈值来自人类二倍体关系推断。对棉花等自交、多倍体植物，只作为相对相似性和异常样本复核线，不能直接声称“一级/二级亲属”。</div><p class='muted'>KING-robust对群体分层较稳健，但父母来自显著不同亚群时可能低估亲缘；建议与PCA、材料谱系、批次和育种来源共同解释。PLINK PI_HAT与KING不是同一估计量。</p><p class='muted'>方法：Manichaikul et al. (2010), Bioinformatics 26:2867–2873, doi:10.1093/bioinformatics/btq559；PLINK 2 KING-robust文档。</p></section></main></body></html>"""
        values={"__GENERATED__":html.escape(result["generated_at"]),"__SAMPLES__":result["summary"]["samples"],"__PAIRS__":result["summary"]["pairs"],"__SOURCES__":result["summary"]["marker_sources"],"__IBS__":result["summary"]["pairs_with_ibs"],"__KING__":result["summary"]["pairs_with_king"],"__REVIEW__":result["summary"]["pairs_at_king_3rd_reference"],"__SOURCE_ROWS__":source_rows,"__IBS_SVG__":(run_dir/"integrated_ibs_heatmap.svg").read_text(encoding="utf-8"),"__KING_SVG__":(run_dir/"king_robust_heatmap.svg").read_text(encoding="utf-8"),"__SCATTER_SVG__":(run_dir/"king_ibs0_scatter.svg").read_text(encoding="utf-8"),"__PAIR_ROWS__":rows}
        for key,value in values.items(): template=template.replace(key,str(value))
        return template

    def status(self,run_id):
        with self._lock: job=self._jobs.get(str(run_id))
        if not job: raise VCFError("样本相似性任务不存在或服务已重启")
        return self._public(job)

    def cancel(self,run_id):
        with self._lock: job=self._jobs.get(str(run_id))
        if not job: raise VCFError("样本相似性任务不存在")
        job["cancel"].set(); job["message"]="正在取消"; return self._public(job)

    def artifact(self,run_id,name):
        with self._lock: job=self._jobs.get(str(run_id))
        if not job: raise VCFError("样本相似性报告不存在")
        allowed={item["name"] for item in job.get("artifacts",[])}
        if name not in allowed or Path(name).name!=name: raise VCFError("样本相似性报告文件不存在")
        path=(Path(job["run_dir"])/name).resolve()
        if path.parent!=Path(job["run_dir"]).resolve() or not path.is_file(): raise VCFError("样本相似性报告文件不存在")
        return path

    @staticmethod
    def _public(job): return {key:value for key,value in job.items() if key not in {"thread","cancel","run_dir"}}
