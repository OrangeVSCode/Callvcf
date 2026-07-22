"""Plant-focused species, crop and ploidy defaults for GPA-Accelerator quality reports.

The crop presets are transparent starting points.  They never lock a value:
the UI sends any user edits back to :func:`resolve_profile`, which records the
source of every effective field and threshold in the report manifest.
"""

from copy import deepcopy

from vcf_service import VCFError


BASE_THRESHOLDS = {
    "sample_missing_warn": 0.05,
    "sample_missing_critical": 0.10,
    "site_missing_warn": 0.10,
    "site_missing_critical": 0.20,
    "median_dp_min_warn": 5.0,
    "median_dp_max_warn": 100.0,
    "median_gq_warn": 20.0,
    "median_gq_critical": 10.0,
    "ab_lower": 0.25,
    "ab_upper": 0.75,
    "maf_structure": 0.05,
    "hwe_p_warn": 1e-6,
    "het_rate_warn": None,
    "het_rate_critical": None,
    # Optional crop-aware reference lines.  None means descriptive only.
    "site_qual_warn": None,
    "expected_titv_min": None,
}


ANIMAL_REQUIRED_THRESHOLDS = (
    "sample_missing_warn", "sample_missing_critical",
    "site_missing_warn", "site_missing_critical",
    "median_dp_min_warn", "median_dp_max_warn",
    "median_gq_warn", "median_gq_critical",
    "ab_lower", "ab_upper", "maf_structure", "hwe_p_warn",
)


PROFILE_CATALOG = [
    {
        "id": "generic_diploid_plant",
        "name": "通用二倍体植物（自然/异交群体）",
        "kingdom": "plant", "ploidy": 2, "genotype_ploidy": "auto",
        "subgenomes": 1, "mating_system": "outcrossing",
        "notes": "HWE仅建议在遗传背景相对一致的群内解释。",
        "overrides": {},
    },
    {
        "id": "plant_inbred",
        "name": "通用自交植物/纯系",
        "kingdom": "plant", "ploidy": 2, "genotype_ploidy": "auto",
        "subgenomes": 1, "mating_system": "selfing",
        "notes": "对残余杂合更敏感；仍应结合材料世代和育种方式。",
        "overrides": {"het_rate_warn": 0.05, "het_rate_critical": 0.10, "hwe_p_warn": None},
    },
    {
        "id": "cotton_population",
        "name": "棉花自然/育种群体（异源四倍体）",
        "kingdom": "plant", "ploidy": 4, "genotype_ploidy": "auto",
        "subgenomes": 2, "mating_system": "mixed",
        "notes": "关闭二倍体HWE硬判定，强调A/D亚基因组、深度与假杂合信号。",
        "overrides": {"hwe_p_warn": None, "ab_lower": 0.15, "ab_upper": 0.85},
    },
    {
        "id": "cotton_inbred",
        "name": "棉花自交系/纯系（异源四倍体）",
        "kingdom": "plant", "ploidy": 4, "genotype_ploidy": "auto",
        "subgenomes": 2, "mating_system": "selfing",
        "notes": "以5%残余杂合作为默认提醒线，并监测高深度假杂合。",
        "overrides": {
            "het_rate_warn": 0.05, "het_rate_critical": 0.10,
            "hwe_p_warn": None, "ab_lower": 0.15, "ab_upper": 0.85,
        },
    },
    {
        "id": "generic_polyploid_plant",
        "name": "通用多倍体植物",
        "kingdom": "plant", "ploidy": 4, "genotype_ploidy": "auto",
        "subgenomes": 1, "mating_system": "unknown",
        "notes": "倍性可修改；不使用二倍体HWE与0.5单峰AB假设作为硬规则。",
        "overrides": {"hwe_p_warn": None, "ab_lower": 0.10, "ab_upper": 0.90},
    },
    {
        "id": "plant_haploid_organelle",
        "name": "植物单倍体/细胞器材料",
        "kingdom": "plant", "ploidy": 1, "genotype_ploidy": 1,
        "subgenomes": 1, "mating_system": "clonal",
        "notes": "适用于植物单倍体、叶绿体或线粒体VCF。",
        "overrides": {"het_rate_warn": 0.005, "het_rate_critical": 0.02, "hwe_p_warn": None},
    },
    {
        "id": "custom",
        "name": "自定义物种与阈值（非植物从这里进入）",
        "kingdom": "other", "ploidy": 2, "genotype_ploidy": "auto",
        "subgenomes": 1, "mating_system": "unknown",
        "notes": "植物可自由定制；动物必须填写核心阈值、GT倍性和繁殖方式并确认。",
        "overrides": {},
    },
]


# Values below are suggested starting points distilled from 植物细化规划.md.
# They are deliberately not presented as universal biological pass/fail rules.
CROP_CATALOG = [
    {
        "id": "rice", "name": "水稻", "scientific_name": "Oryza sativa",
        "profile_id": "plant_inbred", "reference_hint": "IRGSP-1.0 或项目采用的水稻参考版本",
        "genome_size_mb": 370, "ploidy": 2, "genotype_ploidy": "auto", "subgenomes": 1,
        "mating_system": "selfing",
        "notes": "自交作物；杂合率和参考一致性应重点复核。Ti/Tv 2.3仅作大规模SNP集的经验参考。",
        "threshold_overrides": {"median_dp_min_warn": 10, "site_qual_warn": 30, "expected_titv_min": 2.3},
    },
    {
        # Keep the historical id so saved browser settings remain compatible.
        "id": "wheat", "name": "普通/面包小麦（六倍体）", "scientific_name": "Triticum aestivum",
        "profile_id": "generic_polyploid_plant", "reference_hint": "IWGSC RefSeq（请记录具体版本）",
        "genome_size_mb": 17000, "ploidy": 6, "genotype_ploidy": "auto", "subgenomes": 3,
        "mating_system": "selfing",
        "notes": "六倍体A/B/D亚基因组；重点检查染色体命名、亚基因组映射、深度和等位平衡。",
        "threshold_overrides": {"median_dp_min_warn": 10, "median_dp_max_warn": 200, "site_qual_warn": 30, "ab_lower": 0.10, "ab_upper": 0.90, "hwe_p_warn": None},
    },
    {
        "id": "wheat_durum", "name": "硬粒小麦（四倍体）", "scientific_name": "Triticum turgidum ssp. durum",
        "profile_id": "generic_polyploid_plant", "reference_hint": "Svevo（请记录具体组装版本）",
        "genome_size_mb": 12000, "ploidy": 4, "genotype_ploidy": "auto", "subgenomes": 2,
        "mating_system": "selfing",
        "notes": "四倍体A/B亚基因组；不能沿用普通小麦的D亚基因组规则，需核对染色体命名与参考版本。",
        "threshold_overrides": {"median_dp_min_warn": 10, "median_dp_max_warn": 200, "site_qual_warn": 30, "ab_lower": 0.10, "ab_upper": 0.90, "hwe_p_warn": None},
    },
    {
        "id": "wheat_wild_emmer", "name": "野生二粒小麦（四倍体）", "scientific_name": "Triticum turgidum ssp. dicoccoides",
        "profile_id": "generic_polyploid_plant", "reference_hint": "Zavitan WEW v2.0 或项目采用的野生二粒小麦版本",
        "genome_size_mb": 10100, "ploidy": 4, "genotype_ploidy": "auto", "subgenomes": 2,
        "mating_system": "selfing",
        "notes": "四倍体A/B亚基因组野生材料；天然群体结构和真实杂合可能高于纯系，杂合率阈值应结合材料来源调整。",
        "threshold_overrides": {"median_dp_min_warn": 10, "median_dp_max_warn": 200, "site_qual_warn": 30, "ab_lower": 0.10, "ab_upper": 0.90, "hwe_p_warn": None, "het_rate_warn": 0.10, "het_rate_critical": 0.20},
    },
    {
        "id": "wheat_einkorn", "name": "一粒小麦（二倍体）", "scientific_name": "Triticum monococcum",
        "profile_id": "plant_inbred", "reference_hint": "请选择与材料一致的一粒小麦参考组装并记录版本",
        "genome_size_mb": 5600, "ploidy": 2, "genotype_ploidy": "auto", "subgenomes": 1,
        "mating_system": "selfing",
        "notes": "二倍体A基因组；不要套用四倍体或六倍体小麦的亚基因组与等位平衡解释。",
        "threshold_overrides": {"median_dp_min_warn": 10, "median_dp_max_warn": 150, "site_qual_warn": 30, "hwe_p_warn": None},
    },
    {
        # Keep the historical id so saved browser settings remain compatible.
        "id": "cotton", "name": "陆地棉（异源四倍体）", "scientific_name": "Gossypium hirsutum",
        "profile_id": "cotton_inbred", "reference_hint": "TM-1 或项目采用的陆地棉参考基因组及版本",
        "genome_size_mb": 2500, "ploidy": 4, "genotype_ploidy": "auto", "subgenomes": 2,
        "mating_system": "selfing",
        "notes": "异源四倍体A/D亚基因组（AD1）；关注高深度假杂合和亚基因组偏倚。自然群体可改为棉花群体Profile。",
        "threshold_overrides": {"median_dp_min_warn": 10, "median_dp_max_warn": 150, "site_qual_warn": 30},
    },
    {
        "id": "cotton_barbadense", "name": "海岛棉（异源四倍体）", "scientific_name": "Gossypium barbadense",
        "profile_id": "cotton_inbred", "reference_hint": "Hai7124、Xinhai21 或项目采用的海岛棉参考版本",
        "genome_size_mb": 2500, "ploidy": 4, "genotype_ploidy": "auto", "subgenomes": 2,
        "mating_system": "selfing",
        "notes": "异源四倍体A/D亚基因组（AD2）；不能直接套用陆地棉坐标，需确认参考组装、染色体别名和注释版本一致。",
        "threshold_overrides": {"median_dp_min_warn": 10, "median_dp_max_warn": 150, "site_qual_warn": 30},
    },
    {
        "id": "cotton_herbaceum", "name": "草棉（二倍体）", "scientific_name": "Gossypium herbaceum",
        "profile_id": "plant_inbred", "reference_hint": "A1-Wagad、A1-africanum 或项目采用的草棉参考版本",
        "genome_size_mb": 1700, "ploidy": 2, "genotype_ploidy": "auto", "subgenomes": 1,
        "mating_system": "selfing",
        "notes": "二倍体A1基因组；草棉不是树棉（G. arboreum），也不使用四倍体棉花的A/D亚基因组和假杂合解释。",
        "threshold_overrides": {"median_dp_min_warn": 10, "median_dp_max_warn": 150, "site_qual_warn": 30, "hwe_p_warn": None},
    },
    {
        "id": "soybean", "name": "大豆", "scientific_name": "Glycine max",
        "profile_id": "plant_inbred", "reference_hint": "Williams 82 或项目采用的大豆参考版本",
        "genome_size_mb": 1000, "ploidy": 2, "genotype_ploidy": "auto", "subgenomes": 1,
        "mating_system": "selfing",
        "notes": "典型自交作物；残余杂合提醒线应结合材料世代调整。",
        "threshold_overrides": {"median_dp_min_warn": 10, "site_qual_warn": 30},
    },
    {
        "id": "maize", "name": "玉米", "scientific_name": "Zea mays",
        "profile_id": "generic_diploid_plant", "reference_hint": "B73 RefGen（请记录具体版本）",
        "genome_size_mb": 2300, "ploidy": 2, "genotype_ploidy": "auto", "subgenomes": 1,
        "mating_system": "outcrossing",
        "notes": "异交背景下杂合率较高，不启用固定杂合率上限；群体结构、HWE和LD需分群解释。",
        "threshold_overrides": {"sample_missing_warn": 0.10, "sample_missing_critical": 0.20, "median_dp_min_warn": 10, "site_qual_warn": 30, "het_rate_warn": None, "het_rate_critical": None},
    },
    {
        "id": "millet_sorghum", "name": "谷子/高粱等禾谷类", "scientific_name": "Setaria / Sorghum",
        "profile_id": "plant_inbred", "reference_hint": "请选择与材料一致的谷子、高粱或其他禾谷类参考版本",
        "genome_size_mb": 650, "ploidy": 2, "genotype_ploidy": "auto", "subgenomes": 1,
        "mating_system": "selfing",
        "notes": "该预设覆盖约500–800 Mb的自交或半自交禾谷类；物种和繁殖方式必须按项目修订。",
        "threshold_overrides": {"median_dp_min_warn": 10, "site_qual_warn": 30},
    },
    {
        "id": "arabidopsis", "name": "拟南芥", "scientific_name": "Arabidopsis thaliana",
        "profile_id": "plant_inbred", "reference_hint": "TAIR10 或项目采用的拟南芥参考版本",
        "genome_size_mb": 135, "ploidy": 2, "genotype_ploidy": "auto", "subgenomes": 1,
        "mating_system": "selfing",
        "notes": "小基因组自交模式植物；默认从更严格的缺失率和残余杂合提醒线开始。",
        "threshold_overrides": {"sample_missing_warn": 0.03, "sample_missing_critical": 0.05, "site_missing_warn": 0.05, "site_missing_critical": 0.10, "median_dp_min_warn": 10, "site_qual_warn": 30, "het_rate_warn": 0.02, "het_rate_critical": 0.05},
    },
    {
        "id": "tomato", "name": "番茄", "scientific_name": "Solanum lycopersicum",
        "profile_id": "plant_inbred", "reference_hint": "SL/Heinz 或项目采用的番茄参考版本",
        "genome_size_mb": 950, "ploidy": 2, "genotype_ploidy": "auto", "subgenomes": 1,
        "mating_system": "selfing",
        "notes": "栽培与野生材料差异大；中等杂合提醒线需按群体来源调整。MAF属于下游分析参数，不等同原始VCF质控。",
        "threshold_overrides": {"median_dp_min_warn": 10, "site_qual_warn": 30, "het_rate_warn": 0.10, "het_rate_critical": 0.20},
    },
    {
        "id": "potato", "name": "栽培马铃薯（四倍体）", "scientific_name": "Solanum tuberosum",
        "profile_id": "generic_polyploid_plant", "reference_hint": "DM/PGSC或项目采用的马铃薯参考版本",
        "genome_size_mb": 840, "ploidy": 4, "genotype_ploidy": "auto", "subgenomes": 1,
        "mating_system": "clonal",
        "notes": "按常见四倍体栽培材料起步；二倍体材料须把倍性和GT编码倍性改为2。",
        "threshold_overrides": {"median_dp_min_warn": 10, "median_dp_max_warn": 150, "site_qual_warn": 30, "hwe_p_warn": None},
    },
    {
        "id": "pepper", "name": "辣椒", "scientific_name": "Capsicum annuum",
        "profile_id": "plant_inbred", "reference_hint": "请填写与材料一致的辣椒参考基因组及版本",
        "genome_size_mb": 3200, "ploidy": 2, "genotype_ploidy": "auto", "subgenomes": 1,
        "mating_system": "selfing",
        "notes": "较大二倍体基因组；杂合率阈值应随自交程度、杂交材料或地方品种来源调整。",
        "threshold_overrides": {"median_dp_min_warn": 10, "site_qual_warn": 30, "het_rate_warn": 0.10, "het_rate_critical": 0.20},
    },
    {
        "id": "cucurbit", "name": "瓜类（通用）", "scientific_name": "Cucurbitaceae",
        "profile_id": "generic_diploid_plant", "reference_hint": "请选择黄瓜、西瓜、甜瓜等对应物种的参考版本",
        "genome_size_mb": 400, "ploidy": 2, "genotype_ploidy": "auto", "subgenomes": 1,
        "mating_system": "mixed",
        "notes": "瓜类繁殖系统差异明显；该值仅作二倍体起点，必须按具体物种与材料类型修改。",
        "threshold_overrides": {"median_dp_min_warn": 10, "site_qual_warn": 30, "het_rate_warn": None, "het_rate_critical": None},
    },
    {
        "id": "tobacco", "name": "栽培烟草", "scientific_name": "Nicotiana tabacum",
        "profile_id": "generic_polyploid_plant", "reference_hint": "请填写项目采用的烟草参考基因组及版本",
        "genome_size_mb": 4500, "ploidy": 4, "genotype_ploidy": "auto", "subgenomes": 2,
        "mating_system": "selfing",
        "notes": "异源四倍体；重点检查亚基因组、深度、等位平衡和假杂合，不使用二倍体HWE硬判定。",
        "threshold_overrides": {"median_dp_min_warn": 10, "median_dp_max_warn": 150, "site_qual_warn": 30, "ab_lower": 0.15, "ab_upper": 0.85, "hwe_p_warn": None, "het_rate_warn": 0.10, "het_rate_critical": 0.20},
    },
]


def profile_catalog():
    return [{k: deepcopy(v) for k, v in item.items() if k != "overrides"} for item in PROFILE_CATALOG]


def _profile_by_id(profile_id):
    return next((deepcopy(item) for item in PROFILE_CATALOG if item["id"] == profile_id), None)


def crop_catalog():
    """Return crop presets with their fully resolved, editable threshold values."""
    public = []
    for crop in CROP_CATALOG:
        item = deepcopy(crop)
        profile = _profile_by_id(item["profile_id"])
        thresholds = deepcopy(BASE_THRESHOLDS)
        thresholds.update(profile.get("overrides", {}))
        thresholds.update(item.pop("threshold_overrides", {}))
        item["thresholds"] = thresholds
        public.append(item)
    return public


def _number(value, name, allow_none=False):
    if value in (None, "") and allow_none:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        raise VCFError("阈值 {} 不是有效数字".format(name))


def _same_value(left, right):
    if left in (None, "") and right in (None, ""):
        return True
    try:
        return abs(float(left) - float(right)) < 1e-12
    except (TypeError, ValueError):
        return str(left).strip() == str(right).strip()


def resolve_profile(config=None):
    config = dict(config or {})
    crop_id = str(config.get("crop_id") or "").strip()
    crop = next((deepcopy(item) for item in CROP_CATALOG if item["id"] == crop_id), None)
    if crop_id and crop is None:
        raise VCFError("未知作物预设：{}".format(crop_id))

    profile_id = str(config.get("profile_id") or (crop or {}).get("profile_id") or "generic_diploid_plant")
    if crop and profile_id != crop["profile_id"]:
        raise VCFError("作物预设与Profile不一致；请重新选择作物，或清除作物预设后手动设置")
    profile = _profile_by_id(profile_id)
    if profile is None:
        raise VCFError("未知质量评估 Profile：{}".format(profile_id))

    if crop:
        profile.update({
            "ploidy": crop["ploidy"], "genotype_ploidy": crop["genotype_ploidy"],
            "subgenomes": crop["subgenomes"], "mating_system": crop["mating_system"],
            "species_name": "{}（{}）".format(crop["name"], crop["scientific_name"]),
        })
    requested_kingdom = str(config.get("kingdom") or profile.get("kingdom") or "other")
    if crop and requested_kingdom != "plant":
        raise VCFError("作物预设只能用于植物分析")
    if profile_id != "custom" and requested_kingdom != "plant":
        raise VCFError("GPA-Accelerator内置Profile仅用于植物；动物或其他非植物数据必须选择“自定义物种与阈值”")
    profile["kingdom"] = requested_kingdom if profile_id == "custom" else "plant"

    crop_source = "作物预设：{}".format(crop["name"]) if crop else None
    field_sources = {
        "species_name": crop_source if crop else "用户填写",
        "ploidy": crop_source or profile["name"],
        "genotype_ploidy": crop_source or profile["name"],
        "subgenomes": crop_source or profile["name"],
        "mating_system": crop_source or profile["name"],
    }
    crop_defaults = {
        "species_name": profile.get("species_name"), "ploidy": profile.get("ploidy"),
        "genotype_ploidy": profile.get("genotype_ploidy"), "subgenomes": profile.get("subgenomes"),
        "mating_system": profile.get("mating_system"),
    }

    if config.get("mating_system"):
        profile["mating_system"] = str(config["mating_system"])
        if not _same_value(profile["mating_system"], crop_defaults.get("mating_system")):
            field_sources["mating_system"] = "用户自定义（基于{}）".format(crop["name"]) if crop else "用户自定义"
    if config.get("species_name"):
        profile["species_name"] = str(config["species_name"]).strip()
        if not _same_value(profile["species_name"], crop_defaults.get("species_name")):
            field_sources["species_name"] = "用户自定义（基于{}）".format(crop["name"]) if crop else "用户自定义"
    for key in ("ploidy", "subgenomes"):
        if config.get(key) not in (None, ""):
            try:
                profile[key] = int(config[key])
            except (TypeError, ValueError):
                raise VCFError("{} 必须是整数".format("倍性" if key == "ploidy" else "亚基因组数量"))
            if not 1 <= profile[key] <= 16:
                raise VCFError("{} 必须在 1–16 之间".format("倍性" if key == "ploidy" else "亚基因组数量"))
            if not _same_value(profile[key], crop_defaults.get(key)):
                field_sources[key] = "用户自定义（基于{}）".format(crop["name"]) if crop else "用户自定义"
    gt_ploidy = config.get("genotype_ploidy", profile.get("genotype_ploidy", "auto"))
    if gt_ploidy in (None, "", "auto"):
        profile["genotype_ploidy"] = "auto"
    else:
        try:
            profile["genotype_ploidy"] = int(gt_ploidy)
        except (TypeError, ValueError):
            raise VCFError("VCF GT编码倍性必须是auto或1–16的整数")
        if not 1 <= profile["genotype_ploidy"] <= 16:
            raise VCFError("VCF GT编码倍性必须在1–16之间")
    if not _same_value(profile["genotype_ploidy"], crop_defaults.get("genotype_ploidy")):
        field_sources["genotype_ploidy"] = "用户自定义（基于{}）".format(crop["name"]) if crop else "用户自定义"

    thresholds = deepcopy(BASE_THRESHOLDS)
    sources = {key: "通用基础阈值" for key in thresholds}
    for key, value in profile.get("overrides", {}).items():
        thresholds[key] = value
        sources[key] = profile["name"]
    if crop:
        # Selecting a crop adopts the complete resolved starting set, including
        # values inherited unchanged from the generic profile.
        sources = {key: crop_source for key in thresholds}
        for key, value in crop.get("threshold_overrides", {}).items():
            thresholds[key] = value
    crop_thresholds = deepcopy(thresholds)

    custom = config.get("thresholds") or {}
    if not isinstance(custom, dict):
        raise VCFError("自定义阈值必须是键值对象")
    if profile["kingdom"] == "animal":
        if profile_id != "custom":
            raise VCFError("动物数据只能使用“自定义物种与阈值”Profile")
        if not profile.get("species_name"):
            raise VCFError("动物自定义分析必须填写物种名称")
        if profile.get("genotype_ploidy") == "auto":
            raise VCFError("动物自定义分析必须明确选择VCF GT编码倍性，不能使用自动识别")
        if profile.get("mating_system") in {None, "", "unknown"}:
            raise VCFError("动物自定义分析必须明确选择繁殖/材料类型")
        missing = [key for key in ANIMAL_REQUIRED_THRESHOLDS if custom.get(key) in (None, "")]
        if missing:
            raise VCFError("动物自定义分析还缺少核心阈值：{}".format("、".join(missing)))
        if config.get("animal_parameters_confirmed") is not True:
            raise VCFError("动物自定义分析必须确认：参数由使用者依据物种、群体设计和测序方案自行设定")
    for key in thresholds:
        if key in custom and custom[key] not in ("",):
            thresholds[key] = _number(custom[key], key, allow_none=True)
            if not crop or not _same_value(thresholds[key], crop_thresholds.get(key)):
                sources[key] = "用户自定义（基于{}）".format(crop["name"]) if crop else "用户自定义"

    proportion_keys = [
        "sample_missing_warn", "sample_missing_critical", "site_missing_warn",
        "site_missing_critical", "ab_lower", "ab_upper", "maf_structure",
        "hwe_p_warn", "het_rate_warn", "het_rate_critical",
    ]
    for key in proportion_keys:
        value = thresholds.get(key)
        if value is not None and not 0 <= value <= 1:
            raise VCFError("阈值 {} 必须在 0–1 之间".format(key))
    for key in ("median_dp_min_warn", "median_dp_max_warn", "median_gq_warn", "median_gq_critical", "site_qual_warn", "expected_titv_min"):
        value = thresholds.get(key)
        if value is not None and value < 0:
            raise VCFError("阈值 {} 不能小于0".format(key))
    if thresholds["sample_missing_warn"] > thresholds["sample_missing_critical"]:
        raise VCFError("样本缺失率 warning 不能高于 critical")
    if thresholds["site_missing_warn"] > thresholds["site_missing_critical"]:
        raise VCFError("位点缺失率 warning 不能高于 critical")
    if thresholds["median_gq_critical"] > thresholds["median_gq_warn"]:
        raise VCFError("GQ critical 阈值不能高于 warning 阈值")
    if thresholds["ab_lower"] >= thresholds["ab_upper"]:
        raise VCFError("AB 下限必须小于上限")
    if thresholds.get("het_rate_warn") is not None and thresholds.get("het_rate_critical") is not None:
        if thresholds["het_rate_warn"] > thresholds["het_rate_critical"]:
            raise VCFError("杂合率 warning 不能高于 critical")

    profile.pop("overrides", None)
    if crop:
        profile.update({
            "crop_id": crop["id"], "crop_name": crop["name"],
            "crop_scientific_name": crop["scientific_name"],
            "reference_hint": crop["reference_hint"], "genome_size_mb": crop["genome_size_mb"],
            "crop_notes": crop["notes"], "preset_disclaimer": "建议起始值；已按当前界面参数解析，可自由修改且需结合项目设计复核。",
        })
    else:
        profile.update({"crop_id": None, "crop_name": None, "reference_hint": None, "genome_size_mb": None})
    profile["analysis_scope"] = "non_plant_custom" if profile["kingdom"] != "plant" else ("plant_crop_preset" if crop else ("plant_custom" if profile_id == "custom" else "plant_builtin"))
    profile["parameter_responsibility"] = "user" if profile["analysis_scope"] == "non_plant_custom" else "callvcf_plant_preset_or_user_override"
    profile["field_sources"] = field_sources
    profile["thresholds"] = thresholds
    profile["threshold_sources"] = sources
    profile["interpretation"] = {
        "hwe_enabled": thresholds.get("hwe_p_warn") is not None and profile["genotype_ploidy"] in {2, "auto"},
        "diploid_ab_model": profile["genotype_ploidy"] in {2, "auto"},
        "absolute_het_enabled": thresholds.get("het_rate_warn") is not None,
        "cohort_relative_het": thresholds.get("het_rate_warn") is None,
        "titv_reference_enabled": thresholds.get("expected_titv_min") is not None,
        "qual_reference_enabled": thresholds.get("site_qual_warn") is not None,
    }
    return profile
