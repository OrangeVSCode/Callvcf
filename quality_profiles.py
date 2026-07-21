"""Species/ploidy-aware defaults for the CallVCF quality report.

Profiles are deliberately conservative.  They are starting points, not claims
that a single threshold is biologically correct for every cohort.  Every value
can be overridden by the user and the effective source is returned in reports.
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
}


PROFILE_CATALOG = [
    {
        "id": "generic_diploid_animal",
        "name": "通用二倍体动物（异交）",
        "kingdom": "animal",
        "ploidy": 2,
        "genotype_ploidy": "auto",
        "subgenomes": 1,
        "mating_system": "outcrossing",
        "notes": "杂合率采用队列内稳健离群值，不设跨物种绝对上限。",
        "overrides": {},
    },
    {
        "id": "generic_diploid_plant",
        "name": "通用二倍体植物（自然/异交群体）",
        "kingdom": "plant",
        "ploidy": 2,
        "genotype_ploidy": "auto",
        "subgenomes": 1,
        "mating_system": "outcrossing",
        "notes": "HWE仅建议在遗传背景相对一致的群内解释。",
        "overrides": {},
    },
    {
        "id": "plant_inbred",
        "name": "通用自交植物/纯系",
        "kingdom": "plant",
        "ploidy": 2,
        "genotype_ploidy": "auto",
        "subgenomes": 1,
        "mating_system": "selfing",
        "notes": "对残余杂合更敏感；仍应结合材料世代和育种方式。",
        "overrides": {"het_rate_warn": 0.05, "het_rate_critical": 0.10},
    },
    {
        "id": "cotton_population",
        "name": "棉花自然/育种群体（异源四倍体）",
        "kingdom": "plant",
        "ploidy": 4,
        "genotype_ploidy": "auto",
        "subgenomes": 2,
        "mating_system": "mixed",
        "notes": "关闭二倍体HWE硬判定，强调A/D亚基因组、深度与假杂合信号。",
        "overrides": {"hwe_p_warn": None, "ab_lower": 0.15, "ab_upper": 0.85},
    },
    {
        "id": "cotton_inbred",
        "name": "棉花自交系/纯系（异源四倍体）",
        "kingdom": "plant",
        "ploidy": 4,
        "genotype_ploidy": "auto",
        "subgenomes": 2,
        "mating_system": "selfing",
        "notes": "以5%残余杂合作为默认提醒线，并监测高深度假杂合。",
        "overrides": {
            "het_rate_warn": 0.05,
            "het_rate_critical": 0.10,
            "hwe_p_warn": None,
            "ab_lower": 0.15,
            "ab_upper": 0.85,
        },
    },
    {
        "id": "generic_polyploid_plant",
        "name": "通用多倍体植物",
        "kingdom": "plant",
        "ploidy": 4,
        "genotype_ploidy": "auto",
        "subgenomes": 1,
        "mating_system": "unknown",
        "notes": "倍性可修改；不使用二倍体HWE与0.5单峰AB假设作为硬规则。",
        "overrides": {"hwe_p_warn": None, "ab_lower": 0.10, "ab_upper": 0.90},
    },
    {
        "id": "haploid_microbe",
        "name": "通用单倍体微生物/细胞器",
        "kingdom": "microbe",
        "ploidy": 1,
        "genotype_ploidy": 1,
        "subgenomes": 1,
        "mating_system": "clonal",
        "notes": "任何多等位GT都作为混合、污染或倍性设定问题的候选信号。",
        "overrides": {"het_rate_warn": 0.005, "het_rate_critical": 0.02, "hwe_p_warn": None},
    },
    {
        "id": "custom",
        "name": "自定义物种与阈值",
        "kingdom": "other",
        "ploidy": 2,
        "genotype_ploidy": "auto",
        "subgenomes": 1,
        "mating_system": "unknown",
        "notes": "以通用阈值为起点，所有字段均可覆盖。",
        "overrides": {},
    },
]


def profile_catalog():
    return [{k: deepcopy(v) for k, v in item.items() if k != "overrides"} for item in PROFILE_CATALOG]


def _number(value, name, allow_none=False):
    if value in (None, "") and allow_none:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        raise VCFError("阈值 {} 不是有效数字".format(name))


def resolve_profile(config=None):
    config = dict(config or {})
    profile_id = str(config.get("profile_id") or "generic_diploid_plant")
    profile = next((deepcopy(x) for x in PROFILE_CATALOG if x["id"] == profile_id), None)
    if profile is None:
        raise VCFError("未知质量评估 Profile：{}".format(profile_id))

    for key in ("kingdom", "mating_system"):
        if config.get(key):
            profile[key] = str(config[key])
    if config.get("species_name"):
        profile["species_name"] = str(config["species_name"]).strip()
    for key in ("ploidy", "subgenomes"):
        if config.get(key) not in (None, ""):
            try:
                profile[key] = int(config[key])
            except (TypeError, ValueError):
                raise VCFError("{} 必须是整数".format("倍性" if key == "ploidy" else "亚基因组数量"))
            if not 1 <= profile[key] <= 16:
                raise VCFError("{} 必须在 1–16 之间".format("倍性" if key == "ploidy" else "亚基因组数量"))
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

    thresholds = deepcopy(BASE_THRESHOLDS)
    sources = {key: "通用基础阈值" for key in thresholds}
    for key, value in profile.get("overrides", {}).items():
        thresholds[key] = value
        sources[key] = profile["name"]

    custom = config.get("thresholds") or {}
    for key in thresholds:
        if key in custom and custom[key] not in ("",):
            thresholds[key] = _number(custom[key], key, allow_none=True)
            sources[key] = "用户自定义"

    proportion_keys = [
        "sample_missing_warn", "sample_missing_critical", "site_missing_warn",
        "site_missing_critical", "ab_lower", "ab_upper", "maf_structure",
        "hwe_p_warn", "het_rate_warn", "het_rate_critical",
    ]
    for key in proportion_keys:
        value = thresholds.get(key)
        if value is not None and not 0 <= value <= 1:
            raise VCFError("阈值 {} 必须在 0–1 之间".format(key))
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
    profile["thresholds"] = thresholds
    profile["threshold_sources"] = sources
    profile["interpretation"] = {
        "hwe_enabled": thresholds.get("hwe_p_warn") is not None and profile["genotype_ploidy"] in {2, "auto"},
        "diploid_ab_model": profile["genotype_ploidy"] in {2, "auto"},
        "absolute_het_enabled": thresholds.get("het_rate_warn") is not None,
        "cohort_relative_het": thresholds.get("het_rate_warn") is None,
    }
    return profile
