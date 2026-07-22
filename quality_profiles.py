"""Plant-focused species/ploidy defaults for the CallVCF quality report.

Profiles are deliberately conservative.  They are starting points, not claims
that a single threshold is biologically correct for every plant cohort. Animal
analysis is available only through the explicit custom-parameter route.
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
        "id": "plant_haploid_organelle",
        "name": "植物单倍体/细胞器材料",
        "kingdom": "plant",
        "ploidy": 1,
        "genotype_ploidy": 1,
        "subgenomes": 1,
        "mating_system": "clonal",
        "notes": "适用于植物单倍体、叶绿体或线粒体VCF；任何多等位GT都作为混合或倍性设定问题候选信号。",
        "overrides": {"het_rate_warn": 0.005, "het_rate_critical": 0.02, "hwe_p_warn": None},
    },
    {
        "id": "custom",
        "name": "自定义物种与阈值（非植物从这里进入）",
        "kingdom": "other",
        "ploidy": 2,
        "genotype_ploidy": "auto",
        "subgenomes": 1,
        "mating_system": "unknown",
        "notes": "植物可自由定制；动物必须填写核心阈值、GT倍性和繁殖方式并确认，不提供动物默认参数。",
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

    requested_kingdom = str(config.get("kingdom") or profile.get("kingdom") or "other")
    if profile_id != "custom" and requested_kingdom != "plant":
        raise VCFError("CallVCF内置Profile仅用于植物；动物或其他非植物数据必须选择“自定义物种与阈值”")
    profile["kingdom"] = requested_kingdom if profile_id == "custom" else "plant"
    if config.get("mating_system"):
        profile["mating_system"] = str(config["mating_system"])
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
    profile["analysis_scope"] = "non_plant_custom" if profile["kingdom"] != "plant" else ("plant_custom" if profile_id == "custom" else "plant_builtin")
    profile["parameter_responsibility"] = "user" if profile["analysis_scope"] == "non_plant_custom" else "callvcf_plant_profile_or_user_override"
    profile["thresholds"] = thresholds
    profile["threshold_sources"] = sources
    profile["interpretation"] = {
        "hwe_enabled": thresholds.get("hwe_p_warn") is not None and profile["genotype_ploidy"] in {2, "auto"},
        "diploid_ab_model": profile["genotype_ploidy"] in {2, "auto"},
        "absolute_het_enabled": thresholds.get("het_rate_warn") is not None,
        "cohort_relative_het": thresholds.get("het_rate_warn") is None,
    }
    return profile
