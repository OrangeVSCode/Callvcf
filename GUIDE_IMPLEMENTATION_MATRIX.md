# 四份质量评估指南实施矩阵

本文件对应《总规划》《vcf认识》《补充规划1》《补充规划2》。状态不是营销标签，而是报告运行时口径：`complete` 表示当前输入可计算，`conditional` 表示必须提供额外数据，`not_applicable` 表示不能由 VCF 单独可靠推断。

## 已实现

- Header、reference/source/FILTER/INFO/FORMAT/contig 审计，排序、畸形记录、相邻重复和 INDEL 最简表达检查。
- SNP、INDEL、MNP、SVTYPE/BND/CNV 分类，FILTER、Ti/Tv、MAF、位点缺失率和 1 Mb 密度热点。
- QUAL、QD、MQ、FS、SOR、MQRankSum、ReadPosRankSum 的覆盖率和分位数。
- 样本 call rate、缺失率、DP、GQ、AB、杂合率、GT 倍性和相位率。
- 基于规范化双等位、缺失率/MAF过滤且LD剪枝面板的 HWE、PCA、IBD/PI_HAT、IBS、基因型不一致率、差异分、最近邻、近交系数F、LD衰减和ROH。
- SV 长度、IMPRECISE、CIPOS/CIEND 和支持证据字段覆盖率。
- HTML、JSON、TSV、ZIP、运行清单、告警表、推荐过滤表；HTML可打印/另存PDF。
- bcftools 两阶段安全执行器：索引、排序副本、标准化/拆分副本、统计标签副本、表达式过滤副本；永不覆盖源文件或已有目标。

## 有额外输入时可用

- GFF3/GTF 基因位置、基因结构轨道；ANN/CSQ/BCSQ 或外部表的功能后果；蛋白结构域；GWAS `.ps` 与表型联动：在 Lead 位点高级分析页使用。
- LDBlockShow 正式 block、R²/D′ 区域图：安装后在 Lead 位点高级分析页使用。
- FASTA 的逐位点 REF 一致性与 build 验证、样本表驱动的品种内/批次分析、BED 重复区/黑名单重叠：报告目前明确标为 `conditional`，下一阶段接入专用输入框与规则。

## 不由 VCF 单独下结论

- 污染比例、测序/比对质量、CNV read-depth 证据需要 BAM/CRAM 或专用 sidecar。
- 临床致病性、外部人群频率和疾病数据库仅适用于配置了相应物种数据库的插件。
- UMAP 只适合辅助浏览，不替代 PCA/亲缘估计，当前不作为质量判定。

## 解释限制

- 棉花、自交系和多倍体的 HWE 默认仅描述；杂合率优先使用队列稳健离群和材料类型解释。
- PLINK 1.9 输出的是 IBD/PI_HAT 与 IBS DST，不冒充 KING kinship。高相似样本对是复核证据，不自动删除样本。
- 智能扫描时总记录、类型、FILTER、排序和密度使用全文件；基因型和高成本分布使用跨全基因组确定性抽样。
