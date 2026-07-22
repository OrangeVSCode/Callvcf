# GPA-Accelerator Windows 试用说明

1. 双击 `GPA-Accelerator.exe`，等待浏览器自动打开本地页面。
2. 在“选择 VCF”中载入 `demo/tiny.vcf`。
3. 打开“位点 × 表型”，选择 `demo/association_phenotype.csv`，位点填 `1:100`，表型填 `SPAD`，匹配方式选择“前缀/家族”，可勾选“同时计算合并表型”，然后开始分析。
4. 打开“表型质控与统计”，选择 `demo/phenotype_demo.csv`；默认勾选遗传力，可查看单次观测、环境内均值和多年多点均值H²、G×E、分环境结果及置信度。
5. 打开“样本相似性与亲缘”，先只选择 `demo/demo.kin0` 即可快速试用外部KING表导入；也可把 `demo/tiny.vcf` 放入SNP框测试VCF计算。首次按VCF计算KING时，在“本地注释与软件资源”中分别安装PLINK 1.9和PLINK 2。
6. 实际数据可直接使用 VCF/VCF.GZ/BCF 和 XLSX/CSV/TSV/TXT/PS/PHEN 表型文件；样本 ID 会自动对齐。

## 样本相似性与亲缘

- SNP、INDEL、SV三个输入框可只填一个，也可任意填写两个或三个；每类独立计算后再生成综合结果。
- `.kin0` 可单独导入，也可与VCF同时使用；若文件只有`IBS0`而没有`IBS`列，软件不会虚构综合IBS相似性。
- 输出包括逐类型样本对、综合样本对、最近邻、复核候选、IBS/KING方阵、热图、关系图、HTML报告与ZIP。
- 页面上的KING关系线来自人类二倍体资料，对棉花等植物只用于异常材料复核，不自动删除材料。

## EMMAX

完整 Windows 包已经包含官方 EMMAX `emmax-intel-binary-20120210` 和 `emmax-kin-intel64`，程序首次启动会离线部署，无需再次下载。可在“位点 × 表型”的 EMMAX 卡片查看状态。

官方 EMMAX 是 Ubuntu x86-64 程序，因此 Windows 上真正执行 EMMAX 混合模型前仍需完成一次 Ubuntu/WSL 初始化：

1. Windows 终端运行 `wsl --install`（若已安装可跳过）。
2. 打开 Ubuntu，创建 Linux 用户名和密码并进入一次终端。
3. 重新打开 GPA-Accelerator；状态由“已随包部署”变为“已就绪”。

普通“候选位点 × 表型”快速统计不依赖 EMMAX，可直接使用。全基因组/批量关联的亲缘关系和群体结构校正才需要 EMMAX。

完成 Ubuntu 初始化且 PLINK 状态为“已就绪”后，在“位点 × 表型”下展开“全基因组 EMMAX”：

1. 复用上方表型文件、表型名/前缀和输出目录。
2. 默认选择推荐的 BN kinship；需要时可指定已有 `.kinf` 和协变量文件。
3. 点击“运行全基因组 EMMAX”，页面会显示后台阶段与进度。
4. 完成后可直接查看 Manhattan、Q-Q、λGC、伪遗传力和 Top hits，并下载每个表型的 `emmax_results.tsv.gz` 及报告 ZIP。

默认不保留大型 TPED/TFAM/BED 中间文件，以节省磁盘空间；如有审计需要再勾选保留。

所有数据均在本机处理。关闭页面前可点击右上角“关闭工具”。
