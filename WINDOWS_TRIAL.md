# GPA-Accelerator Windows 试用说明

1. 双击 `GPA-Accelerator.exe`，等待浏览器自动打开本地页面。
2. 在“选择 VCF”中载入 `demo/tiny.vcf`。
3. 打开“位点 × 表型”，选择 `demo/association_phenotype.csv`，位点填 `1:100`，表型填 `SPAD`，匹配方式选择“前缀/家族”，可勾选“同时计算合并表型”，然后开始分析。
4. 实际数据可直接使用 VCF/VCF.GZ/BCF 和 XLSX/CSV/TSV/TXT/PS/PHEN 表型文件；样本 ID 会自动对齐。

## EMMAX

完整 Windows 包已经包含官方 EMMAX `emmax-intel-binary-20120210` 和 `emmax-kin-intel64`，程序首次启动会离线部署，无需再次下载。可在“位点 × 表型”的 EMMAX 卡片查看状态。

官方 EMMAX 是 Ubuntu x86-64 程序，因此 Windows 上真正执行 EMMAX 混合模型前仍需完成一次 Ubuntu/WSL 初始化：

1. Windows 终端运行 `wsl --install`（若已安装可跳过）。
2. 打开 Ubuntu，创建 Linux 用户名和密码并进入一次终端。
3. 重新打开 GPA-Accelerator；状态由“已随包部署”变为“已就绪”。

普通“候选位点 × 表型”快速统计不依赖 EMMAX，可直接使用。全基因组/批量关联的亲缘关系和群体结构校正才需要 EMMAX。

所有数据均在本机处理。关闭页面前可点击右上角“关闭工具”。
