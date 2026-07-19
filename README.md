# CallVCF

CallVCF 是一个面向大规模 VCF 的本地交互查询工具。它在浏览器中提供界面，在服务器端调用 `bcftools`，不会把基因型数据上传到第三方服务。

## 功能

- 自动识别 plain VCF、`vcf.gz`/BGZF、BCF、索引状态、样本数、contig 和观察到的 SNP/INDEL/SV 类型。
- 批量检查 `chr:pos` 位点是否存在；同一位置的多条记录分别展示。
- 查询一个或一批位点在所有样本中的基因型分布：`0/0`、`0/1`、`1/1`、缺失及多等位/其他，并列出具体样本。
- 查询指定样本在指定位置的 GT、实际等位基因和类别矩阵。
- 统计指定样本在选定位点或整个 VCF 中的各种基因型数量，并按 SNP、INDEL、SV 分类。
- 对单个或多个 Lead 位点计算窗口内成对 `r²`，按阈值提取连锁变异/SNP 集，并以最左至最右连锁变异跨度定义工作连锁区间。
- 同时绘制独立的 `R²` 与 `|D′|` 三角 LD 热图；未定相二倍体基因型的 `D′` 通过 EM 估计。两张图均可叠加 GFF3/GTF 基因结构轨道并导出高清 PNG 或矢量 PDF。
- 通过 GFF3/GTF 判断位点位于 CDS、外显子、基因内部或基因上下游；读取 VCF `ANN`、`CSQ`、`BCSQ` 功能注释。
- 匹配本地功能注释表和蛋白结构域表，保留具体基因、转录本、HGVS 与结构域记录。
- 读取单个或整目录 EMMAX/GWAS `.ps` 文件，汇总连锁区域内每个表型的最小 P 值与峰值 Marker。
- 可调用本地 LDBlockShow，传入 VCF、连锁区间、GFF 与区域 GWAS 数据，保留 `.blocks.gz`、`.site.gz`、SVG/PNG 等原始结果。
- 页面内置可选工具管理器：一键下载官方 PLINK 1.9 与 LDBlockShow，自动保存到当前用户的 CallVCF 工具目录；第三方工具不上传用户数据。
- 有索引时随机访问；没有索引时自动切换到流式目标筛选。
- 支持服务器目录内 VCF/BCF 文件发现，无需复制大型文件。

## 运行条件

- Python 3.8 或更高版本（仅使用标准库）
- `bcftools` 1.10 或更高版本
- 推荐同时安装 `bgzip` 和 `tabix`

Ubuntu/Debian：

```bash
sudo apt-get install bcftools tabix
```

Conda：

```bash
conda install -c bioconda bcftools
```

## 启动

### Windows 本地页面

安装 Python 3.8+ 后，直接双击 `Start_CallVCF.bat`。工具会自动打开浏览器，并只监听本机 `127.0.0.1`。纯文本 VCF 和 VCF.GZ 可直接读取；安装 `bcftools` 后可加速带索引的大型 VCF，并支持 BCF。

也可以运行：

```powershell
python start_local.py
```

GitHub 的 `build-windows` 工作流可以构建单文件 `CallVCF.exe`；手动触发工作流后从 Actions artifact 下载，发布 Release 时也会自动附加 EXE。

### Linux / 服务器

```bash
git clone https://github.com/OWNER/Callvcf.git
cd Callvcf
python3 app.py --host 127.0.0.1 --port 8765
```

也可以在后台启动：

```bash
chmod +x start.sh stop.sh
./start.sh
./stop.sh
```

浏览器打开 `http://127.0.0.1:8765`。

如果工具运行在远程服务器，建议保持服务只监听本机，再建立 SSH 隧道：

```bash
ssh -L 8765:127.0.0.1:8765 user@server
```

随后在本地浏览器打开 `http://127.0.0.1:8765`。

## 使用方式

1. 输入服务器上的 VCF/BCF 绝对路径，点击“载入并识别”。
2. 位点输入支持 `5:2538289` 或 `5 2538289`，可批量粘贴。
3. 根据需要选择：
   - 位点是否存在
   - 位点基因型分布
   - 样本 × 位点
   - 样本全局统计
   - Lead 位点高级分析
   - 样本 Lead-SNP 画像

样本统计留空位点列表时会扫描整个 VCF。对于数百万记录的文件，这一步可能需要数分钟。

## Lead 位点高级分析

页面中的每项功能都可以独立勾选：

1. “仅指定 Lead”模式填写一个 Lead 位点和左右窗口；“区间内多 Lead 联集”模式填写指定区间及该区间内的 Lead-SNP 列表。
2. 工具计算 Lead 与区间内二等位变异的成对 `r²`。多 Lead 模式返回与任一 Lead 达到阈值的 SNP/INDEL/SV 联集、每个 Lead 的连锁跨度和每个位点对应的 Lead 列表。
3. 页面同步绘制区域 GWAS 信号、GFF3/GTF 基因结构、各 Lead 连锁跨度，以及独立的 `R²` 和 `|D′|` 三角热图；每张图均可输出 2×/4×/6× 高清 PNG 或矢量 PDF，位点过多时按连锁强度抽样。
4. “基因结构轨道”和“相对基因位置”需要 GFF3/GTF；“突变功能”优先读取 VCF INFO 中的 `ANN`/`CSQ`/`BCSQ`，也可补充至少含 `chr`、`pos` 的 TSV/CSV。
5. “蛋白结构域”需要本地 TSV/CSV，表头建议含 `gene`、`gene_id`、`protein` 或 `protein_id`。
6. “区域表型关联”接受一个 `.ps` 文件或包含多个 `.ps` 的目录，可识别常见 `marker/chr/pos/pvalue` 表头和 `chr:pos` Marker。
7. “调用 LDBlockShow”可直接使用页面中的“一键安装”。Windows 上 LDBlockShow 通过 WSL 运行；PLINK 1.9 使用官方 Windows 64 位版。也可手动指定可执行文件。输出默认放在 VCF 同目录的 `CallVCF_LDBlockShow`，也可指定其他目录。
8. LDBlockShow 图形可选择 D′、R² 或两者，并可选择 PLINK Gabriel、Solid Spine、自定义阈值或不划分 block。PLINK Gabriel 模式使用 LDBlockShow 官方包内配套的 PLINK。

注意：页面按成对 `r²` 阈值跨度给出的“工作连锁区间”，与 Gabriel、solid-spine 等正式 haplotype block 算法不是同一概念。勾选 LDBlockShow 后会同时获得软件自身的 block 判定结果。

LDBlockShow 的基本调用形式为：

```bash
LDBlockShow -InVCF input.vcf.gz -OutPut result_prefix -Region chr:start-end -SeleVar 2
```

大型数据建议使用 BGZF 压缩且已建立索引的 VCF。若 Lead 窗口内超过 20,000 条变异，工具会提示缩小窗口，避免浏览器和内存被一次查询占满。

## 样本 Lead-SNP 画像

1. 在页面上方选择一个或多个样本，在画像页粘贴 Lead-SNP 列表，并选择含 `pvalue` 的 GWAS `.ps` 文件或目录。
2. 若结果同时含 `beta` 和效应等位基因（支持 `effect_allele/EA/A1/allele1`），页面会显示样本效应剂量、群体均值、显著程度和效应方向。
3. 可选逐行填写 `表型名=HIGH` 或 `表型名=LOW`。只有指定了有利方向的表型才计算综合趋势指数；指数按样本相对群体的效应剂量、β 方向和 `-log10(P)` 加权，范围为 -100 到 100。
4. 画像会区分“等位基因效应方向”和“样本相对趋势”。例如不利等位基因本身方向为不利，但样本携带量低于群体均值时，样本相对趋势可以是有利。

综合趋势指数用于同一 VCF 群体内的探索性比较，不应替代多环境验证、群体结构校正或育种值估计。

## 可选第三方工具

CallVCF 本体不会把 PLINK 或 LDBlockShow 二进制直接合并进单文件 EXE；用户点击“一键安装”后，软件从官方地址下载到 `%LOCALAPPDATA%\CallVCF\tools`（Linux 为 `~/.local/share/CallVCF/tools`）。这样可以独立更新、保留原始许可证，并避免无谓增大安装包。

- PLINK 1.9：GPL-3.0，官方稳定版 beta 7.11（2025-08-19）。
- LDBlockShow：MIT，官方维护仓库 `hewm2008/LDBlockShow`；官方仅支持 Linux/Unix/macOS，Windows 需要 WSL。

详见 `THIRD_PARTY_NOTICES.md`。

## VCF 建议

随机位点查询推荐使用 BGZF 压缩并建立索引：

```bash
bgzip input.vcf
tabix -p vcf input.vcf.gz
```

未索引文件仍可查询，但工具需要顺序读取文件。

## 安全说明

- 后端使用参数数组调用 `bcftools`，不通过 shell 拼接 VCF 路径或查询内容。
- 默认只监听 `127.0.0.1`，建议通过 SSH 隧道访问。
- 工具只读取用户指定的 VCF/BCF，不修改源文件。

## 测试

```bash
python3 -m unittest discover -s tests -v
```

## License

MIT
