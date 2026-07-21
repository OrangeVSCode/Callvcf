# CallVCF

CallVCF 是一个面向大规模 VCF 的本地交互查询工具。它在浏览器中提供界面，在服务器端调用 `bcftools`，不会把基因型数据上传到第三方服务。

## 功能

- 自动识别 plain VCF、普通 gzip VCF、BGZF VCF、BCF、索引状态、样本数、contig 和观察到的 SNP/INDEL/SV 类型；识别依据为文件内容，不只看扩展名。
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
- `.vcf.gz` 的位点查询、样本统计、LD/热图、基因注释、样本画像、PLINK 与 LDBlockShow 路径均直接读取压缩源文件，不生成同体积的解压 VCF 副本。
- bcftools 可用且 `.tbi/.csi` 索引有效时随机访问；没有索引时自动切换到内存管道或压缩流式筛选，不再落盘生成临时 BCF 子集。
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

安装 Python 3.8+ 后，直接双击 `Start_CallVCF.bat`。工具会自动打开浏览器，并只监听本机 `127.0.0.1`。纯文本 VCF、普通 `.vcf.gz` 和 BGZF `.vcf.gz` 均可直接读取；安装 `bcftools` 后可利用 `.tbi/.csi` 索引加速大型 VCF，并支持 BCF。

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
   - VCF 质量评估与报告

样本统计留空位点列表时会扫描整个 VCF。对于数百万记录的文件，这一步可能需要数分钟。

## VCF 质量评估与报告

“VCF 质量评估”页签按文件、位点、样本、样本对和群体五层生成本地深度质量报告。当前内置通用二倍体动物、通用二倍体植物、自交植物/纯系、通用多倍体植物、棉花群体、棉花自交系和单倍体微生物 Profile，也可修改物种名、生物学倍性、亚基因组数量、繁殖方式及每一项阈值。

可选加载同版本未压缩参考 FASTA 与 `.fai`、样本分组/批次 TSV/CSV、重复区或黑名单 BED/BED.GZ。报告会进行 REF 与 contig 长度核验、Group/Batch 分层比较、区域重叠审计；棉花/多倍体还输出 A/D 亚基因组杂合率和疑似假杂合窗口。VCF 中已有 ANN/CSQ/BCSQ 时会统计注释覆盖与后果类别。未提供这些外部输入时，模块明确显示为 `conditional`，不会伪装成已完成。

样本缺失率×杂合率图使用 P98 主显示区间，并将超范围样本以贴边三角形和样本ID标出；这样既保留极端值证据，也不会让一两个异常样本把其余样本压缩在坐标原点。图中同步显示当前 Profile 的缺失率/杂合率提醒线。

报告除 Header、FILTER、缺失率、MAF、Ti/Tv、DP/GQ/AB/杂合率外，还审计 QUAL/QD/MQ/FS/SOR/RankSum 分布、相位率、MNP 与 SV 子类型、SV 长度/不确定性/支持证据、INDEL 最简表达、相邻重复记录和 1 Mb 变异密度热点。下载包包含逐样本、逐位点（最多20万条评估位点）、SV、密度窗口、推荐过滤和模块可用性 TSV。

生物学倍性与 VCF 的 GT 编码倍性分开设置。例如棉花在生物学上是异源四倍体，但很多 VCF 仍以 `0/0、0/1、1/1` 二倍体形式编码。GT 编码倍性默认自动识别，也可手工固定，避免产生假的倍性不符告警。

扫描模式分为：

- 智能抽样：遍历全部记录，准确统计总数、SNP/INDEL/SV 类型、FILTER 与染色体分布；跨全基因组确定性抽取目标数量位点展开逐样本 GT/DP/GQ/AD 评估。
- 完整扫描：每条记录均进入逐样本指标计算，适合较小 VCF 或需要完整精度的场景。

每次运行会生成独立目录和以下下载文件：`report.html`、`report_summary.json`、`sample_metrics.tsv`、`warnings.tsv`、`run_manifest.json` 和便携的 `CallVCF_QC_report.zip`。HTML 可离线打开，并通过浏览器打印或另存为 PDF。

质量评估可选用 PLINK 1.9 继续计算 HWE、PCA、亲缘关系（IBD/PI_HAT）、LD 衰减和 ROH。群体模块先建立二等位、缺失率和 MAF 过滤后的标记面板；PCA 与亲缘关系使用 LD 剪枝面板，LD 衰减改用未剪枝的 QC 面板，避免把真实短距离相关性预先删掉。PCA 输出稳健距离离群样本，亲缘模块输出疑似样本对、每个样本的最近邻/相似度摘要与连通分量。棉花、自交材料和多倍体 Profile 的 HWE 只作描述，不会因群体结构或繁殖方式导致全样本报警；只有适用 Profile 的 HWE 与高相似样本对会以保守方式进入分层告警，且不会自动删除样本。

质量评估不会修改源 VCF。页面中的安全自动修复执行器使用 bcftools，先生成 15 分钟有效的审计计划，再执行操作。建立 CSI 索引和排序副本需普通确认；标准化/拆分、补全 AC/AN/AF/MAF/NS/F_MISSING/HWE/ExcHet、精确去重、样本子集、按 DP/GQ 等表达式掩蔽低质量 GT、按表达式过滤位点都必须输入本次计划专属的二次确认短语。所有副本先写随机 `.partial` 文件并校验索引，成功后才原子移动到新路径；输入文件或已有目标永不覆盖，执行期间检测到输入变化也会中止。每次副本修复同时保存前后 `bcftools stats` 与 JSON 修复清单，便于复核和追溯。

“智能质控并进行前后复评”把这些步骤组成一个受控流程：依据物种 Profile、SNP/INDEL/SV 构成以及 DP/GQ/QD/MQ/FS/SOR 等字段的实际覆盖情况生成保守默认参数；用户可切换到自定义模式逐项修改。默认只掩蔽明确低质量 GT、补全统计标签并按位点缺失率过滤，不机械删除稀有位点、多等位位点、非 PASS 记录或高缺失样本。执行后自动重新评估原始 VCF 与新副本，报告评分变化、位点/样本数量变化、已解决/新增/持续告警及下一步建议，并输出独立 HTML/JSON 对比报告。样本排除、标准化和其他改变记录集合的动作仍需明确勾选与二次确认。

Windows 版可在“本地软件资源”中一键安装 bcftools。由于官方主要面向 Unix 环境，CallVCF 使用已初始化的 Ubuntu/WSL 作为受控运行后端，并保留安装来源与许可证记录；这比捆绑来源不明的第三方 `bcftools.exe` 更可审计。首次使用前需完成 Ubuntu 用户创建。

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

大型数据建议使用 BGZF 压缩且已建立索引的 VCF。普通 gzip 也能完成所有内置查询，但不能由 tabix 随机定位，因此查询远端位点时仍需从文件开头顺序解压数据流。若 Lead 窗口内超过 20,000 条变异，工具会提示缩小窗口，避免浏览器和内存被一次查询占满。

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

未索引文件仍可查询，工具会直接顺序读取 `.vcf.gz` 数据流，不会先解压成 `.vcf`。如果文件已有 `.tbi/.csi` 但页面显示“检测到索引 · 当前流式读取”，说明当前使用纯 Python 后端；安装/配置 bcftools 后才会显示并使用“已索引 · 随机访问”。

未索引 VCF 无法从文件头直接得知总记录数，因此初次载入时“记录数”可能显示 `—`。点击该卡片中的“统计总数”，CallVCF 会完整顺序读取一次压缩数据流并缓存结果，全程不生成解压副本。以项目测试用的 250 MB、370 样本 SNP BGZF 为例，本机统计 3,002,929 条记录约需 10 秒；实际时间取决于磁盘和压缩率。

## 安全说明

- 后端使用参数数组调用 `bcftools`，不通过 shell 拼接 VCF 路径或查询内容。
- 默认只监听 `127.0.0.1`，建议通过 SSH 隧道访问。
- 工具只读取用户指定的 VCF/BCF，不修改源文件。

## 测试

```bash
python3 -m unittest discover -s tests -v
```

若本机有大型真实 SNP、INDEL、SV 压缩 VCF，可运行只读验收脚本。脚本只把首个 contig 中最多 400 条记录写成临时 `.vcf.gz`，完成后自动删除，不会生成整份解压文件：

```powershell
python tests/real_compressed_smoke.py `
  --snp "D:\path\SNP.vcf.gz" `
  --indel "D:\path\INDEL.vcf.gz" `
  --sv "D:\path\SV.vcf.gz"
```

该验收覆盖内容级 BGZF 识别、样本/contig 表头、SNP/INDEL/SV 分类、位点存在性、基因型分布、样本矩阵、样本统计、LD 和热图数据构建。

质量评估引擎可再运行：

```powershell
python tests/real_quality_smoke.py `
  --snp "D:\path\SNP.vcf.gz" `
  --indel "D:\path\INDEL.vcf.gz" `
  --sv "D:\path\SV.vcf.gz"
```

## License

MIT
