# CallVCF

CallVCF 是一个面向大规模 VCF 的本地交互查询工具。它在浏览器中提供界面，在服务器端调用 `bcftools`，不会把基因型数据上传到第三方服务。

## 功能

- 自动识别 plain VCF、`vcf.gz`/BGZF、BCF、索引状态、样本数、contig 和观察到的 SNP/INDEL/SV 类型。
- 批量检查 `chr:pos` 位点是否存在；同一位置的多条记录分别展示。
- 查询一个或一批位点在所有样本中的基因型分布：`0/0`、`0/1`、`1/1`、缺失及多等位/其他，并列出具体样本。
- 查询指定样本在指定位置的 GT、实际等位基因和类别矩阵。
- 统计指定样本在选定位点或整个 VCF 中的各种基因型数量，并按 SNP、INDEL、SV 分类。
- 对单个 Lead 位点计算窗口内成对 `r²`，按阈值提取连锁变异/SNP 集，并以最左至最右连锁变异跨度定义工作连锁区间。
- 通过 GFF3/GTF 判断位点位于 CDS、外显子、基因内部或基因上下游；读取 VCF `ANN`、`CSQ`、`BCSQ` 功能注释。
- 匹配本地功能注释表和蛋白结构域表，保留具体基因、转录本、HGVS 与结构域记录。
- 读取单个或整目录 EMMAX/GWAS `.ps` 文件，汇总连锁区域内每个表型的最小 P 值与峰值 Marker。
- 可调用本地 LDBlockShow，传入 VCF、连锁区间、GFF 与区域 GWAS 数据，保留 `.blocks.gz`、`.site.gz`、SVG/PNG 等原始结果。
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

样本统计留空位点列表时会扫描整个 VCF。对于数百万记录的文件，这一步可能需要数分钟。

## Lead 位点高级分析

页面中的每项功能都可以独立勾选：

1. 填写一个 Lead 位点、左右窗口、`r²` 阈值和最少共同有效样本数。
2. “LD 区段与连锁变异集”计算 Lead 与窗口内二等位变异的成对 `r²`。页面可视化全部有效点，并可复制位点或下载 TSV。
3. “相对基因位置”需要 GFF3/GTF；“突变功能”优先读取 VCF INFO 中的 `ANN`/`CSQ`/`BCSQ`，也可补充至少含 `chr`、`pos` 的 TSV/CSV。
4. “蛋白结构域”需要本地 TSV/CSV，表头建议含 `gene`、`gene_id`、`protein` 或 `protein_id`。
5. “区域表型关联”接受一个 `.ps` 文件或包含多个 `.ps` 的目录，可识别常见 `marker/chr/pos/pvalue` 表头和 `chr:pos` Marker。
6. “调用 LDBlockShow”需要填写可执行文件路径。官方程序面向 Linux/Unix/macOS；Windows 页面可勾选“通过 WSL 调用”，并填写 WSL 内的命令或路径。输出默认放在 VCF 同目录的 `CallVCF_LDBlockShow`，也可指定其他目录。

注意：页面按成对 `r²` 阈值跨度给出的“工作连锁区间”，与 Gabriel、solid-spine 等正式 haplotype block 算法不是同一概念。勾选 LDBlockShow 后会同时获得软件自身的 block 判定结果。

LDBlockShow 的基本调用形式为：

```bash
LDBlockShow -InVCF input.vcf.gz -OutPut result_prefix -Region chr:start-end -SeleVar 2
```

大型数据建议使用 BGZF 压缩且已建立索引的 VCF。若 Lead 窗口内超过 20,000 条变异，工具会提示缩小窗口，避免浏览器和内存被一次查询占满。

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
