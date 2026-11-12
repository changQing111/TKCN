# TKCN

[English README](README.md)

TKCN（Targeted K-mer Copy Number Profiling，靶向 k-mer 拷贝数分析）使用目标拷贝之间共有、且在基因组目标区域中特异的 k-mer，估算重复基因或多拷贝基因家族的拷贝数。本文以小麦 `TaFT-B1` 为例，介绍完整流程：候选拷贝鉴定 → 序列提取 → 内存求共有 k-mer 交集 → Bowtie1 特异性筛选 → panel 质控 → 重测序 reads 计数 → 拷贝数估算。

## 方法概述

对于目标拷贝组 `T`，定义：

```text
K(copy) = 一个 copy FASTA 中 canonical ACGT k-mer 的集合
target_shared = 目标组中所有拷贝 K(copy) 的交集
S(target,j) = 在参考 j 中至少命中一个目标区域，且无任何非目标命中的 shared k-mer
target_specific = 所有参考基因组 S(target,j) 的交集
```

`shared_kmers_memory.py` 直接从每个 copy FASTA 枚举 canonical k-mer（序列及其反向互补序列中字典序较小者），在单个 copy 内去重后逐个在内存中取交集；包含非 ACGT 碱基的窗口会被跳过，无需为每个 copy 构建 KMC 数据库。

`bowtie_shared_specific.py` 将多个目标的 shared k-mer 自动合并成带目标标识的联合 FASTA，对每个未屏蔽的完整参考基因组各执行一次 Bowtie1 `-f -S -v 0 -a` 比对，检查正反链的全部精确命中。在每个参考中，保留的 k-mer 必须至少命中一个目标区域，且没有任何非目标区域命中；每次命中的完整区间必须落在同一目标的某一个区域内。最后对各参考的保留集合取交集，得到 shared-specific panel。该步骤替代原来的基因组屏蔽、背景 k-mer 并集及基因组差集运算，参考基因组的 Bowtie1 索引可重复使用。

## 依赖软件

需要安装以下软件：

| 软件 | 用途 |
| --- | --- |
| `miniprot` | 使用蛋白质序列搜索候选基因拷贝 |
| `samtools` | 按基因组区域提取 FASTA 序列 |
| `Bowtie1` (`bowtie`, `bowtie-build`) | 全命中精确比对和参考索引 |
| `KMC3`（`kmc`、`kmc_tools`） | 统计重测序 reads 并提取 panel 计数；copy FASTA 交集及 Bowtie1 特异性筛选无需 KMC |
| `bash`、GNU `awk`、`python3`、`SciPy` | 运行辅助脚本、NB 模型和 panel 过滤；内存交集与 Bowtie1 脚本仅使用 Python 标准库 |

所有脚本均支持 `-h/--help`。如果程序不在 `PATH` 中，可通过环境变量指定其位置：

```bash
export MINIPROT_BIN=/path/to/miniprot
export SAMTOOLS_BIN=/path/to/samtools
export BOWTIE_BIN=/home/changq/my_data/bio_soft/bowtie-1.3.1-linux-x86_64/bowtie
export BOWTIE_BUILD_BIN=/home/changq/my_data/bio_soft/bowtie-1.3.1-linux-x86_64/bowtie-build
export KMC_TOOLS_BIN=/path/to/kmc_tools
```

## 输入文件

### 蛋白质序列和参考基因组

准备目标蛋白质 FASTA，以及每个材料各自的参考基因组 FASTA。例如：

```text
test/TaFT-B1_pep.fasta
/data/genome/CS.fa
/data/genome/JM47.fa
/data/genome/S4185.fa
```

蛋白质比对结果必须手动或通过脚本过滤：保留目标染色体（`TaFT-B1` 示例中为 `Chr7B`）上的命中，并排除同源染色体上的命中和明显错误的匹配。

### 区域文件

`extract_sequences_by_region.sh` 使用每行一个区域的格式，坐标为从 1 开始的闭区间：

```text
Chr7B:10013330-10014587
Chr7B:10025677-10026934
Chr7B:10038024-10039281
```

这些区域文件用于提取序列。Bowtie 筛选使用下文带 target_id、genome_id 和坐标系统的区域表；应列出每个参考中属于该目标的全部 copy。遗漏目标 copy 会使相关 k-mer 被判为非目标命中。

### 拷贝分组表

建议维护一个以制表符分隔的 `all_candidate_copy_list.tsv` 文件：

```text
copy_id sample group fasta
CS_Chr7B_copy01 CS 7B copy_fasta/CS_Chr7B_copy01.fasta
JM47_Chr7B_copy01 JM47 7B copy_fasta/JM47_Chr7B_copy01.fasta
CS_Chr7A_copy01 CS 7A copy_fasta/CS_Chr7A_copy01.fasta
```

第一行为表头。在以下示例中，`group` 也可以表示亚基因组、功能类别或单倍型分组。每个拷贝必须分别保存为独立的 FASTA 文件，且每个文件恰好包含一条非空序列记录；这些 FASTA 文件直接用于内存交集计算。

## 快速开始：鉴定候选拷贝

以下命令均在 `TKCN/` 目录中运行。首先，为每个材料构建 miniprot 索引并生成 GFF：

```bash
cd TKCN

bash src/run_miniport.sh \
  -q test/TaFT-B1_pep.fasta \
  -v CS \
  -g /data/genome/CS.fa \
  -o result/TaFT-B1_CS.gff \
  -t 16
```

设置 `GENOME_DIR=/data/genome` 后可以省略 `-g`。检查多个材料的 GFF 文件后，为每个材料准备候选区域文件，并移除同源染色体上不需要的命中。

## 提取拷贝序列

该脚本会同时生成一个包含所有指定区域的 FASTA，以及每个区域各自独立的 FASTA。`-p` 用于指定拆分序列的文件名前缀；其父目录必须已存在：

```bash
mkdir -p test/TaFT-B1_split_fasta result

bash src/extract_sequences_by_region.sh \
  -f /data/genome/CS.fa \
  -r test/TaFT-B1_CS_region.txt \
  -o result/TaFT-B1_CS_region.fa \
  -p test/TaFT-B1_split_fasta/CS
```

输出文件名示例如下：

```text
test/TaFT-B1_split_fasta/CS_Chr7B_10013330_10014587.fasta
```

将这些拆分 FASTA 文件的路径加入 `all_candidate_copy_list.tsv`。

## 构建目标特异性 panel

### 在内存中计算目标组共有 k-mer

将当前目标组的 copy FASTA 单独放入一个目录，然后运行 `shared_kmers_memory.py`。本示例中，分组表的 FASTA 路径相对于运行命令的 `TKCN/` 目录；用 `readlink -f` 将软链接源转换为绝对路径：

```bash
mkdir -p target_7B/input
while IFS=$'\t' read -r copy_id sample group fasta; do
  [[ -z "$copy_id" || "$copy_id" == "copy_id" || "$copy_id" == \#* ]] && continue
  [[ "$group" == "7B" ]] || continue
  ln -s "$(readlink -f "$fasta")" "target_7B/input/${copy_id}.fasta"
done < all_candidate_copy_list.tsv

python3 src/shared_kmers_memory.py \
  --input-dir target_7B/input \
  --pattern '*.fasta' \
  --target-id TaFTB1_7B \
  -k 31 \
  --out-dir target_7B/shared
```

`--pattern` 为不递归匹配的 glob，默认 `*.fasta`；每个匹配文件必须恰好包含一个 copy，长度至少为 `k`。每个目标使用独立、干净的输入目录，每次运行使用尚不存在的输出目录。多个目标分别计算，进入联合 Bowtie1 筛选前统一使用相同的 `k`。内存需求取决于单个 copy 的去重 k-mer 集合及其运行交集，无需构建全基因组 k-mer 数据库。

```text
target_7B/shared/TaFTB1_7B_shared.k31.tsv  无表头 kmer/count 两列，count 恒为 1
target_7B/shared/TaFTB1_7B_shared.k31.fa   共有 canonical k-mer FASTA
target_7B/shared/per_copy.tsv             每个 copy 的大小、跳过窗口数及运行交集大小
target_7B/shared/report.json              输入路径与哈希、k、copy 数、耗时和峰值内存
```

TSV 中的 `1` 仅表示存在，不是 copy 内出现次数，也不是重测序深度。若交集为空，仍会输出空的 TSV/FASTA，继续前应检查报告。该 TSV 可直接作为 Bowtie1 脚本的两列输入，无需 KMC 建库或 dump 转换。

### Bowtie1 批量筛选共有特异性 k-mer

准备以下三个**制表符分隔**的表（示例中的空格需替换成 Tab）。路径相对于各自清单文件所在目录。每个目标一行，`shared_tsv` 为 `shared_kmers_memory.py` 输出的两列 TSV；所有目标使用相同 k。脚本会从这些 TSV 自动生成 `all_targets.shared.fa`，无需手动拼接 FASTA。

`targets.tsv`:

```text
target_id shared_tsv
TaFTB1_7B target_7B/shared/TaFTB1_7B_shared.k31.tsv
```

`references.tsv`:

```text
genome_id fasta index_prefix
CS /data/genome/CS.fa /data/bowtie_index/CS
JM47 /data/genome/JM47.fa /data/bowtie_index/JM47
```

`target_regions.tsv`:

```text
target_id genome_id reference_name start end coordinate_system
TaFTB1_7B CS Chr7B 10013330 10014587 1-based-closed
TaFTB1_7B CS Chr7B 10025677 10026934 1-based-closed
TaFTB1_7B JM47 Chr7B 20000 21500 1-based-closed
```

区域坐标仅用于演示，请替换为各参考中全部目标 copy 的实际坐标。支持 `1-based-closed` 和 `0-based-half-open`，染色体名称必须与参考一致；必须提供每个目标×参考的区域。已有完整 Bowtie1 `.ebwt` 或 `.ebwtl` 索引会被复用；缺失时由 `bowtie-build` 自动建立。`index_prefix` 可留空（保留该列 Tab），此时在输出目录中建索引。

```bash
python3 src/bowtie_shared_specific.py \
  --targets targets.tsv \
  --references references.tsv \
  --regions target_regions.tsv \
  --out-dir bowtie_specific -t 16
```

输出 `bowtie_specific/<target>.shared-specific.tsv`（无表头 `kmer count`，保留输入计数）及 FASTA，可直接进入下述质量过滤。`classification.tsv` 包含每个参考的 `target_only`、`non_target_only`、`mixed`、`unmapped` 状态和保留标志；`<genome>.hits.tsv` 记录每条命中的 0-based half-open 区间和链方向。`summary.tsv` 汇总数量，`SUCCESS` 标志全部参考成功完成；失败目录不可用作最终结果。输出目录必须不存在；用新目录重跑。默认删除成功解析后的 SAM，可用 `--keep-sam` 保留。没有任何命中也会被剔除；多次目标内命中允许保留。旧的屏蔽/背景并集/差集脚本仅保留供旧流程使用，新流程不调用它们。

### Panel 质量过滤

将 Bowtie1 输出的 shared-specific TSV 过滤为适合 reads 计数的 TSV 和 FASTA panel：

```bash
python3 src/filter_target_specific_kmers.py \
  bowtie_specific/TaFTB1_7B.shared-specific.tsv \
  --out-prefix target_7B/TaFTB1_7B.panel \
  --anchor-fasta test/TaFT-B1.fasta \
  --min-spacing 10 \
  --id-prefix TaFTB1_7B
```

输出文件：

```text
TaFTB1_7B.panel.tsv         保留的 k-mer 及其质量指标
TaFTB1_7B.panel.fa          用于 reads 计数的 FASTA panel
TaFTB1_7B.panel.reject.tsv  被过滤的 k-mer 及过滤原因
```

默认情况下，脚本会根据 GC 含量、长同聚物、低 Shannon 熵和简单重复过滤 k-mer。可通过 `--min-gc`、`--max-gc`、`--max-homopolymer`、`--min-entropy` 和 `--max-repeat-period` 调整阈值。

## Reads 计数与拷贝数估算

构建 panel 后，使用 KMC 统计每个重测序样本中的 k-mer，并准备以下两类计数表：

```text
single_copy_gene_kmer_num/<marker>/*.realcount.tsv
multi_copy_gene_kmer_num/TaFT-B1_copy_num/*.realcount.tsv
```

每个计数表应为无表头的 KMC dump TSV，至少包含两列，第二列为数值型 `count`；文件名中的样本名称将成为输出中的 `sample`。下面是使用 KMC 统计 reads 的示例：

```bash
mkdir -p sample_kmc sample_count kmc_tmp

kmc -k31 -ci1 -cs255 -fm -t4 \
  target_7B/TaFTB1_7B.panel.fa \
  target_7B/TaFTB1_7B.panel.k31 kmc_tmp

printf '%s\n' sample_R1.fastq.gz sample_R2.fastq.gz > sample.reads.list

kmc -k31 -ci1 -cs1000000000 -fq -t16 \
  @sample.reads.list \
  sample_kmc/sample.k31 kmc_tmp

kmc_tools simple \
  sample_kmc/sample.k31 \
  target_7B/TaFTB1_7B.panel.k31 \
  intersect \
  sample_kmc/sample.target.k31 -ocleft -cs1000000000

kmc_tools transform sample_kmc/sample.target.k31 dump -s \
  sample_count/sample.TaFT-B1.realcount.tsv
```

此处 KMC 用于 copy FASTA 交集和 Bowtie1 筛选完成后的 reads 计数。必须将 reads 数据库放在左侧，并用 `-ocleft` 保留其计数；默认的最小值规则会将计数限制为 panel 中的出现次数，无需再重复 panel FASTA。所有步骤应使用相同的 canonical `k=31`，并设置足够的计数容量；先前过滤或计数饱和造成的损失无法恢复。dump 无表头；需要完整 panel 计数表时，在 reads 数据库保留 singleton（`-ci1`）的前提下，将未出现的 panel k-mer 补为零。

如果项目已有 reads 计数流程，可直接使用其生成的 `*.realcount.tsv` 文件，无需再次计数。

### 1. 估算单拷贝标记的 NB 深度

首先，为每个样本中的每个标记估算 NB/Poisson 模型参数：

```bash
cd TKCN

python3 src/estimate_depth_from_realcount_nb.py marker \
  ../single_copy_gene_kmer_num/*/*.realcount.tsv \
  > ../single_copy_gene_kmer_num/single_copy_marker_nb.tsv
```

然后使用多个标记，联合估算每个样本的单拷贝深度：

```bash
python3 src/estimate_depth_from_realcount_nb.py joint \
  --marker-stats ../single_copy_gene_kmer_num/single_copy_marker_nb.tsv \
  > ../single_copy_gene_kmer_num/single_copy_depth_nb_joint.tsv
```

`single_copy_depth_nb_joint.tsv` 中的 `nb_depth_estimate` 是后续拷贝数估算使用的一份拷贝参考深度。

### 2. 使用 NB 模型估算目标基因拷贝数

推荐使用 NB 模型，并通过似然值在相邻整数之间进行强制判定：

```bash
bash src/estimate_copy_from_realcount_depth.sh \
  --target ../multi_copy_gene_kmer_num/TaFT-B1_copy_num \
  --depth ../single_copy_gene_kmer_num/single_copy_depth_nb_joint.tsv \
  --target-method nb \
  --copy-call-method nb_likelihood \
  --out ../copy_num_estimate/TaFT-B1_copy_estimate_nb_nb_forced.tsv
```

核心计算为：

```text
copy_estimate = target_nb_mean_count / single_copy_depth
```

在输出结果中，`copy_call` 是推荐的整数拷贝数，`copy_estimate` 是连续估计值；`copy_call_loglik_delta` 越大，通常表示整数判定越稳定。当目标计数的方差接近均值时，脚本使用 `poisson_limit`，但仍可进行 NB 似然判定。

### 3. 可选的模糊判定

如果不希望对接近半整数的结果进行强制判定，可使用 `--copy-half-window`：

```bash
bash src/estimate_copy_from_realcount_depth.sh \
  --target ../multi_copy_gene_kmer_num/TaFT-B1_copy_num \
  --depth ../single_copy_gene_kmer_num/single_copy_depth_nb_joint.tsv \
  --target-method nb \
  --copy-half-window 0.1 \
  --out ../copy_num_estimate/TaFT-B1_copy_estimate_nb_nb_ambiguous.tsv
```

例如，`--copy-half-window 0.1` 会将接近 `2.40-2.60` 的值标记为 `ambiguous_half_integer`，并将 `copy_call` 设为 `NA`。主分析推荐使用强制判定结果；模糊判定结果可用于敏感性检查和低置信度注释。

NB 输入格式、输出列和模型细节以本目录中的 `src/estimate_depth_from_realcount_nb.py` 和 `src/estimate_copy_from_realcount_depth.sh --help` 为准；这两个脚本在本仓库中是自包含实现。

## 脚本索引

| 脚本 | 功能 |
| --- | --- |
| `src/run_miniport.sh` | 为一个材料构建 miniprot 索引并生成 GFF |
| `src/extract_sequences_by_region.sh` | 按区域提取合并 FASTA 和拆分后的拷贝 FASTA |
| `src/mask_fasta_regions.sh` | 旧流程：使用 `N` 屏蔽候选区域 |
| `src/batch_kmc_union.sh` | 旧流程：对输入目录中的 KMC 数据库求并集 |
| `src/batch_kmc_inter.sh` | 旧流程：对各 copy 的 KMC 数据库求交集 |
| `src/shared_kmers_memory.py` | 直接从独立 copy FASTA 在内存中求 canonical k-mer 交集 |
| `src/build_copy_background.sh` | 旧流程：构建非目标拷贝背景和完整背景 |
| `src/kmc_subtract_dump.sh` | 旧流程：计算 KMC 差集并导出 TSV |
| `src/kmc_inter_dump.sh` | 计算两个 KMC 数据库的交集并导出 TSV |
| `src/bowtie_shared_specific.py` | Bowtie1 全部精确命中、区域筛选及跨参考交集，构建 shared-specific panel |
| `src/filter_target_specific_kmers.py` | 过滤 k-mer 并输出 panel TSV/FASTA 文件 |
| `src/estimate_depth_from_realcount_nb.py` | 估算标记的 NB 参数和联合单拷贝深度 |
| `src/estimate_copy_from_realcount_depth.sh` | 根据目标计数和单拷贝深度估算整数及连续拷贝数 |

## 注意事项

- 一个 KMC 数据库由 `.kmc_pre` 和 `.kmc_suf` 两个文件组成；移动或复制时必须成对处理。
- 用于计算目标共有 k-mer 的输入目录只能包含当前目标组，每个 copy FASTA 一条序列记录。不要混入其他目标的 copy。
- 应使用独立材料或 reads 数据集验证 k-mer panel。低深度样本以及拷贝数估计值接近半整数的样本应标记为低置信度。
- 本目录中的脚本同时处理候选区域、k-mer panel 构建和 NB 拷贝数估算。生成的 `panels/` 数据保留在本地，并已从 Git 中排除。

## 测试数据

`test/TaFT-B1.fasta` 是 `filter_target_specific_kmers.py` 使用的锚定 FASTA，用于在 panel 输出中标注 k-mer 位置并执行 `--min-spacing` 稀疏化。它不是完整的参考基因组，不能作为 `extract_sequences_by_region.sh` 的 `-f` 输入。`test/` 中还包含蛋白质序列、区域列表以及样本/材料列表，可用于检查输入格式。运行流程之前，请确保已安装 `samtools`、`bowtie`、`bowtie-build`、`kmc`、`kmc_tools`、`miniprot`、GNU `awk` 和 Python 包 `scipy`，并根据本地环境调整命令路径。
