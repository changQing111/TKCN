# TKCN

TKCN (Targeted K-mer Copy Number Profiling) estimates copy numbers for duplicated genes or multi-copy gene families using k-mers that are shared by the target copies and specific relative to the background. Using wheat `TaFT-B1` as an example, this document describes the complete workflow: candidate-copy identification, sequence extraction, KMC set operations, and panel construction.

## Method overview

For a target copy group `T`, define:

```text
target_shared = intersection of k-mers across every copy in the target group
S(target,j) = shared k-mers with target hits and no off-target hits in reference j
target_specific = intersection of S(target,j) across all reference genomes
```

Map to complete, unmasked references using Bowtie1 `-f -S -v 0 -a`. Every full alignment interval must lie inside one region belonging to the query target. KMC intersections across target copies remain unchanged.

## Dependencies

The following software is required:

| Software | Purpose |
| --- | --- |
| `miniprot` | Search for candidate gene copies using protein sequences |
| `samtools` | Extract FASTA sequences by genomic region |
| `Bowtie1` (`bowtie`, `bowtie-build`) | Exact all-hit mapping and reference indexing |
| `KMC3` (`kmc`, `kmc_tools`) | Build k-mer databases and perform set operations |
| `bash`, `awk`, `python3` | Run helper scripts and filter panels |

All scripts support `-h/--help`. If a program is not on `PATH`, override its location with environment variables:

```bash
export MINIPROT_BIN=/path/to/miniprot
export SAMTOOLS_BIN=/path/to/samtools
export BOWTIE_BIN=/home/changq/my_data/bio_soft/bowtie-1.3.1-linux-x86_64/bowtie
export BOWTIE_BUILD_BIN=/home/changq/my_data/bio_soft/bowtie-1.3.1-linux-x86_64/bowtie-build
export KMC_TOOLS_BIN=/path/to/kmc_tools
```

## Input files

### Protein sequences and reference genomes

Prepare a target protein FASTA and one reference-genome FASTA for each accession. For example:

```text
test/TaFT-B1_pep.fasta
/data/genome/CS.fa
/data/genome/JM47.fa
/data/genome/S4185.fa
```

Protein-alignment results must be filtered manually or by script: retain hits on the target chromosome (`Chr7B` in the `TaFT-B1` example), and exclude homoeologous-chromosome hits and obvious incorrect matches.

### Region files

`extract_sequences_by_region.sh` and `mask_fasta_regions.sh` use one region per line in 1-based, closed-coordinate format:

```text
Chr7B:10013330-10014587
Chr7B:10025677-10026934
Chr7B:10038024-10039281
```

These region files are used for sequence extraction. Bowtie filtering uses the target/genome-specific region table below. Include all target copies in each reference; omitted target copies cause their hits to be classified as off-target.

### Copy-group table

We recommend maintaining a tab-separated `all_candidate_copy_list.tsv` file:

```text
copy_id sample group fasta
CS_Chr7B_copy01 CS 7B copy_fasta/CS_Chr7B_copy01.fasta
JM47_Chr7B_copy01 JM47 7B copy_fasta/JM47_Chr7B_copy01.fasta
CS_Chr7A_copy01 CS 7A copy_fasta/CS_Chr7A_copy01.fasta
```

The first line is a header. In the examples below, `group` can instead represent a subgenome, functional class, or haplotype group. Each copy must have its own FASTA file and KMC database; do not build one database by combining different copies.

## Quick start: identify candidate copies

Run the following commands from `TKCN/`. First, build a miniprot index and generate a GFF for each accession:

```bash
cd TKCN

bash src/run_miniport.sh \
  -q test/TaFT-B1_pep.fasta \
  -v CS \
  -g /data/genome/CS.fa \
  -o result/TaFT-B1_CS.gff \
  -t 16
```

You may omit `-g` by setting `GENOME_DIR=/data/genome`. After inspecting GFF files for multiple accessions, prepare a candidate-region file for each accession and remove unwanted hits on homoeologous chromosomes.

## Extract copy sequences

This script generates both one FASTA containing all requested regions and one separate FASTA for each region. `-p` specifies the prefix for split sequences; its parent directory must already exist:

```bash
mkdir -p test/TaFT-B1_split_fasta result

bash src/extract_sequences_by_region.sh \
  -f /data/genome/CS.fa \
  -r test/TaFT-B1_CS_region.txt \
  -o result/TaFT-B1_CS_region.fa \
  -p test/TaFT-B1_split_fasta/CS
```

Output filenames resemble:

```text
test/TaFT-B1_split_fasta/CS_Chr7B_10013330_10014587.fasta
```

Add the paths of these split FASTA files to `all_candidate_copy_list.tsv`.

## Build k-mer databases

### Build one database for each candidate copy

```bash
mkdir -p kmc_batch_db kmc_batch_tmp
```

```bash
while IFS=$'\t' read -r copy_id sample group fasta; do
  [[ -z "$copy_id" || "$copy_id" == "copy_id" || "$copy_id" == \#* ]] && continue
  kmc -k31 -ci1 -fm -t8 "$fasta" \
    "kmc_batch_db/${copy_id}.k31" kmc_batch_tmp
done < all_candidate_copy_list.tsv
```

## Build a target-specific panel

### Target shared k-mers

Select the copy databases for the target group and calculate their intersection:

```bash
mkdir -p target_7B_kmc
awk -v group=7B 'NR > 1 && $3 == group {print $1}' \
  all_candidate_copy_list.tsv > target_7B_kmc/copies.list

mkdir -p target_7B_kmc/input
while read -r copy_id; do
  ln -sf "../../kmc_batch_db/${copy_id}.k31.kmc_pre" \
    "target_7B_kmc/input/${copy_id}.k31.kmc_pre"
  ln -sf "../../kmc_batch_db/${copy_id}.k31.kmc_suf" \
    "target_7B_kmc/input/${copy_id}.k31.kmc_suf"
done < target_7B_kmc/copies.list

bash src/batch_kmc_inter.sh \
  -i target_7B_kmc/input \
  -o target_7B_kmc/TaFTB1_7B.target_shared.k31
```

`batch_kmc_inter.sh` intersects all KMC databases in the input directory. Therefore, first create an input directory containing only the current target group. Alternatively, calculate the intersection explicitly in steps with `kmc_tools`:

```bash
kmc_tools simple \
  kmc_batch_db/CS_Chr7B_copy01.k31 \
  kmc_batch_db/JM47_Chr7B_copy01.k31 \
  intersect target_7B_kmc/step1.k31

kmc_tools simple \
  target_7B_kmc/step1.k31 \
  kmc_batch_db/S4185_Chr7B_copy01.k31 \
  intersect target_7B_kmc/TaFTB1_7B.target_shared.k31
```

### Batch shared-specific filtering with Bowtie1

Prepare three **tab-separated** tables (replace displayed spaces with tabs). Paths are relative to each manifest directory. Supply one row per target and the two-column KMC intersection dump; all targets must use the same k.

`targets.tsv`:

```text
target_id shared_tsv
TaFTB1_7B target_7B_kmc/TaFTB1_7B.target_shared.k31.tsv
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

Replace example coordinates with all actual target-copy intervals in each reference. Supported coordinate systems: `1-based-closed` and `0-based-half-open`. Chromosome names must match the reference; every target/reference pair requires regions. Complete Bowtie1 `.ebwt` or `.ebwtl` indexes are reused; missing indexes are built. A blank `index_prefix` field builds under the output directory.

```bash
python3 src/bowtie_shared_specific.py \
  --targets targets.tsv \
  --references references.tsv \
  --regions target_regions.tsv \
  --out-dir bowtie_specific -t 16
```

Outputs are `bowtie_specific/<target>.shared-specific.tsv` (headerless `kmer count`, preserving input counts) and FASTA, ready for QC below. `classification.tsv` reports `target_only`, `non_target_only`, `mixed`, and `unmapped` status plus retention per reference. `<genome>.hits.tsv` records every alignment with 0-based half-open coordinates and strand. `summary.tsv` reports counts; `SUCCESS` marks completion across all references. Failed output directories are not final results. Use a new output directory for every run. SAM files are removed after successful parsing unless `--keep-sam` is set. Unmapped queries are excluded; multiple target-only hits are allowed. Masking/background/subtraction scripts remain for legacy use only.

### Panel quality filtering

Filter the KMC dump table into a TSV and FASTA panel suitable for read counting:

```bash
python3 src/filter_target_specific_kmers.py \
  bowtie_specific/TaFTB1_7B.shared-specific.tsv \
  --out-prefix target_7B_kmc/TaFTB1_7B.panel \
  --anchor-fasta test/TaFT-B1.fasta \
  --min-spacing 10 \
  --id-prefix TaFTB1_7B
```

Outputs:

```text
TaFTB1_7B.panel.tsv         Retained k-mers and quality metrics
TaFTB1_7B.panel.fa          FASTA panel for read counting
TaFTB1_7B.panel.reject.tsv  Rejected k-mers and their reasons
```

By default, the script filters k-mers by GC content, long homopolymers, low Shannon entropy, and simple repeats. Adjust thresholds with `--min-gc`, `--max-gc`, `--max-homopolymer`, `--min-entropy`, and `--max-repeat-period`.

## Read counting and copy-number estimation

After building the panel, count k-mers in each resequencing sample with KMC and prepare two kinds of count tables:

```text
single_copy_gene_kmer_num/<marker>/*.realcount.tsv
multi_copy_gene_kmer_num/TaFT-B1_copy_num/*.realcount.tsv
```

Each count table must contain at least `kmer` and `count` columns. The sample name in each filename becomes the output `sample`. Example KMC read counting:

```bash
mkdir -p sample_kmc sample_count

printf '%s\n' sample_R1.fastq.gz sample_R2.fastq.gz > sample.reads.list

kmc -k31 -ci1 -fq -t16 \
  @sample.reads.list \
  sample_kmc/sample.k31 kmc_tmp

kmc_tools simple \
  sample_kmc/sample.k31 \
  target_7B_kmc/TaFTB1_7B.panel.repeat10000.k31 \
  intersect \
  sample_kmc/sample.target.k31

kmc_tools transform sample_kmc/sample.target.k31 dump -s \
  sample_count/sample.TaFT-B1.realcount.tsv
```

`TaFTB1_7B.panel.repeat10000.k31` must be built from a panel FASTA repeated enough times that the repeat count exceeds the expected read k-mer depth; otherwise, genuine counts will not be retained in the KMC intersection.

If the project already has a read-counting workflow, you can use its generated `*.realcount.tsv` files directly instead of counting again.

### 1. Estimate NB depth for single-copy markers

First estimate NB/Poisson model parameters for every marker in every sample:

```bash
cd TKCN

python3 src/estimate_depth_from_realcount_nb.py marker \
  ../single_copy_gene_kmer_num/*/*.realcount.tsv \
  > ../single_copy_gene_kmer_num/single_copy_marker_nb.tsv
```

Then jointly estimate the single-copy depth for each sample using multiple markers:

```bash
python3 src/estimate_depth_from_realcount_nb.py joint \
  --marker-stats ../single_copy_gene_kmer_num/single_copy_marker_nb.tsv \
  > ../single_copy_gene_kmer_num/single_copy_depth_nb_joint.tsv
```

`nb_depth_estimate` in `single_copy_depth_nb_joint.tsv` is the one-copy reference used for subsequent copy-number estimation. For marker input-directory layout and optional parameters, see `../docs/nb-copy-number-estimation-method.md`.

### 2. Estimate target-gene copy number with NB models

We recommend the NB model with likelihood-based forced calls between adjacent integers:

```bash
bash src/estimate_copy_from_realcount_depth.sh \
  --target ../multi_copy_gene_kmer_num/TaFT-B1_copy_num \
  --depth ../single_copy_gene_kmer_num/single_copy_depth_nb_joint.tsv \
  --target-method nb \
  --copy-call-method nb_likelihood \
  --out ../copy_num_estimate/TaFT-B1_copy_estimate_nb_nb_forced.tsv
```

The core calculation is:

```text
copy_estimate = target_nb_mean_count / single_copy_depth
```

In the output, `copy_call` is the recommended integer copy number, `copy_estimate` is the continuous estimate, and a larger `copy_call_loglik_delta` generally indicates a more stable integer call. When the variance of target counts is close to the mean, the script uses `poisson_limit`, while still allowing NB likelihood calling.

### 3. Optional ambiguous calls

If you do not want to force calls near half-integers, use `--copy-half-window`:

```bash
bash src/estimate_copy_from_realcount_depth.sh \
  --target ../multi_copy_gene_kmer_num/TaFT-B1_copy_num \
  --depth ../single_copy_gene_kmer_num/single_copy_depth_nb_joint.tsv \
  --target-method nb \
  --copy-half-window 0.1 \
  --out ../copy_num_estimate/TaFT-B1_copy_estimate_nb_nb_ambiguous.tsv
```

For example, `--copy-half-window 0.1` labels values near `2.40-2.60` as `ambiguous_half_integer`, with `copy_call` set to `NA`. The forced-call output is recommended for the primary analysis; use the ambiguous-call output for sensitivity checks and low-confidence annotations.

For NB input formats, output columns, and model details, see `../docs/nb-copy-number-estimation-method.md`. The root-level scripts documented there have identically named entry points in `TKCN/src`; these call the project's shared implementation to prevent divergence between two algorithms.

## Script index

| Script | Function |
| --- | --- |
| `src/run_miniport.sh` | Build a miniprot index and generate a GFF for one accession |
| `src/extract_sequences_by_region.sh` | Extract a combined FASTA and split copy FASTA files by region |
| `src/mask_fasta_regions.sh` | Legacy: Mask candidate regions with `N` |
| `src/batch_kmc_union.sh` | Take the union of KMC databases in an input directory |
| `src/batch_kmc_inter.sh` | Take the intersection of KMC databases in an input directory |
| `src/build_copy_background.sh` | Legacy: Build the non-target-copy background and complete background |
| `src/kmc_subtract_dump.sh` | Legacy: Calculate a KMC set difference and export a TSV |
| `src/bowtie_shared_specific.py` | Bowtie1 shared-specific panel construction |
| `src/filter_target_specific_kmers.py` | Filter k-mers and write panel TSV/FASTA files |
| `src/estimate_depth_from_realcount_nb.py` | Estimate marker NB parameters and joint single-copy depth |
| `src/estimate_copy_from_realcount_depth.sh` | Estimate integer and continuous copy numbers from target counts and single-copy depth |

## Notes

- A KMC database consists of `.kmc_pre` and `.kmc_suf` files; always move or copy them together.
- The input for target shared k-mers must contain only the current target group. Do not put all copy databases directly in one intersection directory.
- Validate k-mer panels using independent accessions or read datasets. Flag low-depth samples and copy estimates near half-integers as low confidence.
- The scripts in this directory handle candidate regions and k-mer-panel construction only. NB copy-number statistics scripts and tests are located in the repository root.

## Test data

`test/TaFT-B1.fasta` is the anchor FASTA used by `filter_target_specific_kmers.py` to annotate k-mer positions in panel output and to apply `--min-spacing` thinning. It is not a complete reference genome and cannot be used as the `-f` input to `extract_sequences_by_region.sh`. `test/` also contains protein sequences, region lists, and sample/accession lists that can be used to check input formats. Before running the workflow, ensure that `samtools`, `bowtie`, `bowtie-build`, `kmc`, `kmc_tools`, and `miniprot` are installed, and adjust commands for local paths as needed.
