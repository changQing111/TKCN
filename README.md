# TKCN

TKCN (Targeted K-mer Copy Number Profiling) estimates copy numbers for duplicated genes or multi-copy gene families using k-mers shared by the target copies and specific to their genomic regions. Using wheat `TaFT-B1` as an example, the workflow is: candidate-copy identification → sequence extraction → in-memory shared k-mer intersection → Bowtie1 specificity filtering → panel QC → resequencing read counting → copy-number estimation.

## Method overview

For a target copy group `T`, define:

```text
K(copy) = set of canonical ACGT k-mers from one copy FASTA
target_shared = intersection of K(copy) across every copy in the target group
S(target,j) = shared k-mers with target hits and no off-target hits in reference j
target_specific = intersection of S(target,j) across all reference genomes
```

`shared_kmers_memory.py` enumerates each copy's canonical k-mers (the lexicographically smaller of a sequence and its reverse complement), removes duplicates within each copy, and intersects the sets in memory. It skips windows containing non-ACGT bases. No per-copy KMC databases are needed.

`bowtie_shared_specific.py` combines the shared k-mers from all targets into one FASTA with target identifiers, then maps it once to each complete, unmasked reference using Bowtie1 `-f -S -v 0 -a`. Both strands and all exact hits are considered. In each reference, a retained k-mer must have at least one target hit and no off-target hits; every full alignment interval must fit within one region belonging to that target. Only k-mers passing in every reference form the shared-specific panel. This replaces genome masking, background k-mer unions, and genome set subtraction; reusable Bowtie1 indexes provide the reference search step.

## Dependencies

The following software is required:

| Software | Purpose |
| --- | --- |
| `miniprot` | Search for candidate gene copies using protein sequences |
| `samtools` | Extract FASTA sequences by genomic region |
| `Bowtie1` (`bowtie`, `bowtie-build`) | Exact all-hit mapping and reference indexing |
| `KMC3` (`kmc`, `kmc_tools`) | Count resequencing reads and retrieve panel counts; not needed for copy FASTA intersections or Bowtie1 specificity filtering |
| `bash`, GNU `awk`, `python3`, `SciPy` | Run helper scripts, panel filtering, and NB models; shared-kmer and Bowtie1 scripts use the Python standard library only |

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

`extract_sequences_by_region.sh` uses one region per line in 1-based, closed-coordinate format:

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

The first line is a header. In the examples below, `group` can instead represent a subgenome, functional class, or haplotype group. Each copy must have its own FASTA file containing exactly one nonempty sequence record. These files are the direct inputs to the in-memory intersection.

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

## Build a target-specific panel

### Target shared k-mers in memory

Collect only the target group's copy FASTAs in one directory, then run `shared_kmers_memory.py`. In this example, manifest paths are relative to the `TKCN/` working directory; `readlink -f` makes symlink sources absolute:

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

`--pattern` is a nonrecursive glob (`*.fasta` by default); each matching file must contain one copy at least `k` bases long. Use a fresh input directory for each target and a new output directory for each run. Repeat this step separately for other targets, always using the same `k` before joint Bowtie1 filtering. Memory use depends on the unique k-mer sets of individual copies and their running intersection, rather than whole-genome k-mer databases.

```text
target_7B/shared/TaFTB1_7B_shared.k31.tsv  Headerless kmer/count table, count always 1
target_7B/shared/TaFTB1_7B_shared.k31.fa   Shared canonical k-mers in FASTA format
target_7B/shared/per_copy.tsv             Per-copy sizes, skipped windows, running intersection
target_7B/shared/report.json              Input paths/hashes, k, copy count, timing and peak RSS
```

The TSV count is a presence indicator, not within-copy multiplicity or resequencing depth. Empty intersections produce empty TSV/FASTA files; inspect the reports before proceeding. The TSV is directly compatible with the Bowtie1 script's two-column input format, so no KMC database or dump conversion is required.

### Batch shared-specific filtering with Bowtie1

Prepare three **tab-separated** tables (replace displayed spaces with tabs). Paths are relative to each manifest directory. Supply one row per target and the two-column TSV from `shared_kmers_memory.py`; all targets must use the same k. The script generates `all_targets.shared.fa` automatically from these TSVs; do not concatenate FASTAs manually.

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

Filter the Bowtie1 shared-specific TSV into a TSV and FASTA panel suitable for read counting:

```bash
python3 src/filter_target_specific_kmers.py \
  bowtie_specific/TaFTB1_7B.shared-specific.tsv \
  --out-prefix target_7B/TaFTB1_7B.panel \
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

Each count table must be a headerless KMC dump with at least two columns (`kmer`, numeric `count`). The sample name in each filename becomes the output `sample`. Example KMC read counting:

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

Here KMC is used for read counting, after the copy FASTA intersection and Bowtie1 filtering. Keep the reads database on the left and use `-ocleft` to retain its counts; the default minimum rule would cap counts at the panel multiplicity. No repeated panel FASTA is needed. Use the same canonical `k=31` convention throughout and sufficient counter capacity; earlier filtering or count saturation cannot be recovered. The dump is headerless; for a complete panel count table, add absent panel k-mers with count zero when the reads database preserves singleton counts (`-ci1`).

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

`nb_depth_estimate` in `single_copy_depth_nb_joint.tsv` is the one-copy reference used for subsequent copy-number estimation.

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

For NB input formats, output columns, and model details, consult `src/estimate_depth_from_realcount_nb.py` and `src/estimate_copy_from_realcount_depth.sh --help`; these are self-contained implementations in this repository.

## Script index

| Script | Function |
| --- | --- |
| `src/run_miniport.sh` | Build a miniprot index and generate a GFF for one accession |
| `src/extract_sequences_by_region.sh` | Extract a combined FASTA and split copy FASTA files by region |
| `src/mask_fasta_regions.sh` | Legacy: Mask candidate regions with `N` |
| `src/batch_kmc_union.sh` | Legacy: Take the union of KMC databases in an input directory |
| `src/batch_kmc_inter.sh` | Legacy: Take the intersection of per-copy KMC databases |
| `src/shared_kmers_memory.py` | Intersect canonical k-mer sets directly from individual copy FASTAs in memory |
| `src/build_copy_background.sh` | Legacy: Build the non-target-copy background and complete background |
| `src/kmc_subtract_dump.sh` | Legacy: Calculate a KMC set difference and export a TSV |
| `src/bowtie_shared_specific.py` | Bowtie1 shared-specific panel construction |
| `src/filter_target_specific_kmers.py` | Filter k-mers and write panel TSV/FASTA files |
| `src/estimate_depth_from_realcount_nb.py` | Estimate marker NB parameters and joint single-copy depth |
| `src/estimate_copy_from_realcount_depth.sh` | Estimate integer and continuous copy numbers from target counts and single-copy depth |

## Notes

- A KMC database consists of `.kmc_pre` and `.kmc_suf` files; always move or copy them together.
- The input for target shared k-mers must contain only the current target group, with one record per copy FASTA. Do not mix other targets into the same input directory.
- Validate k-mer panels using independent accessions or read datasets. Flag low-depth samples and copy estimates near half-integers as low confidence.
- The scripts in this directory cover candidate regions, k-mer-panel construction, and NB copy-number estimation. Generated `panels/` data are kept locally and excluded from Git.

## Test data

`test/TaFT-B1.fasta` is the anchor FASTA used by `filter_target_specific_kmers.py` to annotate k-mer positions in panel output and to apply `--min-spacing` thinning. It is not a complete reference genome and cannot be used as the `-f` input to `extract_sequences_by_region.sh`. `test/` also contains protein sequences, region lists, and sample/accession lists that can be used to check input formats. Before running the workflow, ensure that `samtools`, `bowtie`, `bowtie-build`, `kmc`, `kmc_tools`, and `miniprot` are installed, and adjust commands for local paths as needed.
