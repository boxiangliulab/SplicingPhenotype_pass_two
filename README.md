# bulk sQTL workflow for eQTLGen: 0. splicing phenotype

We implement the sQTL workflow as a set of sequential per-cohort stages. This is the cookbook for **Stage 0: splicing phenotype preparation** (`0_SplicingPhenotype`).

Starting from RNA-seq **FASTQ**, this stage STAR-aligns the reads **against a shared splice-junction database** (see below), quantifies both junction (intron excision) and intron-retention (IR) splice sites using [ISSAC](https://github.com/boxiangliulab/ISSAC/tree/ISSAC_analysis), then QQ-normalises both modalities and **combines them into one raw combined splice-site matrix**. Junction and IR splice sites are processed with their own filters but concatenated into a single splice-site matrix (features carry a `junction_` / `IR_` ID prefix) so that downstream stages QC, compute splicing PCs, and HASE-encode the genotype exactly once. **No sample QC is done here** — that moves into Stage 1 (`1_DataQC`) to run alongside genotype QC, mirroring eQTLGen where the phenotype QC lives in `1_DataQC`.

> **Shared splice-junction database.** The STAR alignment (Step 1) is run against a precomputed junction table, `data/GTExIntervalBios.junc.tab` (STAR `--sjdbFileChrStartEnd`), built by **pooling the splice-junction files across INTERVAL, BIOS, and GTEx**. Every cohort therefore aligns to the **same harmonised set of known junctions**, so the junctions quantified here are directly comparable across cohorts in the downstream meta-analysis — rather than each cohort discovering its own idiosyncratic junction set. This is also why all cohorts **start from FASTQ** (so the shared database is actually used) rather than from pre-aligned BAMs.

The deliverables of this stage — `{cohort_name}_combined_splice_site.txt.gz` (splicing) and `{cohort_name}_gene_counts.txt.gz` (expression) — are consumed directly by the next stage (`1_DataQC`), not uploaded to central analysis.

For this cookbook, we follow the style of [eQTLGen Phase II cookbook](https://eqtlgen.github.io/eqtlgen-web-site/eQTLGen-p2-cookbook.html) and aim to require minimal additional user actions.



## Eligibility

If you have RNA-seq data on your cohort(s), you will be able to run this analysis. **All cohorts start from RNA-seq FASTQ** — the STAR alignment (Step 1) is mandatory, so that every cohort aligns against the shared splice-junction database (INTERVAL + BIOS + GTEx) and produces the per-gene counts reused in Step 3a.



## Prerequisites

If you have run the workflow in eQTLGen Phase II cookbook, you should have all the following prerequisites. However, the current pipelines and configurations are tailored for use with the **PBS scheduler** by default. If your HPC uses a different scheduler, modify the job scripts to request the appropriate resources accordingly.

- High-Performance Computing (HPC) environment

- Singularity / Apptainer

- Nextflow



## Getting Started

### Setup

Make an analysis root folder for this project, e.g. `eQTLGen_splicing_phase1`. Please ensure the folder is on a data storage disk (not your home directory) to avoid issues.

```
mkdir eQTLGen_splicing_phase1
```

Inside `eQTLGen_splicing_phase1`, download the wrapper for this stage. Extract all stages into this one common parent folder so they end up as siblings — later stages read earlier outputs via `../../` relative paths, and the shared `tools/` launcher plus `singularitycache/` / `nextflowcache/` are shared:

```
cd eQTLGen_splicing_phase1
wget -O eQTLGen_splicing_0_SplicingPhenotype_ISSAC_bulk_wrapper.tar.gz     "https://www.dropbox.com/scl/fi/1a1p1irvtim6rubyvw5f3/eQTLGen_splicing_0_SplicingPhenotype_ISSAC_bulk_wrapper.tar.gz?rlkey=a9h36dtp61h3n3kt7uks7rc3j&dl=1"
wget -O eQTLGen_splicing_0_SplicingPhenotype_ISSAC_bulk_wrapper.tar.gz.md5 "https://www.dropbox.com/scl/fi/yhul9kzl1n6vsuaa1d0ck/eQTLGen_splicing_0_SplicingPhenotype_ISSAC_bulk_wrapper.tar.gz.md5?rlkey=6a953w50z96tnhbzc9mvim0gm&dl=1"
```

Please check the md5sum of the file and extract the content.

```
# This should print eQTLGen_splicing_0_SplicingPhenotype_ISSAC_bulk_wrapper.tar.gz: OK
md5sum -c eQTLGen_splicing_0_SplicingPhenotype_ISSAC_bulk_wrapper.tar.gz.md5
tar -xzvf eQTLGen_splicing_0_SplicingPhenotype_ISSAC_bulk_wrapper.tar.gz
```

You should be able to see the stage folder `0_SplicingPhenotype` (with `data/` and `singularity_img/`) plus the shared `tools/`. Below is the full content (displayed by command `tree ./`):

```
.
|-- 0_SplicingPhenotype
|   |-- data
|   |   |-- GRCh38.primary_assembly.genome.fa.gz
|   |   |-- gencode.v48.annotation.bed12
|   |   `-- gencode.v48.annotation.gtf
|   |-- singularity_img
|   |   |-- eqtlgen_splicing_star_sort.sif
|   |   |-- eqtlgen_splicing_junction_strandness.sif
|   |   `-- eqtlgen_splicing_phenotype_prep.sif
|   `-- splicing_phenotype                        <-- main work folder
|       |-- bin
|       |-- conf
|       |-- nextflow.config
|       |-- star_alignment_with_qc_pass_two.nf
|       |-- phenotype_prep.nf
|       |-- expression_quant.nf
|       |-- submit_1_star_alignment_template_pass_two.sh
|       |-- submit_2_strandness.sh
|       |-- submit_3_phenotype_prep_template.sh
|       `-- submit_3a_expression_quant_template.sh
`-- tools
    `-- nextflow
```

Test the nextflow executable. You need to load Java >=11 before testing (e.g. `module load [Java >=11 module name in your HPC]`).
We encountered some issues when using Java 21, with various versions of nextflow, it should work, but if it doesn't please try and older version of Java (Java 17 worked perfectly)

```
cd tools
./nextflow run info
# If this binary works on your machine, you will see:
# N E X T F L O W  ~  version 21.10.6
```

(Skip this step if nextflow works) If the nextflow does not work, please download and self-install the Nextflow executable and test. You might need to load Java >= 11 before running self-install.

```
wget -qO- https://get.nextflow.io | bash
./nextflow run info
```

Enter the main work directory, and prepare to run the pipelines.

```
cd ../0_SplicingPhenotype/splicing_phenotype
```

Inside the pipeline folder (`splicing_phenotype`) are the template scripts using the name format `submit_*_template.sh`. You should adjust these according to your data and environment (e.g. specify the path to your input folder, add required HPC modules), and save them. We specifically marked the sections to edit with the tag: "USER CONFIGURATION REQUIRED".

To ease the process, we have pre-filled some paths in the templates, assuming that you use the recommended folder structure. **You only ever edit the lines inside the `USER CONFIGURATION REQUIRED` block, and within it only the values shown as `[ ... ]` placeholders.** Everything else (container paths, reference paths, resources, `-profile`) is pre-filled to the bundle's relative layout.

You should submit each pipeline from inside the pipeline folder.



## Data

The starting point is always **RNA-seq raw FASTQ files**. Run the four steps in order:

1. `submit_1_star_alignment_template_pass_two.sh` — STAR alignment against the shared junction database (+ per-gene counts),
2. `submit_2_strandness.sh` — infer library strand-specificity (see [Strand-specificity](#strand-specificity)),
3. `submit_3_phenotype_prep_template.sh` — junction + IR → combined splice-site matrix,
4. `submit_3a_expression_quant_template.sh` — merge STAR gene counts → gene-count matrix (can run in parallel with Step 3).

The genome assembly used for alignment is **hg38** (bundled `data/` FASTA + `gencode.v48` GTF). STAR also writes the XS tag on spliced reads (needed for junction strand assignment when the library is unstranded) — handled automatically because alignment is always done here, not supplied as pre-made BAMs.



## Configure and run

### Step 1: STAR alignment (from FASTQ)

Modify the job file `submit_1_star_alignment_template_pass_two.sh` in the main work directory. Search for the section marked "USER CONFIGURATION REQUIRED", and at least provide the following:

1. `FASTQ_DIR` — the absolute path of your input FASTQ directory.

2. `READ_LENGTH` — read length / maximum read length of the sequencing read mates (e.g. `150` for paired-end 150). `sjdbOverhang` is set to `READ_LENGTH - 1`.

3. `SINGLE_END` / `FASTQ_PATTERN` — set `SINGLE_END=true` and `FASTQ_PATTERN='*.fastq.gz'` for single-end data (defaults are paired-end).

4. `NEXTFLOW_PROFILES` — always keep `singularity`; change `pbs` to your scheduler.

`SJDB_JUNCTIONS` is **pre-filled** to the shared junction database `../data/GTExIntervalBios.junc.tab` (STAR `--sjdbFileChrStartEnd`, built across INTERVAL + BIOS + GTEx) — do not change it unless you have a specific reason; keeping it is what makes junctions comparable across cohorts.

```
qsub submit_1_star_alignment_template_pass_two.sh   # or equivalent command for your scheduler
```

This runs FASTQ QC & trimming (fastp) → STAR alignment (against the shared junction database, with `--quantMode GeneCounts` so per-sample `ReadsPerGene.out.tab` is produced for Step 3a) → QC reports, writing BAMs to `../output/bam_alignment/aligned_bam/` and a QC summary to `../output/bam_alignment/star_alignment_qc_summary.tsv`. The script prints an `awk` command to review PASS / WARNING / FAIL samples against uniquely-mapped and total-mapped thresholds.

### Step 2: (Optional) Infer strandness from BAM files

If you are not sure what strand parameter to use in the phenotype pipeline, run the helper script described in [Strand-specificity](#strand-specificity).

### Step 3: Extract phenotypes from BAM files (junction + IR → combined splice-site matrix)

Modify the job file `submit_3_phenotype_prep_template.sh` in the main work directory. Search for the section marked "USER CONFIGURATION REQUIRED", and at least provide the following:

1. `COHORT_NAME` — your cohort name.

2. `STRAND` — strand parameter: `XS` (unstranded), `RF` (reverse stranded), or `FR` (forward stranded).

3. `BAM_DIR` — pre-filled to the Step 1 STAR output (`../output/bam_alignment/aligned_bam`); no change needed.

4. `-profile` in the Nextflow command (at the bottom of the script) and the corresponding config file.

```
qsub submit_3_phenotype_prep_template.sh   # or equivalent command for your scheduler
```

This one pipeline (`phenotype_prep.nf`, one container with ISSAC-bulk compiled in) extracts and clusters junctions, computes IR, QQ-normalises both modalities, and **combines** them into the raw combined splice-site matrix consumed by Stage 1. Optional filter parameters (`--read_count_threshold`, `--min_intron_len`, `--sd_threshold`, IR toggles, etc.) are listed at the bottom of the template and only need changing in special cases.

### Step 3a: Expression quantification (gene counts → expression matrix)

Stage 1 now also needs a **raw gene-count expression matrix** (expression is used as the LD / permuted-data phenotype downstream). STAR (in `star_alignment_with_qc.nf`) already emits per-sample gene counts via `--quantMode GeneCounts` (`<sample>.ReadsPerGene.out.tab`), so this step just **merges** them — no re-quantification. Run `submit_3a_expression_quant_template.sh`; under "USER CONFIGURATION REQUIRED" provide `COHORT_NAME`, `COUNTS_DIR` (pre-filled to the STAR `aligned_bam` dir), `STRAND` (must match your library — selects the right `ReadsPerGene` column), and `-profile`.

```
qsub submit_3a_expression_quant_template.sh   # or equivalent for your scheduler
```

This runs `expression_quant.nf` (thin DSL2 companion to `phenotype_prep.nf`, same container) → `output/expression_matrix/{cohort_name}_gene_counts.txt.gz`. It can run in parallel with Step 3 (both read the same aligned BAMs). STAR GeneCounts were validated against featureCounts (concordant per gene), so reuse is safe.



## Strand-specificity

If you are not sure whether your RNA-seq library is stranded or not, please run the short diagnostic script `submit_2_strandness.sh` after Step 1 (it runs on the STAR-aligned BAMs). We use `infer_experiment.py` from [RSeQC](https://rseqc.sourceforge.net/). It samples several representative BAMs across the cohort, averages the results, checks consistency, and prints a recommended `--strand` value.

```
# BAM_DIR is pre-filled to the Step 1 STAR output; adjust headers/profile for your scheduler
qsub submit_2_strandness.sh   # or equivalent command for your scheduler
```

You can check the results in `prepare/strandness.txt` after the job completes. This usually takes less than 1 minute. For paired-end RNA-seq:

       Usually, strand-specific RNA-seq libraries are prepared using a dUTP-based protocol, resulting in reverse-stranded reads (--strand RF). In some cases (e.g. older datasets or non–strand-specific protocols), RNA-seq libraries may be unstranded (--strand XS).
       If the fraction of reads explained by 1+-,1-+,2++,2-- is much larger than that explained by 1++,1--,2+-,2-+ and is close to 1             => RF (reverse-stranded, dUTP-based)
       If 1+-,1-+,2++,2-- is close to 1++,1--,2+-,2-+ and both are close to 0.5             => XS (unstranded)

       Below are two example outputs:

       Example 1: Paired-end, non–strand-specific:
       This is PairEnd Data
       Fraction of reads failed to determine: 0.0172
       Fraction of reads explained by "1++,1--,2+-,2-+": 0.4903
       Fraction of reads explained by "1+-,1-+,2++,2--": 0.4925

       Example 2: Paired-end, strand-specific (reverse-stranded):
       This is PairEnd Data
       Fraction of reads failed to determine: 0.0072
       Fraction of reads explained by "1++,1--,2+-,2-+": 0.0487
       Fraction of reads explained by "1+-,1-+,2++,2--": 0.9441

If you encounter other scenarios or ambiguous results, please let us know. You can also refer to this [table](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/) for more details.



## Output

The pipeline generates an `output/` folder with the raw splice phenotypes and the combined splice-site matrix:

```
output/
├── bam_alignment/                                (from Step 1, if you start from FASTQ)
│   ├── 00_fastq_qc/
│   ├── 00_star_genome_index/
│   ├── aligned_bam/                              {sample_name}.bam + STAR logs
│   └── star_alignment_qc_summary.tsv
├── 01_junc/           {sample_name}.junc
├── 02_stat/           {sample_name}.stat
├── 03_phenotype/      {cohort_name}.filtered, {cohort_name}.prop         (junction phenotypes)
├── 04_IR/             {cohort_name}_IR.filtered, {cohort_name}_IR.prop   (intron retention)
├── phenotype_prep/    *.qqnorm_chr*.bed.gz                               (QQ-normed, gene-annotated BEDs)
├── combined_matrix/
│   └── {cohort_name}_combined_splice_site.txt.gz   <-- DELIVERABLE 1 (combined junction + IR splice-site matrix)
└── expression_matrix/
    └── {cohort_name}_gene_counts.txt.gz            <-- DELIVERABLE 2 (raw gene counts, from Step 3a)
```

The deliverables of this stage are:

- `combined_matrix/{cohort_name}_combined_splice_site.txt.gz` — the combined (junction + IR) splice-site matrix.
- `expression_matrix/{cohort_name}_gene_counts.txt.gz` — the raw gene-count matrix (Step 3a).

Both are consumed directly by the **next stage** (`1_DataQC`): splicing QC (which defines the analysed sample set) + a diagnostic expression QC + covariate preparation, alongside genotype QC. Neither is uploaded to central analysis at this stage — keep them in place (together with `03_phenotype/`, `04_IR/`) so `1_DataQC`'s pre-filled single `--splice_pheno_dir ../../0_SplicingPhenotype/output` can derive them all by `cohort_name`.



## Troubleshooting

You can check the scheduler-generated stdout/stderr log files (e.g. `3_phenotype_prep.o*` in PBS) or inspect `.nextflow.log`.

Please also refer to Section "Running, monitoring and debugging Nextflow pipelines" in [eQTLGen Phase II cookbook](https://eqtlgen.github.io/eqtlgen-web-site/eQTLGen-p2-cookbook.html).



## Citation

Nextflow:
[Di Tommaso, P., Chatzou, M., Floden, E. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316–319 (2017). https://doi.org/10.1038/nbt.3820](https://www.nature.com/articles/nbt.3820)

RSeQC:
[Wang L, Wang S, Li W. RSeQC: quality control of RNA-seq experiments. Bioinformatics. 2012 Aug 15;28(16):2184-5. doi: 10.1093/bioinformatics/bts356. Epub 2012 Jun 27. PMID: 22743226.](https://academic.oup.com/bioinformatics/article/28/16/2184/325191?login=false)

ISSAC:
[https://www.medrxiv.org/content/10.64898/2026.05.06.26352548v2.full.pdf](https://www.medrxiv.org/content/10.64898/2026.05.06.26352548v2.full.pdf)



## Contact

Please contact e1101920@u.nus.edu if you have any other questions.
