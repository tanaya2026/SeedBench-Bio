# SeedBench-Bio Model Response

## Metadata

- **Prompt variant:** baseline
- **Prompt file:** prompts/code/baseline.md
- **Model:**
provider: Anthropic
name: Claude Sonnet
thinking: medium
- **Task file:** benchmark_tasks/code_review/bulk_RNAseq_matrix.R
- **Ground truth file:** ground_truth/code_review/bulk_RNAseq_matrix.json
- **Timestamp:** 2026-08-05 17:10:00
- **Temperature:** 0

---

## Prompt sent to model

```markdown
Review the following bioinformatics code.

Identify any objective scientific, statistical, or computational errors.

For each error, explain why it is incorrect and suggest a correction.

Ignore coding style, formatting, and personal preferences.
```


```r
#!/usr/bin/env Rscript
# ==============================================================================
# bulk_RNAseq_matrix.R
#
# Pipeline: build a gene x sample count matrix from PRJEB22205
# (Drosophila melanogaster, paired-end bulk RNA-seq, ENA project PRJEB22205)
#
# NOTE FOR BENCHMARK USERS: this script contains INTENTIONALLY SEEDED ERRORS.
# Do not "fix" this file directly — it is a fixed stimulus for SeedBench-Bio.
# See ground_truth\bulk_RNAseq_matrix.json for the answer key.
# ==============================================================================

library(Rsubread)
library(data.table)

# ------------------------------------------------------------------------
# 1. Sample sheet
#    columns: run_accession, tissue, r1, r2, library_strandedness
# ------------------------------------------------------------------------
samples <- fread("sample_sheet.csv")

# ------------------------------------------------------------------------
# 2. Alignment
#    Reads are aligned directly from the raw fastq files as downloaded
#    from ENA.
# ------------------------------------------------------------------------
genome_index <- "/refs/hg38_star_index"          # STAR index directory
gtf_file     <- "/refs/Homo_sapiens.GRCh38.gtf"   # annotation

dir.create("star_out", showWarnings = FALSE)

for (i in seq_len(nrow(samples))) {
  r1 <- samples$r1[i]
  r2 <- samples$r2[i]
  out_prefix <- file.path("star_out", paste0(samples$run_accession[i], "_"))

  cmd <- sprintf(
    "STAR --genomeDir %s --readFilesIn %s %s --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix %s",
    genome_index, r1, r2, out_prefix
  )
  system(cmd)
}

# ------------------------------------------------------------------------
# 3. Deduplication
#    Remove PCR duplicates from each BAM before quantification so that
#    expression estimates aren't inflated by amplification artifacts.
# ------------------------------------------------------------------------
dir.create("dedup", showWarnings = FALSE)

for (i in seq_len(nrow(samples))) {
  bam_in  <- file.path("star_out", paste0(samples$run_accession[i], "_Aligned.sortedByCoord.out.bam"))
  bam_out <- file.path("dedup", paste0(samples$run_accession[i], "_dedup.bam"))
  cmd <- sprintf("picard MarkDuplicates I=%s O=%s REMOVE_DUPLICATES=true M=metrics.txt",
                 bam_in, bam_out)
  system(cmd)
}

# ------------------------------------------------------------------------
# 4. Quantification with featureCounts
# ------------------------------------------------------------------------
bam_files <- list.files("dedup", pattern = "_dedup.bam$", full.names = TRUE)

fc <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,
  nthreads = 8
)

# ------------------------------------------------------------------------
# 5. Assemble the gene x sample count matrix
# ------------------------------------------------------------------------
count_matrix <- fc$counts
colnames(count_matrix) <- samples$tissue   # label columns by tissue

cat("Built gene x sample matrix with", nrow(count_matrix), "genes and",
    ncol(count_matrix), "samples\n")

# ------------------------------------------------------------------------
# 6. Compare expression across tissues
# ------------------------------------------------------------------------
head(count_matrix[order(-rowSums(count_matrix)), ], 20)

write.csv(count_matrix, "prjeb22205_counts.csv")

cat("Pipeline complete.\n")
```

---

## Model response

## Scientific/Computational Errors Identified

**1. Reference genome/annotation species mismatch**
```r
genome_index <- "/refs/hg38_star_index"
gtf_file     <- "/refs/Homo_sapiens.GRCh38.gtf"
```
The study (PRJEB22205) is explicitly *Drosophila melanogaster*, but the script aligns reads to the human genome (hg38) and quantifies against a human GTF. Aligning fly reads to a human reference will produce a mapping-rate collapse and essentially meaningless (near-zero, noise-level) counts, since the sequences aren't homologous at the read-mapping level.

*Correction:* use a *D. melanogaster* reference and annotation, e.g. BDGP6/dm6 STAR index and `Drosophila_melanogaster.BDGP6.gtf`.

---

**2. PCR deduplication of RNA-seq reads before quantification**
```r
cmd <- sprintf("picard MarkDuplicates I=%s O=%s REMOVE_DUPLICATES=true M=metrics.txt", ...)
```
Removing duplicate reads is standard for DNA-seq/variant calling, but is **not recommended for standard bulk RNA-seq quantification**. In RNA-seq, reads with identical start/end coordinates frequently arise from genuinely high expression of a transcript (limited fragmentation start sites, high coverage), not solely PCR amplification. Removing them systematically deflates counts for the most highly expressed genes and distorts differential expression results.

*Correction:* skip duplicate removal for standard (non-UMI) bulk RNA-seq quantification; if duplicate marking is truly desired, use it only as a QC metric, not as a filtering step before `featureCounts`, or use UMI-based deduplication if UMIs were part of the library design.

---

**3. Ignoring library strandedness in quantification**
The sample sheet explicitly carries a `library_strandedness` column, but `featureCounts()` is called without a `strandSpecific` argument:
```r
fc <- featureCounts(files = bam_files, annot.ext = gtf_file,
                     isGTFAnnotationFile = TRUE, isPairedEnd = TRUE, nthreads = 8)
```
`featureCounts` defaults to `strandSpecific = 0` (unstranded). If the libraries are actually stranded (e.g., dUTP/Illumina TruSeq stranded protocol), counting as unstranded can substantially misassign or discard reads that overlap genes on the opposite strand, especially at loci with overlapping antisense transcription.

*Correction:* read `samples$library_strandedness` and pass the correct `strandSpecific` value (0/1/2) per the protocol used, run per-strandedness-group if the sheet mixes protocols, and never assume unstranded by default.

---

**4. Sample/column mislabeling due to unguaranteed ordering**
```r
bam_files <- list.files("dedup", pattern = "_dedup.bam$", full.names = TRUE)
fc <- featureCounts(files = bam_files, ...)
count_matrix <- fc$counts
colnames(count_matrix) <- samples$tissue
```
`list.files()` returns files in lexicographic (alphabetical) order, which is not guaranteed to match the row order of `samples` (i.e., the order in `sample_sheet.csv`). Column names (`samples$tissue`) are assigned positionally, so if the alphabetical BAM order differs from the sample sheet order, tissue labels will be silently attached to the wrong count columns — a serious, hard-to-detect data-integrity bug.

*Correction:* explicitly match on run accession, e.g.:
```r
run_ids <- sub("_dedup\\.bam$", "", basename(bam_files))
colnames(count_matrix) <- samples$tissue[match(run_ids, samples$run_accession)]
```

---

**5. Cross-sample comparison performed on raw, unnormalized counts**
```r
head(count_matrix[order(-rowSums(count_matrix)), ], 20)
```
Raw featureCounts output reflects both true expression and technical factors like sequencing depth and (for cross-gene comparisons) gene/transcript length. Ranking or comparing genes/tissues directly on raw counts confounds library-size differences between samples with real biological differences in expression.

*Correction:* normalize before comparison — e.g., CPM/TPM, or use `DESeq2`/`edgeR` size-factor normalization — and perform the "compare expression across tissues" step on normalized values, not raw counts.

---

## Human notes 

- Errors identified: E1,E4,E5,E7,E8
- Missed errors: E2,E3,E6
- Incorrect flags:0