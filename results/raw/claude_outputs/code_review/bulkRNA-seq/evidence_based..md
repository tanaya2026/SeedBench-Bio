# SeedBench-Bio Model Response

## Metadata

- **Prompt variant:** evidence_based
- **Prompt file:** prompts/code/evidence_based.md
- **Model:**
provider: Anthropic
name: Claude Sonnet
thinking: medium
- **Task file:** benchmark_tasks/code_review/bulk_RNAseq_matrix.R
- **Ground truth file:** ground_truth/code_review/bulk_RNAseq_matrix.json
- **Timestamp:** 2026-08-06 12:27:00
- **Temperature:** 0

---
## Prompt sent to model
Review the following code.

Report only objective errors supported by evidence from the code.

For every reported issue provide:

1. The relevant code location.
2. Why the code is objectively incorrect.
3. The bioinformatics, statistical, or computational principle that is violated.
4. A corrected implementation.

If you cannot justify an issue with evidence from the code, do not report it.

Ignore style, formatting, efficiency, and naming conventions.

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

---

## Model response
## Issue 1: PCR‑duplicate removal is applied to RNA‑seq quantification data

**Location:**
```r
cmd <- sprintf("picard MarkDuplicates I=%s O=%s REMOVE_DUPLICATES=true M=metrics.txt",
               bam_in, bam_out)
```
followed by using the `dedup` BAMs as input to `featureCounts`.

**Why incorrect:** Picard `MarkDuplicates` flags/removes reads as duplicates solely based on identical 5′ mapping coordinates. In RNA‑seq, reads from highly expressed transcripts legitimately pile up at the same fragment start positions because of high biological abundance, not PCR amplification bias — there is no molecular tag (e.g., UMI) here to distinguish PCR duplicates from independently sampled fragments of an abundant transcript. Removing these reads preferentially strips signal from highly expressed genes and systematically distorts the count matrix used for expression comparisons.

**Principle violated:** Standard RNA‑seq quantification practice (e.g., ENCODE/GATK/Picard guidance) explicitly recommends *against* position‑based duplicate removal for standard (non‑UMI) RNA‑seq gene expression quantification, because expression level and duplication rate are confounded — this is a statistical/methodological error, not a matter of style.

**Correction:** Quantify directly from the STAR-sorted BAMs without duplicate removal (or, if duplicate marking is desired for QC, mark but do not remove, and never feed the "removed" BAM into featureCounts):

```r
bam_files <- file.path("star_out", paste0(samples$run_accession, "_Aligned.sortedByCoord.out.bam"))

fc <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,
  nthreads = 8
)
```

---

## Issue 2: `library_strandedness` is read but never used in quantification

**Location:**
```r
samples <- fread("sample_sheet.csv")   # includes library_strandedness column
...
fc <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,
  nthreads = 8
)
```

**Why incorrect:** `featureCounts` defaults to `strandSpecific = 0` (unstranded) when the argument is omitted. The sample sheet explicitly carries a `library_strandedness` field per sample, indicating that strand information is relevant/variable, yet it is never passed to `featureCounts`. For a stranded library, counting as unstranded assigns reads to any gene overlapping the fragment on either strand, inflating counts for antisense overlapping genes and misassigning reads in loci with overlapping genes on opposite strands.

**Principle violated:** Read counting must match the actual library protocol; using the wrong strandedness setting for a stranded protocol produces biased/incorrect gene-level counts (a computational‑correctness issue, not merely suboptimal).

**Correction:** Map the sample sheet's strandedness values to featureCounts' `strandSpecific` codes and pass them (per unique value, since featureCounts takes one setting per call or a vector matching `files`):

```r
strand_map <- c(unstranded = 0, forward = 1, reverse = 2)
strand_codes <- strand_map[samples$library_strandedness]

fc <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,
  strandSpecific = strand_codes,
  nthreads = 8
)
```

---

## Issue 3: Sample labels are assigned without matching to the actual column order of the count matrix

**Location:**
```r
bam_files <- list.files("dedup", pattern = "_dedup.bam$", full.names = TRUE)
fc <- featureCounts(files = bam_files, ...)
count_matrix <- fc$counts
colnames(count_matrix) <- samples$tissue
```

**Why incorrect:** `list.files()` returns files in lexicographic (alphabetical) order, which determines the column order of `fc$counts`. `samples$tissue` is in the row order of the original sample sheet (`sample_sheet.csv`). These two orderings are not guaranteed to match unless `run_accession` values happen to sort identically to the sample‑sheet row order. Assigning `samples$tissue` directly onto `count_matrix` columns therefore risks silently mislabeling samples with the wrong tissue.

**Principle violated:** Data provenance/identity integrity — sample metadata must be joined to data by an explicit key, not by positional assumption, when the two are generated by independently ordered processes.

**Correction:** Derive the run accession from each BAM filename and explicitly match it to the sample sheet:

```r
bam_files <- list.files("dedup", pattern = "_dedup.bam$", full.names = TRUE)
bam_accessions <- sub("_dedup\\.bam$", "", basename(bam_files))

fc <- featureCounts(files = bam_files, annot.ext = gtf_file,
                     isGTFAnnotationFile = TRUE, isPairedEnd = TRUE, nthreads = 8)

count_matrix <- fc$counts
match_idx <- match(bam_accessions, samples$run_accession)
colnames(count_matrix) <- samples$tissue[match_idx]
```

---

## Human notes 

- Errors identified: E4, E5, E7
- Missed errors: E1, E2, E3, E6, E8
- Incorrect flags: 0