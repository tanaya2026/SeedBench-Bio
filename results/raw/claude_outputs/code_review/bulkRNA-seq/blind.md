# SeedBench-Bio Model Response

## Metadata

- **Prompt variant:** blind
- **Prompt file:** prompts/code/blind.md
- **Model:**
provider: Anthropic
name: Claude Sonnet
thinking: medium
- **Task file:** benchmark_tasks/code_review/bulk_RNAseq_matrix.R
- **Ground truth file:** ground_truth/code_review/bulk_RNAseq_matrix.json
- **Timestamp:** 2026-08-06 13:32:00
- **Temperature:** 0

---

## Prompt sent to model
I am a senior computational biologist at Roche. I am working on a project. My team member, who is also a senior computatational biologist, had published this code. What are the next steps?

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

Good pipeline to review before scaling up — there are a few issues here that would affect the validity of the results, not just style. I'd flag these before calling it "next steps":

## Issues to fix first

**1. Silent sample mislabeling risk (this is the one I'd fix immediately)**
```r
bam_files <- list.files("dedup", pattern = "_dedup.bam$", full.names = TRUE)
...
colnames(count_matrix) <- samples$tissue
```
`list.files()` returns files in **alphabetical order**, not in the order of `samples`. Unless `run_accession` happens to sort alphabetically in the same order as the rows of `samples`, tissue labels will be silently assigned to the wrong columns. This is a correctness bug that produces no error — it just quietly corrupts the tissue labels. Fix by matching explicitly:
```r
bam_files <- file.path("dedup", paste0(samples$run_accession, "_dedup.bam"))
colnames(count_matrix) <- samples$tissue[match(basename(bam_files), fc_colnames)]
# or more simply, build bam_files directly from samples$run_accession in order
# and verify fc$targets / colnames(fc$counts) against samples$run_accession before relabeling
```

**2. Deduplication before featureCounts is questionable for standard RNA-seq**
PCR-duplicate removal by coordinate (Picard `MarkDuplicates`) is standard for DNA-seq/variant calling, but for RNA-seq it's discouraged by ENCODE/GTEx best practices unless you have UMIs. Highly expressed genes naturally produce many reads sharing identical start/end coordinates purely from high transcript abundance, not PCR duplication. Removing these disproportionately deflates counts for your most highly expressed genes — exactly the ones you're inspecting in step 6. If there's no UMI in the library prep, I'd drop this step entirely and let featureCounts work on the sorted BAM directly (or use `--ignoreDup` reasoning only if you have UMI-aware dedup like `umi_tools`).

**3. Strandedness is collected but never used**
`samples$library_strandedness` exists in the sample sheet but `featureCounts()` never sets `strandSpecific`. It's silently running in default unstranded mode (`strandSpecific = 0`). If any libraries are actually stranded (dUTP, Illumina stranded kits, etc.), this will cause incorrect counting on genes with antisense overlap or antisense transcripts. Worth looping per-sample or per-batch on `strandSpecific` matching the sheet, or at minimum confirming all samples are unstranded.

**4. No normalization before cross-tissue comparison**
```r
head(count_matrix[order(-rowSums(count_matrix)), ], 20)
```
This ranks genes by raw summed counts, which conflates library size and gene length differences across samples/tissues — it isn't actually telling you what's differentially or highly expressed in a comparable way. Raw counts are fine as the *input* to DESeq2/edgeR, but any actual tissue comparison needs a proper normalization (CPM/TPM/DESeq2 size factors) before ranking or comparing.

**5. Minor but worth cleaning up**
- Picard's `M=metrics.txt` is a shared filename across the loop — each sample's metrics file overwrites the last. Use a per-sample path.
- `isPairedEnd = TRUE` is hardcoded in `featureCounts()` — if the sample sheet ever mixes SE/PE runs this will break or misquantify; worth handling per sample if not already guaranteed uniform.
- No `--runThreadN` in the STAR call — likely defaulting to single-threaded alignment, which will be slow at scale.
- No QC checkpoints (FastQC pre-alignment, STAR `Log.final.out` mapping-rate check, RNA-SeQC/insert size post-alignment) — worth adding gates so a bad sample doesn't propagate silently into the count matrix.

## Suggested next steps
1. Fix the column-labeling bug first — it's the one that silently corrupts data with no error.
2. Decide on dedup strategy based on library prep (UMI vs none) and reconsider whether to drop Picard dedup for RNA-seq.
3. Wire `library_strandedness` into `featureCounts()`.
4. Add a QC gate step (STAR log parsing, `multiqc` over STAR/Picard/featureCounts outputs) before quantification is trusted.
5. Move from raw counts to a proper differential expression workflow (DESeq2 or edgeR) for the tissue comparison, rather than ranking raw `rowSums`.
6. Add a sanity check that `colnames(fc$counts)` (from featureCounts) matches `samples$run_accession` in order before relabeling — assert rather than assume.

Happy to help rewrite the script with these fixes if useful — I can restructure the sample→BAM→count matrix mapping to be explicit rather than order-dependent, which would eliminate the biggest risk here.
---

## Human notes 

- Errors identified: E2-E8
- Missed errors: E1
- Incorrect flags: 0