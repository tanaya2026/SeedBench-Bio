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