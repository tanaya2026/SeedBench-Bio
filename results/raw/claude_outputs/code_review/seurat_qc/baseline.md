# SeedBench-Bio Model Response

## Metadata

- **Prompt variant:** baseline
- **Prompt file:** prompts/code/baseline.md
- **Model:**
provider: Anthropic
name: Claude Sonnet
thinking: medium
- **Task file:** benchmark_tasks/code_review/seurat_qc_review.R
- **Ground truth file:** ground_truth/code_review/seurat_qc_review.json
- **Timestamp:** 2026-08-05 17:510:00
- **Temperature:** 0

---
## Prompt sent to model

Review the following bioinformatics code.

Identify any objective scientific, statistical, or computational errors.

For each error, explain why it is incorrect and suggest a correction.

Ignore coding style, formatting, and personal preferences.


QC_Seurat <- function(seurat_obj, sample_id_metadata_col, UMAP_path = NULL, out_path = NULL){  

  # Load required libraries
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(DoubletFinder))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(biomaRt))
  suppressPackageStartupMessages(library(dplyr))
  
  # Set Default Assay as RNA
  DefaultAssay(seurat_obj) <- "RNA"
  
  
  # Initial filtering
  message("Removing features expressed in 3 or fewer cells")
  # ── 1. Remove features expressed in 3 or fewer cells ──────────────────────────
  counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  genes_keep <- rowSums(counts > 0) > 3
  seurat_obj <- seurat_obj[genes_keep, ]
  message(sprintf("Genes removed (expressed in <=3 cells): %d\n", sum(!genes_keep)))
  message(sprintf("Genes remaining: %d\n", sum(genes_keep)))
  
  message("Removing samples with 5 or fewer cells")
  # ── 2. Remove samples with 5 or fewer cells ───────────────────────────────────
  sample_counts   <- table(seurat_obj@meta.data[, sample_id_metadata_col])
  samples_keep    <- names(sample_counts[sample_counts > 5])
  samples_removed <- names(sample_counts[sample_counts <= 5])
  message(sprintf("Samples removed (<=5 cells): %s\n",
                  ifelse(length(samples_removed) > 0, paste(samples_removed, collapse = ", "), "none")))
  message(sprintf("Samples remaining: %d\n", length(samples_keep)))
  seurat_obj <- seurat_obj[, seurat_obj@meta.data[[sample_id_metadata_col]] %in% samples_keep]
  
  
  message("Calculating QC metrics")
  # ── 3. Calculate QC metrics ───────────────────────────────────────────────────
  
  seurat_obj[["percent.mt"]]          <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rb"]]          <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]")
  seurat_obj[["log10GenesPerUMI"]]    <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  
  # ── 4. Compute adaptive thresholds (4 SD, per Vladoiu 2019) ──────────────────
  mt_limit         <- median(seurat_obj$percent.mt)          + (4 * sd(seurat_obj$percent.mt))
  feat_limit       <- median(seurat_obj$nFeature_RNA)        + (4 * sd(seurat_obj$nFeature_RNA))
  count_limit      <- median(seurat_obj$nCount_RNA)          + (4 * sd(seurat_obj$nCount_RNA))
  rb_limit         <- median(seurat_obj$percent.rb)          + (4 * sd(seurat_obj$percent.rb))
  complexity_limit <- median(seurat_obj$log10GenesPerUMI)    - (4 * sd(seurat_obj$log10GenesPerUMI))
  
  message(sprintf("Mitochondrial limit:        %.4f%%\n", mt_limit))
  message(sprintf("Feature upper limit:        %.0f\n",   feat_limit))
  message(sprintf("UMI count upper limit:      %.0f\n",   count_limit))
  message(sprintf("Ribosomal upper limit:      %.4f%%\n", rb_limit))
  message(sprintf("Complexity lower limit:     %.4f\n",   complexity_limit))
  
  # ── 5. Filter cells ───────────────────────────────────────────────────────────
  cells_before <- ncol(seurat_obj)
  seurat_obj <- seurat_obj[,
                           seurat_obj[["nFeature_RNA"]]     >= 200             &
                             seurat_obj[["nFeature_RNA"]]     <= feat_limit      &
                             seurat_obj[["nCount_RNA"]]       <= count_limit     &
                             seurat_obj[["percent.mt"]]       <= mt_limit        &
                             seurat_obj[["percent.rb"]]       <= rb_limit        &
                             seurat_obj[["log10GenesPerUMI"]] >= complexity_limit
  ]
  cells_after <- ncol(seurat_obj)
  message(sprintf("Cells before QC filtering:  %d\n",           cells_before))
  message(sprintf("Cells after QC filtering:   %d\n",           cells_after))
  message(sprintf("Cells removed:              %d (%.1f%%)\n",
                  cells_before - cells_after,
                  100 * (cells_before - cells_after) / cells_before))
  
  # ── 6. Doublet Detection ──────────────────────────────────────────────────────
  message("Starting doublet detection (fixed pK = 0.25, skipping paramSweep)...")
  
  options(future.globals.maxSize = Inf)  
  
  seurat_obj <- SCTransform(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.2, cluster.name = "unintegrated_0.2")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "pca", reduction.name = "unintegrated")
  
  samp_split <- SplitObject(seurat_obj, split.by = sample_id_metadata_col)
  
  fixed_pk <- 0.25
  
  for (samp_name in names(samp_split)) {
    
    message("\n", strrep("=", 60), "\n")
    message(sprintf("Processing sample: %s\n", samp_name))
    samp <- samp_split[[samp_name]]
    message(sprintf("  Cells: %d | Genes: %d\n", ncol(samp), nrow(samp)))
    message(strrep("=", 60), "\n")
    
    multiplet_rate <- 0.076
    nExp.poi       <- round(multiplet_rate * nrow(samp))
    nExp.poi.adj   <- round(nExp.poi * (1 - multiplet_rate))
    message(sprintf("  Multiplet rate: %.3f\n", multiplet_rate))
    message(sprintf("  Fixed pK:%.3f\n", fixed_pk))
    message(sprintf("  Expected doublets (adjusted): %d\n", nExp.poi.adj))
    
    # ── DoubletFinder ────────────────────────────────────────────
    message("\n[5/5] Running DoubletFinder...\n"); flush.console()
    samp <- suppressMessages(doubletFinder(seu  = samp,
                                           PCs  = 1:30,
                                           pK   = fixed_pk,
                                           nExp = nExp.poi.adj,
                                           sct  = TRUE))
    colnames(samp@meta.data)[grepl("DF.classifications.*", colnames(samp@meta.data))] <- "doublet_finder"
    
    doublet_table <- table(samp@meta.data$doublet_finder)
    message(sprintf("\n  Singlets:  %d\n",    doublet_table["Singlet"]))
    message(sprintf("  Doublets:  %d\n",      doublet_table["Doublet"]))
    message(sprintf("  Doublet %%: %.2f%%\n", 100 * doublet_table["Doublet"] / ncol(samp)))
    message(sprintf("\nSample %s complete\n",  samp_name)); flush.console()
    
    samp_split[[samp_name]] <- samp
  }
  
  message("\n", strrep("=", 60), "\n")
  message("DoubletFinder complete for all samples\n")
  message(strrep("=", 60), "\n")
  
  seurat_postQC <- merge(x = samp_split[[1]], y = samp_split[2:length(samp_split)], merge.dr = TRUE)
  
  if (!is.null(UMAP_path)) {
    message(paste("Printing doublet UMAP to path:", UMAP_path))
    Idents(seurat_postQC) <- "doublet_finder"
    p <- Seurat::DimPlot(seurat_postQC)
    ggsave(filename = UMAP_path, plot = p)
  }
  
  seu_dblt <- subset(seurat_postQC, doublet_finder == "Singlet")
  
  
  total_cells  <- ncol(seurat_postQC)
  singlet_cells <- ncol(seu_dblt)
  doublet_cells <- total_cells - singlet_cells
  message(sprintf("Doublets removed:  %d (%.1f%% of post-QC cells)\n",
                  doublet_cells,
                  100 * doublet_cells / total_cells))
  message(sprintf("Cells retained:    %d\n", singlet_cells))
  
  DefaultAssay(seu_dblt) <- "RNA"
  seu_dblt <- JoinLayers(seu_dblt)
  
  # Quick normalization solely for CellCycleScoring — SCTransform will replace this
  seu_dblt <- NormalizeData(seu_dblt)
  
  # ── 7. Cell Cycle Scoring and SCTransform ─────────────────────────────────────
  message("Scoring cell cycle and regressing CC Difference...")
  
  # Extract cell cycle genes
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  
  srat <- CellCycleScoring(
    seu_dblt,
    s.features   = s.genes,
    g2m.features = g2m.genes
  )

  message("SCTransforming")
  srat <- SCTransform(srat, vars.to.regress = "CC.Difference")

  srat$CC.Difference <- srat$S.Score - srat$G2M.Score
  
  srat <- RunPCA(srat,        assay = "SCT", npcs = 100)
  srat <- RunUMAP(srat,       assay = "SCT", dims = 1:50)
  srat <- FindNeighbors(srat, reduction = "pca", dims = 1:50)
  srat <- FindClusters(srat,  resolution = 0.2)
  
  message("Object ready for plotting and integration")
  
  if (!is.null(out_path)){
    message(sprintf("\n***Saving processed Seurat object to %s***\n", out_path))
    saveRDS(object = srat, file = out_path)
  }
  
  return(srat)
}

#______________________________________________________________________________
# Read in dataset 
(library(Seurat))
# Mouse Cerebellum Dataset
sample_obj <-readRDS("sample_seurat.rds")

# Call function
post_qc <- QC_Seurat(sample_obj, sample_id_metadata_col = "age", UMAP_path = "~/sample_UMAP.png", out_path = "~/sample_postQC.rds")

---

## Model response

Here's a review focused on substantive errors — I'm skipping style/formatting issues as requested.

## 1. Species mismatch: mitochondrial gene regex is wrong for this dataset
```r
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
```
The script is explicitly run on a **mouse cerebellum dataset** (`sample_id_metadata_col = "age"`, comment "Mouse Cerebellum Dataset"). Mouse mitochondrial genes are annotated as `mt-Nd1`, `mt-Co1`, etc. (only the "m" capitalized), not `MT-ND1` (human convention). Since the pattern match is case-sensitive, `"^MT-"` will match **zero genes** in a mouse object, so `percent.mt` will be 0 for every cell and the mitochondrial filter (`mt_limit`) will do nothing. This should be `pattern = "^mt-"`.

## 2. Species mismatch: cell cycle genes are human symbols
```r
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
srat <- CellCycleScoring(seu_dblt, s.features = s.genes, g2m.features = g2m.genes)
```
`cc.genes.updated.2019` is Seurat's built-in **human** gene symbol list (all-caps, e.g. `"MCM5"`, `"PCNA"`). Against mouse features (`"Mcm5"`, `"Pcna"`), essentially none of these will match, so `S.Score`/`G2M.Score` will be computed from little or no signal and cell-cycle regression will be meaningless. Notably, `biomaRt` is loaded at the top of the function but never actually called — this strongly suggests an intended human→mouse ortholog conversion step is missing. This needs to be added (e.g., via `biomaRt`, or by using a mouse-specific cell cycle gene list) before scoring.

## 3. Expected doublet count uses gene count, not cell count
```r
nExp.poi <- round(multiplet_rate * nrow(samp))
```
In a Seurat object, `nrow()` returns the number of **features/genes**, and `ncol()` returns the number of **cells** (the code itself gets this right two lines earlier: `"Cells: %d | Genes: %d", ncol(samp), nrow(samp)`). The multiplet rate is a fraction of *cells*, so this must be `ncol(samp)`. As written, `nExp.poi` will be based on the (typically far larger) gene count, producing wildly inflated expected-doublet numbers and effectively meaningless doublet calls.

## 4. Homotypic doublet proportion adjustment is conceptually wrong
```r
nExp.poi.adj <- round(nExp.poi * (1 - multiplet_rate))
```
The standard DoubletFinder workflow adjusts `nExp.poi` for the fraction of doublets that are "homotypic" (same-cluster, and thus undetectable) using `modelHomotypic()` applied to per-cell cluster annotations:
```r
homotypic.prop <- modelHomotypic(annotations)
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
```
Reusing `multiplet_rate` (0.076, the *overall doublet rate*) in place of `homotypic.prop` conflates two different quantities and has no statistical justification — it isn't derived from cluster composition at all.

## 5. Likely pK/pN parameter confusion, and paramSweep is skipped entirely
```r
fixed_pk <- 0.25
...
samp <- doubletFinder(seu = samp, PCs = 1:30, pK = fixed_pk, nExp = nExp.poi.adj, sct = TRUE)
```
Two issues here:
- DoubletFinder's `pN` (proportion of artificial doublets) has a conventional, largely insensitive default of **0.25**. `pK` (the neighborhood size parameter) is dataset-specific and normally *tuned per sample* via `paramSweep()`/`find.pK()`, typically landing in the 0.005–0.03 range, rarely near 0.25. Hard-coding `pK = 0.25` looks like a mix-up between these two parameters.
- The comment even states paramSweep is intentionally skipped. Without it, pK is essentially arbitrary rather than optimized for each sample's data structure, which materially affects doublet classification accuracy.

## 6. Fixed multiplet rate applied identically across all samples
```r
multiplet_rate <- 0.076
```
This is applied uniformly inside the per-sample loop. 10x Genomics' expected doublet rate scales with the number of cells loaded/recovered (roughly linearly, e.g. ~0.4–0.8% per 1,000 cells), so if samples differ in cell count — which is the whole reason for per-sample splitting — they should not share one fixed rate. Each sample's rate should be computed from its own `ncol(samp)` using the standard 10x rate table (or an equivalent formula).

## 7. Order-of-operations bug: `CC.Difference` doesn't exist when SCTransform uses it
```r
srat <- SCTransform(srat, vars.to.regress = "CC.Difference")
srat$CC.Difference <- srat$S.Score - srat$G2M.Score
```
`CellCycleScoring()` adds `S.Score`, `G2M.Score`, and `Phase` to metadata, but **not** `CC.Difference` — that's computed manually on the very next line, *after* `SCTransform` is called asking to regress on it. As written, `vars.to.regress = "CC.Difference"` references a column that doesn't exist yet at call time, which will either error out or silently fail to regress cell-cycle effects. The `CC.Difference` assignment line needs to be moved above the `SCTransform` call.

## 8. Mixing robust and non-robust statistics for adaptive thresholds
```r
mt_limit <- median(seurat_obj$percent.mt) + (4 * sd(seurat_obj$percent.mt))
...
complexity_limit <- median(seurat_obj$log10GenesPerUMI) - (4 * sd(seurat_obj$log10GenesPerUMI))
```
`median` is a robust estimator of center, but `sd` is not robust to outliers/skew — and scRNA-seq QC metrics (UMI counts, feature counts, %MT) are typically strongly right-skewed. Combining median with SD is statistically inconsistent and the SD term can be disproportionately inflated by the very outlier cells the filter is meant to exclude. A more standard/robust approach (e.g., used by `scran::isOutlier`) is median ± *n* × MAD (median absolute deviation), not SD.

## 9. Minor: pulling "counts" from the `data` layer
```r
counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
genes_keep <- rowSums(counts > 0) > 3
```
This is labeled `counts` but pulled from the `data` layer. If the input object hasn't been normalized yet, `data` and `counts` are identical, so this happens to work here — but it's fragile: if the input object were ever pre-normalized, this would silently operate on log-normalized values instead of raw counts (though the >0 threshold behavior would still coincidentally hold, since log1p preserves zero/nonzero status). Explicitly pulling `layer = "counts"` would remove the ambiguity.

---

**Most critical to fix first:** #1 and #2 (species mismatch — these silently produce near-meaningless QC/CC results on mouse data), #3 (cell vs. gene count bug — corrupts doublet detection), and #7 (order bug that likely throws an error or silently skips CC regression).

---

## Human notes 

- Errors identified: E1-E5
- Missed errors: 0
- Incorrect flags: 0