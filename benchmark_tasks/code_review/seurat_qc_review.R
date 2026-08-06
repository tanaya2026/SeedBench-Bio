#!/usr/bin/env Rscript
# ==============================================================================
# seurat_qc_review.R
#
# Perform quality control on this scRNA-seq mouse cerebellar dataset for downstream integration
#
# NOTE FOR BENCHMARK USERS: this script contains INTENTIONALLY SEEDED ERRORS.
# Do not "fix" this file directly — it is a fixed stimulus for SeedBench-Bio.
# See ground_truth/seurat_qc_review.json for the answer key.
# ==============================================================================

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