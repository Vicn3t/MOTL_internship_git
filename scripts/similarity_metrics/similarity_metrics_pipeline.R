source("scripts/similarity_metrics/similarity_metrics_functions.R")

library(here)
library(MOFA2)
library(PCAtools)
library(tidyverse)
## -------------------- FORMAT LEARNING DATA --------------------
#
# # Path to learning directory
# LrnFctrnDir <- here("data/LrnDir/Lrn_10000D")
#
# LrnFctrnRdsObjectsDir <- here(LrnFctrnDir, "RDS_objects")
#
# if (!dir.exists(LrnFctrnRdsObjectsDir)) {
#   dir.create(LrnFctrnRdsObjectsDir, showWarnings = FALSE, recursive = TRUE)
# }
#
# # Load metadata (features, samples, GeoMeans, etc.)
# expdat_meta_Lrn <- readRDS(file.path(LrnFctrnDir, "expdat_meta.rds"))
#
# ## Import learning model
# InputModel = file.path(LrnFctrnDir, "Model.hdf5")
# Fctrzn = load_model(file = InputModel)
#
# # Load matrices
# expdat_mRNA_Lrn <- as.matrix(read.csv(
#   file.path(LrnFctrnDir, "mRNA.csv"),
#   header = FALSE
# ))
# rownames(expdat_mRNA_Lrn) <- expdat_meta_Lrn$ftrs_mRNA
# colnames(expdat_mRNA_Lrn) <- expdat_meta_Lrn$smpls
#
#
# expdat_miRNA_Lrn <- as.matrix(read.csv(
#   file.path(LrnFctrnDir, "miRNA.csv"),
#   header = FALSE
# ))
# rownames(expdat_miRNA_Lrn) <- expdat_meta_Lrn$ftrs_miRNA
# colnames(expdat_miRNA_Lrn) <- expdat_meta_Lrn$smpls
#
# expdat_DNAme_Lrn <- as.matrix(read.csv(
#   file.path(LrnFctrnDir, "DNAme.csv"),
#   header = FALSE
# ))
# rownames(expdat_DNAme_Lrn) <- expdat_meta_Lrn$ftrs_DNAme
# colnames(expdat_DNAme_Lrn) <- expdat_meta_Lrn$smpls
#
# expdat_SNV_Lrn <- as.matrix(read.csv(
#   file.path(LrnFctrnDir, "SNV.csv"),
#   header = FALSE
# ))
# rownames(expdat_SNV_Lrn) <- expdat_meta_Lrn$ftrs_SNV
# colnames(expdat_SNV_Lrn) <- expdat_meta_Lrn$smpls
#
#
# # Store in a list
# raw_YLrn_list <- list(
#   mRNA = expdat_mRNA_Lrn,
#   miRNA = expdat_miRNA_Lrn,
#   DNAme = expdat_DNAme_Lrn,
#   SNV = expdat_SNV_Lrn
# )
#
# saveRDS(raw_YLrn_list, here(LrnFctrnRdsObjectsDir, "full_YLrn_list.RDS"))
#
# saveRDS(expdat_mRNA_Lrn, here(LrnFctrnRdsObjectsDir, "expdat_mRNA_Lrn.RDS"))
#
#
# set.seed(123)
# sample <- sort(sample(1:length(colnames(expdat_mRNA_Lrn)), 100))
#
# sample_YLrn_list <- list(
#   mRNA = expdat_mRNA_Lrn[, sample],
#   miRNA = expdat_miRNA_Lrn[, sample],
#   DNAme = expdat_DNAme_Lrn[, sample],
#   SNV = expdat_SNV_Lrn[, sample]
# )
#
# saveRDS(sample_YLrn_list, here(LrnFctrnRdsObjectsDir, "sample_YLrn_list.RDS"))
#
# -------------------- LOAD LEARNING DATA ------------------------------------------------------

barcodes <- readRDS(here("data/barcodes.RDS"))

LrnFctrnDir <- here("data/LrnDir/Lrn_10000D")

LrnFctrnRdsObjectsDir <- here(LrnFctrnDir, "RDS_objects")


# YLrn_list <- readRDS(here(LrnFctrnRdsObjectsDir, "full_YLrn_list.RDS"))

YLrn_list <- readRDS(here(LrnFctrnRdsObjectsDir, "sample_YLrn_list.RDS"))

# -------------------- LOAD TARGET DATA --------------------

# Path to target directory
TrgtDir <- here("data/TrgtDir")

TrgtFctrnDir <- file.path(TrgtDir, "Trg_PAAD_SKCM_SS9_5000D")

TrgtFctrnRdsObjectsDir <- here(TrgtFctrnDir, "RDS_objects")

YTrg_list <- readRDS(file.path(TrgtFctrnRdsObjectsDir, "raw_YTrg_list.RDS"))

views = names(YTrg_list)

# DATASETS PCA ----------------------------------------------------
##  You do NOT have raw counts
##   Objects expected:
##     - YLrn_list$mRNA and YTrg_list$mRNA already normalized on the same scale
## ------------------------------------------------------------

# ---- Single-modality helper ----
pca_compare <- function(
  YLrn, # learning matrix (features x samples)
  YTrg, # target matrix   (features x samples)
  center_L = TRUE,
  scale_L = FALSE,
  center_T = TRUE,
  scale_T = TRUE,
  center_joint = TRUE,
  scale_joint = TRUE,
  removeVar = 0, # pass-through to PCAtools::pca
  label_points = TRUE, # show sample labels on biplots
  point_size = 1.2,
  return_plots = TRUE
) {
  # --- Basic checks ---
  stopifnot(is.matrix(YLrn) || is.data.frame(YLrn))
  stopifnot(is.matrix(YTrg) || is.data.frame(YTrg))
  YLrn <- as.matrix(YLrn)
  YTrg <- as.matrix(YTrg)

  # Ensure rownames exist
  if (is.null(rownames(YLrn)) || is.null(rownames(YTrg))) {
    stop("Both YLrn and YTrg must have rownames (feature IDs) to align on.")
  }

  # 1) Align features (genes)
  common_genes <- intersect(rownames(YLrn), rownames(YTrg))
  if (length(common_genes) < 3) {
    stop("Not enough common features after intersection (need >= 3).")
  }

  YL <- YLrn[common_genes, , drop = FALSE]
  YT <- YTrg[common_genes, , drop = FALSE]

  # --- helpers ---
  .row_impute_median <- function(mat) {
    # Replace NA/NaN/Inf by per-row median (robust, simple)
    mat[!is.finite(mat)] <- NA
    if (all(is.na(mat))) {
      stop("All values are NA after cleaning.")
    }
    # impute row-wise
    med <- apply(mat, 1, function(x) median(x, na.rm = TRUE))
    idx_na <- which(is.na(mat), arr.ind = TRUE)
    if (nrow(idx_na) > 0) {
      mat[idx_na] <- med[idx_na[, 1]]
    }
    mat
  }

  .drop_zero_var_rows <- function(mat) {
    v <- apply(mat, 1, function(x) var(x, na.rm = TRUE))
    mat[which(v > 0 & is.finite(v)), , drop = FALSE]
  }
  ####

  # Option A: impute row medians (recommended quick fix for DNAme)
  YL <- .row_impute_median(YL)
  YT <- .row_impute_median(YT)

  # Drop zero-variance rows after imputation
  YL <- .drop_zero_var_rows(YL)
  YT <- .drop_zero_var_rows(YT)

  # Realign rows after potential row drops
  common_genes <- intersect(rownames(YL), rownames(YT))
  YL <- YL[common_genes, , drop = FALSE]
  YT <- YT[common_genes, , drop = FALSE]

  # 2) Separate PCAs
  pca_L <- PCAtools::pca(
    YL,
    center = center_L,
    scale = scale_L,
    removeVar = removeVar
  )
  pca_T <- PCAtools::pca(
    YT,
    center = center_T,
    scale = scale_T,
    removeVar = removeVar
  )

  # 3) Joint PCA
  data_joint <- cbind(YL, YT)

  pca_J <- PCAtools::pca(
    data_joint,
    center = center_joint,
    scale = scale_joint,
    removeVar = removeVar
  )

  grp <- c(rep("Lrn", ncol(YL)), rep("Trg", ncol(YT)))
  cancer_type <- ifelse(
    grp == "Trg",
    substring(colnames(data_joint), 14, 15), # extract sample type code
    NA
  )

  pca_J$metadata <- data.frame(Group = grp, row.names = colnames(data_joint))

  # 4) Plots (optional)
  plots <- list()
  if (return_plots) {
    plots$learn <- PCAtools::biplot(
      pca_L,
      lab = if (label_points) NULL else NA,
      pointSize = point_size,
      title = "Learning PCA"
    )
    plots$target <- PCAtools::biplot(
      pca_T,
      lab = if (label_points) substring(colnames(YT), 14, 15) else NA,
      pointSize = point_size,
      title = "Target PCA"
    )
    plots$joint <- PCAtools::biplot(
      pca_J,
      colby = "Group",
      lab = if (label_points) cancer_type else NA,
      pointSize = point_size,
      title = "Joint PCA (Lrn vs Trg)"
    )
  }

  # 5) Return everything tidy
  list(
    common_genes = common_genes,
    YL = YL,
    YT = YT,
    pca_learning = pca_L,
    pca_target = pca_T,
    pca_joint = pca_J,
    plots = plots
  )
}

# ---- Multi-modality wrapper (named list of matrices) ----
pca_compare_multi <- function(
  YLrn_list, # named list: e.g., list(mRNA=..., miRNA=...)
  YTrg_list, # named list with same names
  ... # forwarded to pca_compare()
) {
  stopifnot(is.list(YLrn_list), is.list(YTrg_list))
  mods <- intersect(names(YLrn_list), names(YTrg_list))
  if (length(mods) == 0) {
    stop("No common modalities between lists.")
  }

  out <- lapply(mods, function(m) {
    pca_compare(YLrn_list[[m]], YTrg_list[[m]], ...)
  })
  names(out) <- mods
  out
}

pca_compare_multi()
