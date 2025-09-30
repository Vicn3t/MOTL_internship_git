# libraries and functions sourcing ----------------------------------------
library(SummarizedExperiment)
library(DESeq2)
library(sesame)
library(dplyr)
library(MOFA2)
library(rhdf5)
library(ggplot2)
library(here)
library(stringr)

source('scripts/MOTL/TCGA_preprocessedData_functions.R')

source('scripts/MOTL/TL_VI_functions.R')


TrgtDir <- here("data/TrgtDir")

CenterTrg = FALSE


# Load TGCA Learning data -------------------------------------------------------------

LrnDir = here('data/LrnDir')

LrnFctrnDir = file.path(LrnDir, 'Lrn_5000D_Fctrzn_100K_001TH')
# LrnFctrnDir = file.path(LrnDir, 'Lrn_10000D')

## Import learning data
expdat_meta_Lrn = readRDS(file.path(LrnFctrnDir, "expdat_meta.rds"))
InputModel = file.path(LrnFctrnDir, "Model.hdf5")
Fctrzn = load_model(file = InputModel)

## Initialise values from factorization
viewsLrn = Fctrzn@data_options$views
likelihoodsLrn = Fctrzn@model_options$likelihoods
MLrn = Fctrzn@dimensions$M

Fctrzn@expectations[["Tau"]] = Tau_init(viewsLrn, Fctrzn, InputModel)

Fctrzn@expectations[["TauLn"]] = sapply(
  viewsLrn,
  TauLn_calculation,
  likelihoodsLrn,
  Fctrzn,
  LrnFctrnDir
)
Fctrzn@expectations[["WSq"]] = sapply(
  viewsLrn,
  WSq_calculation,
  Fctrzn,
  LrnFctrnDir
)
Fctrzn@expectations[["W0"]] = sapply(
  viewsLrn,
  W0_calculation,
  CenterTrg,
  Fctrzn,
  LrnFctrnDir
)


# # Create Target Dataset TGCA -----------------------------------------------------
# # TrgtFctrnDir <- file.path(TrgtDir, "Trg_LAML_PAAD_SKCM_SS5_5000D")
#
# TrgtFctrnDir <- file.path(TrgtDir, "Trg_PAAD_SKCM_SS9_5000D")
# expdat_meta_Trgt = readRDS(file.path(TrgtFctrnDir, "expdat_meta.rds"))
#
# ## mRNA
# raw_expdat_mRNA = as.matrix(read.csv(
#   file.path(TrgtFctrnDir, "mRNA.csv"),
#   header = FALSE
# ))
#
# rownames(raw_expdat_mRNA) = word(
#   expdat_meta_Trgt$ftrs_mRNA,
#   1,
#   sep = fixed(".")
# )
#
# colnames(raw_expdat_mRNA) = expdat_meta_Trgt$smpls
#
# expdat_mRNA = mRNA_addVersion(
#   expdat = raw_expdat_mRNA,
#   Lrndat = Fctrzn@expectations$W$mRNA
# )
#
#
# target_views <- c("miRNA", "DNAme", "SNV")
#
# load_target_data <- function(view_name, meta_data, data_dir) {
#   file_path <- file.path(data_dir, paste0(view_name, ".csv"))
#   data <- as.matrix(read.csv(file_path, header = FALSE))
#
#   rownames(data) <- meta_data[[paste0("ftrs_", view_name)]]
#   colnames(data) <- meta_data$smpls
#
#   return(data)
# }
#
# expdat_list <- lapply(target_views, function(v) {
#   load_target_data(v, expdat_meta_Trgt, TrgtFctrnDir)
# })
#
# names(expdat_list) <- target_views
#
# raw_YTrg_list = list(
#   mRNA = expdat_mRNA,
#   miRNA = expdat_list$miRNA,
#   DNAme = expdat_list$DNAme,
#   SNV = expdat_list$SNV
# )
#
#
# TrgtFctrnRdsObjectsDir <- here(TrgtFctrnDir, "RDS_objects")
#
# if (!dir.exists(TrgtFctrnRdsObjectsDir)) {
#   dir.create(TrgtFctrnRdsObjectsDir, showWarnings = FALSE, recursive = TRUE)
# }
#
# saveRDS(raw_YTrg_list, file.path(TrgtFctrnRdsObjectsDir, "raw_YTrg_list.RDS"))
#
# Load Target Dataset TCGA -------------------------------------

TrgtFctrnDir <- file.path(TrgtDir, "Trg_PAAD_SKCM_SS9_5000D")

raw_YTrg_list <- readRDS(file.path(TrgtFctrnDir, "raw_YTrg_list.RDS"))

# raw_YTrg_list <- raw_YTrg_list[c(1, 2)]

# # CreateTarget Dataset ---------------------------------------------------------------
#
# Goal:
#   1) Read raw GEO counts
#   2) Normalize Ensembl IDs (remove version suffix)
#   3) Pad sample IDs (AD1 -> AD01, HC1 -> HC01) and sort columns
#   4) Map Ensembl IDs -> biotype + miRBase ID (and external gene name)
#   5) Split into mRNA (protein_coding) and miRNA (miRNA + mirbase_id)
#   6) Rename miRNA rows by miRBase ID
#
# library(data.table)
# library(biomaRt)
#
# TrgtFctrnDir <- file.path(TrgtDir, "Alzheimer_data")
#
# ## 1) Read raw counts, clean Ensembl IDs, split meta vs counts
# file <- file.path(TrgtFctrnDir, "GSE276756_raw_counts.txt")
# dt <- fread(file)
#
# strip_version <- function(x) sub("[.].*$", "", x)
#
# dt[, Geneid := strip_version(Geneid)]
# setorder(dt, Geneid)
#
# rownames(dt) <- dt$Geneid
# dt = mRNA_addVersion(expdat = dt, Lrndat = Fctrzn@expectations$W$mRNA)
#
# rn <- rownames(dt)
#
# meta_cols <- intersect(
#   c("Geneid", "Chr", "Start", "End", "Strand", "Length"),
#   names(dt)
# )
# count_cols <- setdiff(names(dt), meta_cols)
#
# meta <- as.data.frame(dt[, ..meta_cols], row.names = rn)
# counts <- as.matrix(dt[, ..count_cols])
# rownames(counts) <- rn
# storage.mode(counts) <- "numeric" # fast, safe numeric conversion
#
# # Normalize sample IDs (AD1->AD01, HC1->HC01) and sort columns
# pad_ids <- function(x) {
#   x <- gsub("AD([0-9])$", "AD0\\1", x)
#   x <- gsub("HC([0-9])$", "HC0\\1", x)
#   x
# }
# colnames(counts) <- pad_ids(colnames(counts))
# counts <- counts[, order(colnames(counts)), drop = FALSE]
#
# meta <- meta[apply(counts, 1, var) > 0, ]
# counts <- counts[apply(counts, 1, var) > 0, ]
#
# ## 2) Map Ensembl -> biotype / miRBase / gene name
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#
# map <- getBM(
#   attributes = c("ensembl_gene_id", "gene_biotype", "mirbase_id"),
#   filters = "ensembl_gene_id",
#   values = meta$Geneid,
#   mart = mart
# )
# is_mRNA <- map$ensembl_gene_id[map$gene_biotype == "protein_coding"]
#
# mi_map <- map %>%
#   filter(
#     map$mirbase_id != "" &
#       !is.na(map$mirbase_id) &
#       map$gene_biotype == "miRNA" &
#       !duplicated(map$mirbase_id) &
#       !duplicated(map$ensembl_gene_id)
#   )
#
# mi_map_ord <- mi_map[order(mi_map$ensembl_gene_id), ]
#
#
# ## 3) Split mRNA / miRNA and rename miRNA rows
# is_miRNA <- mi_map_ord$ensembl_gene_id
#
# mRNA_counts <- counts[meta$Geneid %in% is_mRNA, , drop = FALSE]
# miRNA_counts <- counts[meta$Geneid %in% is_miRNA, , drop = FALSE]
# mRNA_meta <- meta[meta$Geneid %in% is_mRNA, , drop = FALSE]
# miRNA_meta <- meta[meta$Geneid %in% is_miRNA, , drop = FALSE]
#
# rownames(miRNA_counts) <- mi_map_ord$mirbase_id
# miRNA_meta$map_Geneid <- mi_map_ord$ensembl_gene_id
# miRNA_meta$map_mirbase_id <- mi_map_ord$mirbase_id
#
# non_0_var_ids <- apply(miRNA_counts, 2, function(x) var(x, na.rm = TRUE) == 0)
#
# mRNA_counts <- mRNA_counts[, !non_0_var_ids]
# miRNA_counts <- miRNA_counts[, !non_0_var_ids]
#
# cat("Dimensions:\n")
# cat(
#   "  mRNA_counts: ",
#   nrow(mRNA_counts),
#   " x ",
#   ncol(mRNA_counts),
#   "\n",
#   sep = ""
# )
# cat(
#   "  miRNA_counts:",
#   nrow(miRNA_counts),
#   " x ",
#   ncol(miRNA_counts),
#   "\n",
#   sep = ""
# )
#
# ## (Optional) 4) Save and reload cleanly (preserve rownames; avoid extra 'X' column)
# write.csv(
#   mRNA_counts,
#   file.path(TrgtFctrnDir, "GSE276756_mRNA_counts.csv"),
#   row.names = TRUE
# )
# write.csv(
#   miRNA_counts,
#   file.path(TrgtFctrnDir, "GSE276756_miRNA_counts.csv"),
#   row.names = TRUE
# )

# # Load Alzheimer matrices -------------------------------------------------
# #
# # To reload:
# TrgtFctrnDir <- file.path(TrgtDir, "Alzheimer_data")
#
# mRNA_counts <- as.matrix(read.csv(
#   file.path(TrgtFctrnDir, "GSE276756_mRNA_counts.csv"),
#   row.names = 1,
#   check.names = FALSE
# ))
# miRNA_counts <- as.matrix(read.csv(
#   file.path(TrgtFctrnDir, "GSE276756_miRNA_counts.csv"),
#   row.names = 1,
#   check.names = FALSE
# ))
# storage.mode(mRNA_counts) <- "integer"
# storage.mode(miRNA_counts) <- "integer"
#
# # Pack into a list for downstream use (requires DNA methylation object if available)
# raw_YTrg_list = list(
#   mRNA = mRNA_counts,
#   miRNA = miRNA_counts
# )
#
# MOTL Factorization -----------------------------------------------------------

smpls = colnames(raw_YTrg_list[[1]])
viewsTrg = names(raw_YTrg_list)
views = viewsLrn[is.element(viewsLrn, viewsTrg)]
likelihoods = likelihoodsLrn[views]

YTrg_list = TargetDataPreparation(
  views = views,
  YTrg_list = raw_YTrg_list,
  Fctrzn = Fctrzn,
  smpls = smpls,
  expdat_meta_Lrn = expdat_meta_Lrn,
  # normalization = 'Lrn',
  # transformation = TRUE
  normalization = FALSE,
  transformation = FALSE
)

TL_param = initTransferLearningParamaters(
  YTrg = YTrg_list,
  views = views,
  expdat_meta_Lrn = expdat_meta_Lrn,
  Fctrzn = Fctrzn,
  likelihoods = likelihoods
)

TL_OutDir = 'MOTL_Fctrzn'
if (!dir.exists(TL_OutDir)) {
  dir.create(TL_OutDir, showWarnings = FALSE, recursive = TRUE)
}

minFactors = 6 ## floor when dropping factors
StartDropFactor = 1 # after which iteration to start dropping factors
FreqDropFactor = 1 # how often to drop factors
StartELBO = 1 # which iteration to start checking ELBO on, excl initiation iteration
FreqELBO = 5 # how often to assess the ELBO
DropFactorTH = 0.01 # factor with lowest max variance, that is less than this, is dropped
MaxIterations = 10000
MinIterations = 2 #
ConvergenceIts = 2 # number of consectutive checks in a row for which the change in elbo is below the threshold
ConvergenceTH = 0.0005 # change in elbo threshold
ss_start_time = Sys.time()

TL_data = transferLearning_function(
  TL_param = TL_param,
  MaxIterations = MaxIterations,
  MinIterations = MinIterations,
  minFactors = minFactors,
  StartDropFactor = StartDropFactor,
  FreqDropFactor = FreqDropFactor,
  StartELBO = StartELBO,
  FreqELBO = FreqELBO,
  DropFactorTH = DropFactorTH,
  ConvergenceIts = ConvergenceIts,
  ConvergenceTH = ConvergenceTH,
  CenterTrg = CenterTrg,
  ss_start_time = ss_start_time,
  outputDir = TL_OutDir
)


ZMu = TL_data$ZMu
ZMu_0 <- TL_data$ZMu_0 # N (intercept facteurs)

## Factorization quality estimation--------------------------------------------------------------
## Function to compute variance explained (R²-like metric) for a
## given omics view (mRNA, miRNA, DNAme, SNV, etc.)

compute_var_explained <- function(view_name, TL_data, ZMu_0, ZMu, YTrg_list) {
  # 0) Check the view likelihood
  likelihood = likelihoodsLrn[view_name]

  # 1) Extract learned weights (W) and intercept weights (W0) for the view
  W <- TL_data$Fctrzn_Lrn_W[[view_name]] # Learned weight matrix
  W0 <- TL_data$Fctrzn_Lrn_W0[[view_name]] # Intercept (bias term)

  # 2) Estimate the target data reconstruction:
  #    Combine intercept + latent factors (ZMu_0, ZMu) with learned weights
  #    Note: cbind(ZMu_0, ZMu) adds intercept column to latent factors
  if (likelihood == "gaussian") {
    T_estimated <- cbind(ZMu_0, ZMu) %*% t(cbind(W0, W))
  }

  if (likelihood == "bernoulli") {
    logit_scores <- cbind(ZMu_0, ZMu) %*% t(cbind(W0, W))
    T_estimated <- plogis(logit_scores)
  }
  # 3) Get the normalized target data for this view
  T_normalized <- t(YTrg_list[[view_name]])
  # 4) Compute total sum of squares (SS_tot)
  #    This measures total variance of the true data (around its mean)
  SS_tot <- sum(
    (T_normalized - mean(T_normalized, na.rm = TRUE))^2,
    na.rm = TRUE
  )

  # 5) Compute residual sum of squares (RSS)
  #    This measures reconstruction error (difference between
  #    true and estimated data)
  RSS <- sum((T_normalized - T_estimated)^2, na.rm = TRUE)

  # 6) Variance explained = 1 - (RSS / SS_tot)
  #    This is equivalent to a global R² for the reconstruction
  var_explained <- 1 - RSS / SS_tot

  return(var_explained)
}

## Apply function to all available views

views <- c("mRNA", "miRNA", "DNAme", "SNV")
VarExpl_results <- sapply(
  views,
  compute_var_explained,
  TL_data = TL_data,
  ZMu_0 = ZMu_0,
  ZMu = ZMu,
  YTrg_list = YTrg_list
)

VarExpl_results
