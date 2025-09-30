# --- REQUIREMENTS ---
library(MOFA2)
library(here)
library(ggplot2)

## -------------------- LOAD TARGET DATA --------------------
TrgtDir <- here("data/TrgtDir")
TrgtFctrnDir <- file.path(TrgtDir, "Trg_PAAD_SKCM_SS9_5000D")
TrgtFctrnRdsObjectsDir <- here(TrgtFctrnDir, "RDS_objects")

raw_YTrg_list <- readRDS(file.path(TrgtFctrnRdsObjectsDir, "raw_YTrg_list.RDS"))

# Use raw_YTrg_list as-is if already normalized/transformed
YTrg_list <- raw_YTrg_list

# Basic checks
stopifnot(is.list(YTrg_list), length(YTrg_list) >= 1)
views <- names(YTrg_list)
message("Detected views: ", paste(views, collapse = ", "))

# IMPORTANT: MOFA expects features as rows and samples as columns.
# If yours are transposed, uncomment the next line:
# YTrg_list <- lapply(YTrg_list, t)

## -------------------- CREATE MOFA OBJECT --------------------
MOFAobject <- create_mofa(YTrg_list)

## -------------------- OPTIONS --------------------
# Data options: scale/center can help align ranges across views
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE
data_opts$center_groups <- TRUE
head(data_opts)

# Model options: set number of factors (cannot exceed number of samples)
n_samples <- ncol(YTrg_list[[1]])
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- min(10, n_samples)
head(model_opts)

# Auto-detect likelihoods per view: Bernoulli for 0/1 matrices, Gaussian otherwise
auto_lik <- lapply(names(YTrg_list), function(vn) {
  X <- YTrg_list[[vn]]
  u <- unique(na.omit(as.vector(X)))
  if (length(u) > 0 && all(u %in% c(0, 1))) "bernoulli" else "gaussian"
})

names(auto_lik) <- names(YTrg_list)
model_opts$likelihoods <- unlist(auto_lik)


# Training options
train_opts <- get_default_training_options(MOFAobject)
#train_opts$convergence_mode <- "medium" # "fast" | "medium" | "slow"
#train_opts$maxiter <- 10000
#train_opts$seed <- 42

## -------------------- PREPARE & RUN --------------------
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path(tempdir(), "model.hdf5")

model <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

## -------------------- SAVE OUTPUTS --------------------
out_dir <- file.path(TrgtFctrnDir, "MOFA_direct")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(model, file.path(out_dir, "MOFA_direct_model.rds"))

factors_df <- get_factors(model, factors = "all", as.data.frame = TRUE)
weights_list <- get_weights(model, as.data.frame = FALSE) # list by view

saveRDS(factors_df, file.path(out_dir, "MOFA_direct_factors_df.rds"))
saveRDS(weights_list, file.path(out_dir, "MOFA_direct_weights_list.rds"))

## -------------------- VISUALS (no sample labels) --------------------
# If you have sample metadata (optional), set it here, e.g.:
# sample_metadata <- data.frame(Sample=colnames(YTrg_list[[1]]),
#                               Group = ..., row.names=colnames(YTrg_list[[1]]))
# model <- set_sample_metadata(model, sample_metadata)

# Factor scatter: MOFA2 does not print sample labels by default
p12 <- plot_factors(
  model,
  factors = c(1, 2),
  color_by = NULL, # set to "Group" if you added sample_metadata$Group
  dot_size = 2,
  legend = TRUE
) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(out_dir, "MOFA_direct_F1_F2.png"),
  p12,
  width = 6,
  height = 5,
  dpi = 300
)

# Variance explained per factor and per view
p_var <- plot_variance_explained(model)
ggsave(
  file.path(out_dir, "MOFA_direct_variance_explained.png"),
  p_var,
  width = 7,
  height = 5,
  dpi = 300
)

# Example: top weights for a specific view (change "mRNA" to a real view name in your data)
if ("mRNA" %in% names(weights_list)) {
  p_w <- plot_weights(model, view = "mRNA", factor = 1, nfeatures = 20)
  ggsave(
    file.path(out_dir, "MOFA_direct_topweights_mRNA_F1.png"),
    p_w,
    width = 7,
    height = 5,
    dpi = 300
  )
}
