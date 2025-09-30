mRNA_addVersion <- function(expdat, Lrndat) {
  #'
  #' get ensembl mRNA versions from learning dataset
  #' and attach to the correspnding mRNA ensembl id in the target dataset
  #'
  #' @param expdat the mRNA matrix from the target dataset with genes in rows
  #' @param Lrndat the mRNA W matrix from the learning dataset factorization with genes in rows
  #' @returns the target mRNA matrix with versions attached

  tmp = as.data.frame(do.call(rbind, strsplit(rownames(Lrndat), "[.]")))
  # match to stripped ids from target set
  tmp = data.frame(V1 = rownames(expdat)) %>%
    dplyr::left_join(tmp, by = c('V1')) %>%
    as.data.frame()
  # rename target dataset features
  rownames(expdat) = paste0(tmp$V1, '.', tmp$V2)

  # return tidied up matrix
  return(expdat)
}


## -------------------- SETUP --------------------

library(matrixStats)

# Safe log2 transform (in case of counts data)
safelog2 <- function(x) log2(x + 1)

# Helper: z-score a matrix with provided mean and sd
zscore_with_params <- function(X, mu, sd) {
  sd[sd == 0 | is.na(sd)] <- 1
  sweep(sweep(X, 1, mu, "-"), 1, sd, "/")
}

# Helper: z-score computed independently on X
zscore_by_rows <- function(X) {
  mu <- rowMeans(X, na.rm = TRUE)
  sd <- rowSds(X, na.rm = TRUE)
  zscore_with_params(X, mu, sd)
}

# Simple Wasserstein distance in 1D (quantile based)
wasserstein1d_vec <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  q <- seq(0, 1, length.out = max(nx, ny))
  qx <- quantile(x, q, na.rm = TRUE, names = FALSE)
  qy <- quantile(y, q, na.rm = TRUE, names = FALSE)
  mean(abs(qx - qy))
}

# KS statistic wrapper
ks_stat_vec <- function(x, y) {
  suppressWarnings(stats::ks.test(x, y)$statistic[[1]])
}

## -------------------- LOAD LEARNING VIEW --------------------
load_lrn_view <- function(LrnFctrnDir, meta, view) {
  # Load a single omics view (e.g. mRNA.csv) from the learning dataset
  f <- file.path(LrnFctrnDir, paste0(view, ".csv"))
  X <- as.matrix(read.csv(f, header = FALSE))
  rownames(X) <- meta[[paste0("ftrs_", view)]]
  colnames(X) <- meta$smpls
  X
}

## -------------------- ALIGN FEATURES --------------------
align_features <- function(Y_lrn, Y_trg) {
  # Restrict both matrices to common features and match row order
  common <- intersect(rownames(Y_lrn), rownames(Y_trg))
  Y_lrn <- Y_lrn[common, , drop = FALSE]
  Y_trg <- Y_trg[common, , drop = FALSE]
  list(lrn = Y_lrn, trg = Y_trg, features = common)
}

## -------------------- UNIFORMIZATION MODES --------------------
uniformize_pair <- function(
  Y_lrn,
  Y_trg,
  mode = c("anchored", "common", "separate")
) {
  # Apply normalization according to selected mode:
  # - "anchored": compute z-score params from learning set, apply to both
  # - "common": compute z-score params from combined (Lrn ∪ Trg), apply to both
  # - "separate": z-score each dataset independently
  mode <- match.arg(mode)
  if (mode == "anchored") {
    mu <- rowMeans(Y_lrn, na.rm = TRUE)
    sd <- rowSds(Y_lrn, na.rm = TRUE)
    Y_lrn_n <- zscore_with_params(Y_lrn, mu, sd)
    Y_trg_n <- zscore_with_params(Y_trg, mu, sd)
  } else if (mode == "common") {
    X_all <- cbind(Y_lrn, Y_trg)
    mu <- rowMeans(X_all, na.rm = TRUE)
    sd <- rowSds(X_all, na.rm = TRUE)
    Y_lrn_n <- zscore_with_params(Y_lrn, mu, sd)
    Y_trg_n <- zscore_with_params(Y_trg, mu, sd)
  } else {
    # "separate"
    Y_lrn_n <- zscore_by_rows(Y_lrn)
    Y_trg_n <- zscore_by_rows(Y_trg)
  }
  list(lrn = Y_lrn_n, trg = Y_trg_n)
}

## -------------------- SIMILARITY METRICS --------------------
similarity_metrics <- function(Y_lrn_n, Y_trg_n, with_distributions = FALSE) {
  # Compute feature-level summary statistics
  mean_lrn <- rowMeans(Y_lrn_n, na.rm = TRUE)
  mean_trg <- rowMeans(Y_trg_n, na.rm = TRUE)
  var_lrn <- rowVars(Y_lrn_n, na.rm = TRUE)
  var_trg <- rowVars(Y_trg_n, na.rm = TRUE)

  # Correlation of feature means and variances
  cor_mean <- suppressWarnings(cor(
    mean_lrn,
    mean_trg,
    method = "pearson",
    use = "complete.obs"
  ))
  cor_var <- suppressWarnings(cor(
    var_lrn,
    var_trg,
    method = "pearson",
    use = "complete.obs"
  ))

  out <- list(cor_mean = cor_mean, cor_var = cor_var)

  if (with_distributions) {
    # Optional: average Wasserstein and KS distances across features
    wass <- vapply(
      seq_len(nrow(Y_lrn_n)),
      function(i) {
        wasserstein1d_vec(Y_lrn_n[i, ], Y_trg_n[i, ])
      },
      numeric(1)
    )
    ks <- vapply(
      seq_len(nrow(Y_lrn_n)),
      function(i) {
        ks_stat_vec(as.numeric(Y_lrn_n[i, ]), as.numeric(Y_trg_n[i, ]))
      },
      numeric(1)
    )

    out$wasserstein_mean <- mean(wass, na.rm = TRUE)
    out$ks_mean <- mean(ks, na.rm = TRUE)
  }
  out
}

## -------------------- END-TO-END PIPELINE --------------------
run_similarity <- function(
  LrnFctrnDir,
  YLrn_list,
  raw_YTrg_list,
  views = c("mRNA", "miRNA", "DNAme", "SNV"),
  uniformization = "anchored",
  with_distributions = FALSE
) {
  # Full similarity workflow for multiple omics views

  # Load metadata of learning dataset
  meta_lrn <- readRDS(file.path(LrnFctrnDir, "expdat_meta.rds"))

  results <- lapply(views, function(view) {
    # 1) Load learning dataset
    Y_lrn <- YLrn_list[[view]]

    # 2) Get target dataset (already loaded into raw_YTrg_list)
    Y_trg <- raw_YTrg_list[[view]]

    # (Optional) apply log2 if these are raw counts
    # Y_lrn <- safelog2(Y_lrn); Y_trg <- safelog2(Y_trg)

    # 3) Align common features
    aligned <- align_features(Y_lrn, Y_trg)
    Y_lrn_a <- aligned$lrn
    Y_trg_a <- aligned$trg

    # 4) Normalize / uniformize
    uni <- uniformize_pair(Y_lrn_a, Y_trg_a, mode = uniformization)

    # 5) Compute metrics
    mets <- similarity_metrics(
      uni$lrn,
      uni$trg,
      with_distributions = with_distributions
    )
    c(mets, list(n_features = nrow(Y_lrn_a), view = view))
  })

  names(results) <- views
  results
}

## COMPUTE METRICS

metrics <- run_similarity(
  LrnFctrnDir,
  sample_YLrn_list,
  raw_YTrg_list,
  views = views,
  uniformization = "anchored",
  with_distributions = TRUE
)

print(metrics)
