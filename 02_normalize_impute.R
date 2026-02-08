## R/02_normalize_impute.R

suppressPackageStartupMessages({
  library(vsn)
  library(MSnbase)
  library(dplyr)
})

filter_by_group_presence <- function(mat, samples, min_non_na_per_group = 2) {
  stopifnot(ncol(mat) == nrow(samples))

  keep <- rep(TRUE, nrow(mat))
  for (grp in levels(samples$condition)) {
    cols <- which(samples$condition == grp)
    non_na <- rowSums(!is.na(mat[, cols, drop = FALSE]))
    keep <- keep & (non_na >= min_non_na_per_group)
  }
  mat[keep, , drop = FALSE]
}

normalize_vsn <- function(mat) {
  # VSN expects non-negative values, LFQ/intensities are non-negative
  vsn::justvsn(mat)
}

impute_minprob <- function(mat_norm, q = 0.01) {
  MSnbase::impute(mat_norm, method = "MinProb", q = q)
}

preprocess_proteins <- function(pg, samples) {
  # override name from 01 file, but keep behavior as a single entry point
  pg_clean <- clean_maxquant_proteinGroups(pg)
  built <- build_lfq_matrix(pg_clean, samples)

  raw_mat <- built$mat
  raw_mat[raw_mat <= 0] <- NA

  mat_filt <- filter_by_group_presence(raw_mat, samples, min_non_na_per_group = 2)

  mat_norm <- normalize_vsn(mat_filt)
  mat_imp  <- impute_minprob(mat_norm, q = 0.01)

  log2_mat <- log2(mat_imp)

  list(
    pg_clean = pg_clean,
    raw_mat  = raw_mat,
    mat_filt = mat_filt,
    mat_norm = mat_norm,
    mat_imp  = mat_imp,
    log2_mat = log2_mat,
    annot    = built$annot %>% filter(uniprot %in% rownames(log2_mat))
  )
}
