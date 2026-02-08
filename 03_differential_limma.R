## R/03_differential_limma.R

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(readr)
})

run_limma <- function(log2_mat, samples, contrasts) {
  stopifnot(all(colnames(log2_mat) == samples$sample))

  design <- model.matrix(~ 0 + condition, data = samples)
  colnames(design) <- levels(samples$condition)

  fit <- lmFit(log2_mat, design)

  # build contrast expressions from contrasts.csv
  expr <- paste0(contrasts$numerator, " - ", contrasts$denominator)
  names(expr) <- contrasts$contrast

  cm <- makeContrasts(contrasts = expr, levels = design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)

  list(fit = fit2, design = design, contrasts = contrasts)
}

write_results <- function(fit_obj, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (contrast_name in fit_obj$contrasts$contrast) {
    tt <- topTable(fit_obj$fit, coef = contrast_name, number = Inf, sort.by = "P")
    out <- file.path(out_dir, paste0(contrast_name, ".csv"))
    readr::write_csv(tibble::rownames_to_column(tt, "feature_id"), out)
  }
}
