## R/01_import_qc.R

suppressPackageStartupMessages({
  library(dplyr)
})

clean_maxquant_proteinGroups <- function(pg) {
  # Common MaxQuant filters
  pg %>%
    filter(is.na(`Potential contaminant`) | `Potential contaminant` != "+") %>%
    filter(is.na(Reverse) | Reverse != "+") %>%
    filter(is.na(`Only identified by site`) | `Only identified by site` != "+")
}

build_lfq_matrix <- function(pg_clean, samples) {
  missing_cols <- setdiff(samples$mq_column, colnames(pg_clean))
  if (length(missing_cols) > 0) {
    stop("These mq_column values were not found in proteinGroups: ",
         paste(missing_cols, collapse = ", "))
  }

  mat <- pg_clean %>%
    dplyr::select(Protein.IDs, all_of(samples$mq_column))

  # ensure numeric
  mat_num <- mat %>%
    mutate(across(all_of(samples$mq_column), ~ as.numeric(.)))

  rownames_mat <- first_uniprot(mat_num$Protein.IDs)
  m <- as.matrix(mat_num[, samples$mq_column, drop = FALSE])
  rownames(m) <- rownames_mat
  colnames(m) <- samples$sample

  list(
    mat = m,
    annot = tibble(
      Protein.IDs = mat_num$Protein.IDs,
      uniprot     = rownames_mat
    )
  )
}

preprocess_proteins <- function(pg, samples) {
  pg_clean <- clean_maxquant_proteinGroups(pg)
  built <- build_lfq_matrix(pg_clean, samples)

  list(
    pg_clean = pg_clean,
    mat      = built$mat,
    annot    = built$annot
  )
}
