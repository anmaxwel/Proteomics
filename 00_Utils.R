## R/utils.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
})

read_inputs <- function(protein_groups_path, sample_sheet_path, contrasts_path) {
  if (!file.exists(protein_groups_path)) stop("Missing: ", protein_groups_path)
  if (!file.exists(sample_sheet_path))   stop("Missing: ", sample_sheet_path)
  if (!file.exists(contrasts_path))      stop("Missing: ", contrasts_path)

  pg <- data.table::fread(protein_groups_path, data.table = FALSE)

  samples <- readr::read_csv(sample_sheet_path, show_col_types = FALSE) %>%
    mutate(
      condition = factor(condition),
      sample    = as.character(sample),
      mq_column = as.character(mq_column)
    )

  contrasts <- readr::read_csv(contrasts_path, show_col_types = FALSE) %>%
    mutate(across(everything(), as.character))

  list(pg = pg, samples = samples, contrasts = contrasts)
}

first_uniprot <- function(protein_ids) {
  # MaxQuant Protein.IDs often looks like "P12345;Q9XXXX;..."
  x <- str_split_fixed(protein_ids, ";", 2)[, 1]
  x <- str_replace(x, "-\\d+$", "")  # drop isoform suffix like P12345-2
  x
}

safe_scale <- function(x) {
  if (sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}
