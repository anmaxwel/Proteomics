## R/05_enrichment.R

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
})

uniprot_to_symbol <- function(uniprot_ids, species = "human") {
  if (species != "human") stop("Only human mapping implemented in this template.")
  map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(uniprot_ids),
    keytype = "UNIPROT",
    columns = c("SYMBOL")
  )
  map <- map %>% filter(!is.na(SYMBOL)) %>% distinct(UNIPROT, SYMBOL)
  map
}

run_enrichment_go <- function(fit, annot, species = "human", out_dir = "results/enrichment") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  map <- uniprot_to_symbol(annot$uniprot, species = species)

  for (contrast_name in fit$contrasts$contrast) {
    tt <- topTable(fit$fit, coef = contrast_name, number = Inf, sort.by = "P") %>%
      tibble::rownames_to_column("uniprot")

    sig <- tt %>%
      filter(adj.P.Val < 0.05) %>%
      inner_join(map, by = "uniprot") %>%
      distinct(SYMBOL) %>%
      pull(SYMBOL)

    if (length(sig) < 10) next

    eg <- bitr(sig, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    if (nrow(eg) < 10) next

    ego <- enrichGO(
      gene = eg$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      readable = TRUE
    )

    if (is.null(ego) || nrow(ego@result) == 0) next

    png(file.path(out_dir, paste0("GO_BP_", contrast_name, ".png")), width = 1200, height = 900)
    print(dotplot(ego, showCategory = 15) + ggtitle(paste0("GO BP: ", contrast_name)))
    dev.off()

    write.csv(ego@result, file.path(out_dir, paste0("GO_BP_", contrast_name, ".csv")), row.names = FALSE)
  }
}
