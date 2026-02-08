## R/04_plots.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tibble)
})

make_all_plots <- function(raw_pg, cleaned, fit, samples, out_dir) {
  dir.create(out_dir, showWarnings = FALSE)

  # PCA
  pca <- prcomp(t(cleaned$log2_mat), scale. = TRUE)
  pca_df <- data.frame(
    sample = rownames(pca$x),
    condition = samples$condition,
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )

  g <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "PCA (log2, VSN, MinProb imputed)")

  ggsave(file.path(out_dir, "PCA.png"), g, width = 7, height = 5, dpi = 300)

  # Sample correlation heatmap
  cor_mat <- cor(cleaned$log2_mat, use = "pairwise.complete.obs")
  png(file.path(out_dir, "sample_correlation.png"), width = 900, height = 800)
  pheatmap(cor_mat)
  dev.off()

  # Volcano for first contrast
  c1 <- fit$contrasts$contrast[1]
  tt <- topTable(fit$fit, coef = c1, number = Inf, sort.by = "P") %>%
    rownames_to_column("feature_id") %>%
    mutate(sig = adj.P.Val < 0.05)

  v <- ggplot(tt, aes(x = logFC, y = -log10(P.Value), color = sig)) +
    geom_point(size = 1.2, alpha = 0.8) +
    theme_minimal() +
    labs(title = paste0("Volcano: ", c1))

  ggsave(file.path(out_dir, paste0("volcano_", c1, ".png")), v, width = 7, height = 5, dpi = 300)

  # Heatmap of top 50 by adjusted p-value
  top_ids <- tt %>% arrange(adj.P.Val) %>% head(50) %>% pull(feature_id)
  hm <- cleaned$log2_mat[top_ids, , drop = FALSE]
  ann <- data.frame(condition = samples$condition)
  rownames(ann) <- samples$sample

  png(file.path(out_dir, paste0("heatmap_top50_", c1, ".png")), width = 1100, height = 900)
  pheatmap(hm, annotation_col = ann, show_rownames = TRUE)
  dev.off()
}
