setwd("C:/Users/wenzh/Desktop/test")
seu_obj <- readRDS("./CD8_2021Nature+2025CC+2023eLife+thisstudy+2023JITC_Harmony_annotation.rds")
unique(seu_obj$patient_id)

# Calculate the number of cells expressing each gene (count > 0)
gene_expression <- rowSums(GetAssayData(seu_obj, assay = "RNA", slot = "counts") > 0)
# Set minimum cells threshold (e.g., 20; increase for stricter filtering)
min_cells <- 160
# Select genes to keep
genes_to_keep <- rownames(seu_obj)[gene_expression >= min_cells]
# Subset the Seurat object to these genes
seu_obj <- seu_obj[genes_to_keep, ]

counts <- GetAssayData(seu_obj, slot = "counts", assay = "RNA")
meta <- seu_obj@meta.data
meta$group <- paste(meta$patient_id, meta$cell_names, sep = "_")
group_levels <- unique(meta$group)

avg_expr_list <- lapply(group_levels, function(g) {
  cells_in_group <- rownames(meta)[meta$group == g]
  if (length(cells_in_group) > 100) {  # Only retain samples with a sufficient number of cells
    rowMeans(counts[, cells_in_group, drop = FALSE])
  } else {
    NULL
  }
})

names(avg_expr_list) <- group_levels
avg_expr_list <- avg_expr_list[!sapply(avg_expr_list, is.null)]
avg_expr_matrix <- do.call(cbind, avg_expr_list)

gene_mean <- rowMeans(avg_expr_matrix)
gene_sd <- apply(avg_expr_matrix, 1, sd)
gene_cv <- gene_sd / gene_mean

gene_stats <- data.frame(
  gene = rownames(avg_expr_matrix),
  mean_expr = gene_mean,
  cv = gene_cv,
  stringsAsFactors = FALSE
)

library(ggplot2)
library(ggrepel)
library(dplyr)

ggplot(gene_stats, aes(x = mean_expr, y = cv)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  labs(title = "Reference Gene Selection", x = "Mean Expression (log10)", y = "Coefficient of Variation (CV)")

ref_genes <- gene_stats %>%
  filter(mean_expr > 20, cv < 0.6) %>%
  arrange(cv)

genes_to_label <- ref_genes$gene
gene_stats$label <- ifelse(gene_stats$gene %in% genes_to_label, gene_stats$gene, NA)

point_color <- "#4A4A4A"        
highlight_color <- "#D64045"  
label_color <- "#D64045"

ggplot(gene_stats, aes(x = mean_expr, y = cv)) +
  geom_point(color = point_color, alpha = 0.5, size = 1.8) +
  geom_point(
    data = subset(gene_stats, gene %in% genes_to_label),
    color = highlight_color,
    size = 3
  ) +
  geom_text_repel(
    data = subset(gene_stats, gene %in% genes_to_label),
    aes(label = label),
    size = 4,
    fontface = "bold",
    color = label_color,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.25,
    segment.color = "grey40"
  ) +
  scale_x_log10() +
  labs(
    title = "Reference Gene Selection",
    x = "Mean Expression (log10)",
    y = "Coefficient of Variation (CV)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),

    axis.line = element_line(size = 0.6),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank()
  )










