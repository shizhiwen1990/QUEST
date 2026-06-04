library(tidyverse)
library(Matrix)
ref_gene_names <- c("MALAT1","B2M","RPL28","RPL10","RPS27","RPS12","TMSB4X","RPLP1")

counts <- GetAssayData(seu_obj, slot = "counts", assay = "RNA")
counts_cpm <- counts %*% Diagonal(x = 1e6 / colSums(counts))
ref_gene_names <- ref_gene_names[ref_gene_names %in% rownames(counts)]
log2_cpm <- log2(counts_cpm + 1)
log2_cpm_sub <- log2_cpm[ref_gene_names, ]

expr_long <- as.matrix(log2_cpm_sub) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "expression")

p <- expr_long %>%
  ggplot(aes(x = expression)) +
  
  # Smooth density only
  geom_density(
    aes(y = after_stat(density * 100)),
    fill = "#A6CEE3",   # 柔和填充
    color = "#1F4E79",   # 深蓝线条
    alpha = 0.6,
    linewidth = 1
  ) +
  
  facet_wrap(~ gene, scales = "fixed", ncol = 4) +
  
  labs(
    x = "log2(CPM+1) expression in single cells",
    y = "Percentage",
    title = "Distribution of Stable Genes in Single Cells"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title    = element_text(size = 12, face = "bold"),
    axis.text     = element_text(size = 10),
    strip.text    = element_text(size = 11, face = "italic"),
    panel.grid    = element_blank(),
    axis.line     = element_line(size = 0.4),
    plot.margin   = margin(10, 10, 10, 10)
  )

p <- expr_long %>%
  ggplot(aes(x = expression)) +
  
  # Smooth density only
  geom_density(
    aes(y = after_stat(density * 100)),
    fill = "#A6CEE3",   # 柔和填充
    color = "#1F4E79",   # 深蓝线条
    alpha = 0.6,
    linewidth = 1
  ) +
  
  facet_wrap(~ gene, scales = "fixed", ncol = 4) +
  
  labs(
    x = "log2(CPM+1) expression in single cells",
    y = "Percentage",
    title = "Distribution of Stable Genes in Single Cells"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title    = element_text(size = 12, face = "bold"),
    axis.text     = element_text(size = 10),
    strip.text    = element_text(size = 11, face = "italic"),
    panel.grid    = element_blank(),
    axis.line     = element_line(size = 0.4),
    plot.margin   = margin(10, 10, 10, 10)
  )

p

