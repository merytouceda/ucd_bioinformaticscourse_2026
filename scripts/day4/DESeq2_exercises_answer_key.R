# ============================================================
# DESeq2 Follow-up Exercises — Answer Key
# UCD Bioinformatics Summer Course — Week 2, Day 2
# ============================================================
# NOTE: This script assumes dds, vsd, res_D7_vs_D0, and
# res_D14_vs_D0 are already in your environment from the
# main tutorial. Run that script first before running this one.
# ============================================================

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)  # install.packages("ggrepel") if needed


# ============================================================
# EASY EXERCISES
# ============================================================

# ------------------------------------------------------------
# Exercise 1: Change the fold-change threshold on the volcano plot
# Try: fc_threshold <- 0.5  or  fc_threshold <- 2
# ------------------------------------------------------------

fc_threshold <- 2   # <-- change this and re-run

volcano_data_ex1 <- bind_rows(
  as.data.frame(res_D7_vs_D0)  %>% rownames_to_column("gene") %>% mutate(contrast = "D7 vs D0"),
  as.data.frame(res_D14_vs_D0) %>% rownames_to_column("gene") %>% mutate(contrast = "D14 vs D0")
) %>%
  filter(!is.na(padj)) %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange >  fc_threshold ~ "Up",
      padj < 0.05 & log2FoldChange < -fc_threshold ~ "Down",
      TRUE ~ "NS"
    )
  )

# Count DEGs at this threshold
volcano_data_ex1 %>%
  filter(sig != "NS") %>%
  count(contrast, sig)

ggplot(volcano_data_ex1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  facet_wrap(~ contrast) +
  labs(
    title = paste0("Volcano plot — |log2FC| > ", fc_threshold),
    x = "log2 Fold Change", y = "-log10 adjusted p-value", color = "Direction"
  ) +
  theme_bw()

# ------------------------------------------------------------
# Exercise 2: Change heatmap size (n = 25 or n = 100)
# ------------------------------------------------------------

n_genes <- 25   # <-- change this and re-run

top_genes_ex2 <- bind_rows(
  as.data.frame(res_D7_vs_D0)  %>% rownames_to_column("gene"),
  as.data.frame(res_D14_vs_D0) %>% rownames_to_column("gene")
) %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  distinct(gene, .keep_all = TRUE) %>%
  slice_head(n = n_genes) %>%
  pull(gene)

mat_ex2        <- assay(vsd)[top_genes_ex2, ]
mat_scaled_ex2 <- t(scale(t(mat_ex2)))

annotation_col <- data.frame(Day = pheno_data$Day)
rownames(annotation_col) <- rownames(pheno_data)

pheatmap(
  mat_scaled_ex2,
  annotation_col = annotation_col,
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  fontsize_row   = 7,
  main           = paste0("Top ", n_genes, " DEGs (scaled VST counts)")
)

# ------------------------------------------------------------
# Exercise 3: Change the padj threshold
# Try: padj_threshold <- 0.01
# ------------------------------------------------------------

padj_threshold <- 0.01   # <-- change this and re-run

# Count how many DEGs survive at the stricter threshold
bind_rows(
  as.data.frame(res_D7_vs_D0)  %>% rownames_to_column("gene") %>% mutate(contrast = "D7 vs D0"),
  as.data.frame(res_D14_vs_D0) %>% rownames_to_column("gene") %>% mutate(contrast = "D14 vs D0")
) %>%
  filter(!is.na(padj)) %>%
  mutate(
    sig = case_when(
      padj < padj_threshold & log2FoldChange >  1 ~ "Up",
      padj < padj_threshold & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  count(contrast, sig)

# ------------------------------------------------------------
# Exercise 4: Label top 10 genes on the volcano plot with ggrepel
# ------------------------------------------------------------

# Identify top 10 most significant genes per contrast
top_labels <- bind_rows(
  as.data.frame(res_D7_vs_D0)  %>% rownames_to_column("gene") %>% mutate(contrast = "D7 vs D0"),
  as.data.frame(res_D14_vs_D0) %>% rownames_to_column("gene") %>% mutate(contrast = "D14 vs D0")
) %>%
  filter(!is.na(padj)) %>%
  group_by(contrast) %>%
  slice_min(padj, n = 10) %>%
  ungroup()

volcano_data_labels <- bind_rows(
  as.data.frame(res_D7_vs_D0)  %>% rownames_to_column("gene") %>% mutate(contrast = "D7 vs D0"),
  as.data.frame(res_D14_vs_D0) %>% rownames_to_column("gene") %>% mutate(contrast = "D14 vs D0")
) %>%
  filter(!is.na(padj)) %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

ggplot(volcano_data_labels, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data      = top_labels,
    aes(label = gene),
    size      = 3,
    color     = "black",
    max.overlaps = 15
  ) +
  facet_wrap(~ contrast) +
  labs(x = "log2 Fold Change", y = "-log10 adjusted p-value", color = "Direction") +
  theme_bw()


# ============================================================
# MEDIUM EXERCISES
# ============================================================

# ------------------------------------------------------------
# Exercise 5: Count DEGs per contrast and direction as a table
# ------------------------------------------------------------

deg_counts <- bind_rows(
  as.data.frame(res_D7_vs_D0)  %>% rownames_to_column("gene") %>% mutate(contrast = "D7 vs D0"),
  as.data.frame(res_D14_vs_D0) %>% rownames_to_column("gene") %>% mutate(contrast = "D14 vs D0")
) %>%
  filter(!is.na(padj)) %>%
  mutate(
    direction = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  count(contrast, direction) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(Total_DEG = Up + Down)

print(deg_counts)

# ------------------------------------------------------------
# Exercise 6: Save significant DEGs to CSV files
# ------------------------------------------------------------

save_degs <- function(res, contrast_label, padj_cut = 0.05, lfc_cut = 1) {
  as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj), padj < padj_cut, abs(log2FoldChange) > lfc_cut) %>%
    arrange(padj) %>%
    mutate(direction = if_else(log2FoldChange > 0, "Up", "Down"))
}

degs_D7_vs_D0  <- save_degs(res_D7_vs_D0,  "D7 vs D0")
degs_D14_vs_D0 <- save_degs(res_D14_vs_D0, "D14 vs D0")

write_csv(degs_D7_vs_D0,  "DEGs_D7_vs_D0.csv")
write_csv(degs_D14_vs_D0, "DEGs_D14_vs_D0.csv")

cat("D7 vs D0: ", nrow(degs_D7_vs_D0),  "DEGs saved\n")
cat("D14 vs D0:", nrow(degs_D14_vs_D0), "DEGs saved\n")

# ------------------------------------------------------------
# Exercise 7: Add a third contrast — D14 vs D7
# ------------------------------------------------------------

res_D14_vs_D7 <- results(dds, contrast = c("Day", "D14", "D7"))
summary(res_D14_vs_D7)

# Volcano for D14 vs D7
as.data.frame(res_D14_vs_D7) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj)) %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "D14 vs D7",
    x = "log2 Fold Change", y = "-log10 adjusted p-value", color = "Direction"
  ) +
  theme_bw()

# ------------------------------------------------------------
# Exercise 8: Colour the PCA plot by replicate to check for batch effects
# ------------------------------------------------------------

# Add a replicate column to pheno_data if it doesn't exist
# (assumes replicates are in the same order within each Day group)
pheno_data_rep <- pheno_data %>%
  group_by(Day) %>%
  mutate(Replicate = factor(row_number())) %>%
  ungroup()

# Re-run VST with updated colData
colData(vsd)$Replicate <- pheno_data_rep$Replicate

pcaData_rep <- plotPCA(vsd, intgroup = c("Day", "Replicate"), returnData = TRUE)
percentVar  <- round(100 * attr(pcaData_rep, "percentVar"))

ggplot(pcaData_rep, aes(x = PC1, y = PC2, color = Replicate, shape = Day)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  labs(title = "PCA coloured by replicate — checking for batch effects") +
  theme_bw()
# If replicates cluster by Day (not by replicate number), there's no batch effect


# ============================================================
# HARD EXERCISES
# ============================================================

# ------------------------------------------------------------
# Exercise 9: Separate heatmaps for D7-specific vs D14-specific DEGs
# ------------------------------------------------------------

# D7-specific: significant in D7 vs D0 but NOT in D14 vs D0
sig_D7 <- as.data.frame(res_D7_vs_D0) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(gene)

sig_D14 <- as.data.frame(res_D14_vs_D0) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(gene)

only_D7  <- setdiff(sig_D7,  sig_D14)  # unique to D7
only_D14 <- setdiff(sig_D14, sig_D7)   # unique to D14
shared   <- intersect(sig_D7, sig_D14) # significant at both time points

cat("D7-specific DEGs: ", length(only_D7),  "\n")
cat("D14-specific DEGs:", length(only_D14), "\n")
cat("Shared DEGs:      ", length(shared),   "\n")

# Helper to plot heatmap for a gene set (capped at max_n genes)
plot_heatmap_subset <- function(genes, title, max_n = 50) {
  genes <- head(genes, max_n)
  if (length(genes) < 2) {
    message("Not enough genes to plot: ", title)
    return(invisible(NULL))
  }
  mat        <- assay(vsd)[genes, ]
  mat_scaled <- t(scale(t(mat)))
  pheatmap(
    mat_scaled,
    annotation_col = annotation_col,
    cluster_rows   = TRUE,
    cluster_cols   = TRUE,
    show_rownames  = TRUE,
    show_colnames  = TRUE,
    fontsize_row   = 7,
    main           = title
  )
}

plot_heatmap_subset(only_D7,  "D7-specific DEGs (top 50)")
plot_heatmap_subset(only_D14, "D14-specific DEGs (top 50)")
plot_heatmap_subset(shared,   "Shared DEGs (top 50)")

# ------------------------------------------------------------
# Exercise 10: Apply lfcShrink and compare MA plots
# ------------------------------------------------------------

# Shrinkage stabilises log2FC estimates for low-count genes
# using the apeglm method (most recommended)
# BiocManager::install("apeglm") if needed

res_D7_shrunk <- lfcShrink(dds, coef = "Day_D7_vs_D0", type = "apeglm")

par(mfrow = c(1, 2))
plotMA(res_D7_vs_D0,  main = "D7 vs D0 — unshrunken", ylim = c(-5, 5))
plotMA(res_D7_shrunk, main = "D7 vs D0 — shrunken",   ylim = c(-5, 5))
par(mfrow = c(1, 1))

# Notice: after shrinkage, noisy low-count genes (left side) have their
# extreme fold changes pulled toward 0. High-count genes are barely affected.
# Shrinkage gives more reliable ranked gene lists for downstream analysis.
