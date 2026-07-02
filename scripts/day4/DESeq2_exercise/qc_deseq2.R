# UCD - bioinformatics summer course
# Week 2 - Day 2: RNA-seq
# Differential gene expression with DESeq2

#BiocManager::install("DESeq2", force = TRUE)
library(DESeq2)
library(tidyverse)
library(pheatmap)

setwd("~/Documents/GitHub/ucd_bioinformaticscourse_2026/scripts/day4/DESeq2_exercise/")

# ==============================
# READ IN DATA
# ==============================
# Read in gene matrix
Gene_matrix <- read.table("../../data/Gene_matrix.txt", sep="\t", header=TRUE)
dim(Gene_matrix)
#[1] 28524     9

# Read in pheno data
pheno_data <- read.table("../../data/metadata_day4.txt", sep="\t", header=TRUE) %>%
  mutate(Day = factor(Day)) %>%
  column_to_rownames("Sample_ID") %>%
  mutate(Sample_ID = rownames(.))
dim(pheno_data)
#[1] 9 3

levels(pheno_data$Day)
#[1] "D0"  "D14" "D7" 

# Check that the sample order matches
all(rownames(pheno_data) == colnames(Gene_matrix))
#[1] TRUE


# ==============================
# DESeq2
# ==============================
# Prep data for DESeq2
# DESeqDataSetFromMatrix requires a plain matrix, so convert
dds <- DESeqDataSetFromMatrix(Gene_matrix, pheno_data, design = ~ Day)

# Filter genes
# Keep genes that have at least 10 counts across the smallest group size
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dim(dds)
#[1] 15922     9

# Run DESeq2
dds <- DESeq(dds)


# ==============================
# PCA
# ==============================
# Variance-stabilizing transformation
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup="Day", returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

# create plot and save to object
pcaPlot <- ggplot(pcaData, aes(x=PC1, y=PC2, color=Day)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

# view
pcaPlot

# save
ggsave("PCA_plot.pdf", pcaPlot, device="pdf")


# ==============================
# DEG results
# ==============================

# Extract results for each contrast
res_D7_vs_D0  <- results(dds, contrast = c("Day", "D7",  "D0"))
res_D14_vs_D0 <- results(dds, contrast = c("Day", "D14", "D0"))

# Summary of DEGs
summary(res_D7_vs_D0)
summary(res_D14_vs_D0)

#------------------------------ MA plots

par(mfrow = c(1, 2))

plotMA(res_D7_vs_D0,  main = "D7 vs D0",  ylim = c(-5, 5))
plotMA(res_D14_vs_D0, main = "D14 vs D0", ylim = c(-5, 5))

par(mfrow = c(1, 1))

#------------------------------ Volcano plots

# Helper function to prepare results for volcano plot
prep_volcano <- function(res, contrast_label) {
  as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(
      contrast = contrast_label,
      sig = case_when(
        padj < 0.05 & log2FoldChange >  1 ~ "Up",
        padj < 0.05 & log2FoldChange < -1 ~ "Down",
        TRUE ~ "NS"
      )
    )
}

volcano_data <- bind_rows(
  prep_volcano(res_D7_vs_D0,  "D7 vs D0"),
  prep_volcano(res_D14_vs_D0, "D14 vs D0")
)

volcanoPlot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  facet_wrap(~ contrast) +
  labs(x = "log2 Fold Change", y = "-log10 adjusted p-value", color = "Direction") +
  theme_bw()

# view
volcanoPlot

# save
ggsave("Volcano_plot.pdf", volcanoPlot, device="pdf")

#------------------------------ Heatmap

# Get top 50 DEGs by padj pooled across both contrasts
top_genes <- bind_rows(
  as.data.frame(res_D7_vs_D0)  %>% rownames_to_column("gene"),
  as.data.frame(res_D14_vs_D0) %>% rownames_to_column("gene")
) %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  distinct(gene, .keep_all = TRUE) %>%
  slice_head(n = 50) %>%
  pull(gene)

# Extract VST counts for top genes
mat <- assay(vsd)[top_genes, ]

# Scale by row (z-score)
mat_scaled <- t(scale(t(mat)))

# Annotation for columns
annotation_col <- data.frame(Day = pheno_data$Day)
rownames(annotation_col) <- rownames(pheno_data)

pheatmap(
  mat_scaled,
  annotation_col  = annotation_col,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  fontsize_row    = 7,
  main            = "Top 50 DEGs (scaled VST counts)"
)

