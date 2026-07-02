# UCD - bioinformatics summer course
# Week 2 - Day 2: RNA-seq
# Differential gene expression with DESeq2

#BiocManager::install("DESeq2", force = TRUE)
# when prompted with update pacakes (a/s/n:) type n on the console
# Try this installation if the previous one takes forever or doesn't work: 
# BiocManager::install("DESeq2", force = TRUE, ask = FALSE, update = FALSE, Ncpus = 4)
library(DESeq2)
library(tidyverse)
library(pheatmap)

setwd("~/Documents/GitHub/ucd_bioinformaticscourse_2026/scripts/day4/DESeq2_exercise/")

# ==============================
# READ IN DATA
# ==============================
# Read in gene matrix
Gene_matrix <- read_table("../../../data/Gene_matrix.txt", sep="\t", header=TRUE)
dim(Gene_matrix)
#[1] 28524     9

# Read in pheno data
pheno_data <- read_table("../../../data/metadata_day4.txt", sep="\t", header=TRUE) %>%
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
# DESeqDataSetFromMatrix requires three matched inputs:
#
# 1. Gene_matrix: a raw (non-normalized) integer count matrix
#    - rows = genes/features, columns = samples
#    - values must be raw counts (not TPM/FPKM/normalized values),
#      since DESeq2's internal size-factor normalization assumes
#      count data with a negative binomial distribution
#
# 2. pheno_data: a data frame of sample metadata (colData)
#    - rows = samples, in the SAME ORDER as the columns of Gene_matrix
#    - rownames(pheno_data) must exactly match colnames(Gene_matrix)
#    - must contain a column named "Day" (or whatever the design variable is),
#      ideally as a factor, since this defines the groups being compared
#
# 3. design = ~ Day: the model formula
#    - tells DESeq2 which column in pheno_data to use to model
#      variation in counts across samples (here, comparing across
#      levels/timepoints of "Day")
#    - if there are other variables to control for (e.g. batch),
#      they should be included as ~ batch + Day, with the variable
#      of interest listed last
#
# Output: dds is a DESeqDataSet object (extends SummarizedExperiment).
# At this point it just bundles the raw counts + metadata + design
# formula together — no normalization, dispersion estimation, or
# testing has happened yet. Specifically:
#   - assay(dds)        -> the raw count matrix (unchanged from Gene_matrix)
#   - colData(dds)       -> pheno_data, with "Day" stored as the design variable
#   - design(dds)        -> the ~ Day formula
#   - counts(dds)        -> another way to access the raw counts
# Running DESeq(dds) later is what actually performs size factor
# estimation, dispersion estimation, and the negative binomial GLM fitting.
# Until you run DESeq(dds), results(dds) will not work.

dds <- DESeqDataSetFromMatrix(Gene_matrix, pheno_data, design = ~ Day)

# Filter genes
# Keep genes that have at least 10 counts across the smallest group size
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dim(dds)
#[1] 15922     9

#############################################
# Run DESeq2
#############################################
# DESeq(dds) runs the full analysis pipeline on the DESeqDataSet in three steps:
#
# 1. estimateSizeFactors(dds)
#    - calculates per-sample normalization factors (median-of-ratios method)
#    - accounts for differences in sequencing depth/library size between samples
#    - stored in sizeFactors(dds)
#
# 2. estimateDispersions(dds)
#    - estimates gene-wise dispersion (variance beyond what the mean predicts)
#    - fits a dispersion trend across all genes, then shrinks each gene's
#      estimate toward that trend (empirical Bayes shrinkage)
#    - this shrinkage is what makes DESeq2 robust with small sample sizes
#    - stored in dispersions(dds) / accessible via plotDispEsts(dds)
#
# 3. nbinomWaldTest(dds)
#    - fits a negative binomial GLM per gene using the design formula
#    - tests coefficients via Wald tests
#    - produces the log2 fold changes and p-values that results() will extract
#
# Output: dds is still a DESeqDataSet, but now fully populated with the
# fitted model. Nothing is normalized/transformed for visualization yet
# (that's what counts(dds, normalized=TRUE) or vst()/rlog() are for).
# What you now have access to:
#   - sizeFactors(dds)     -> normalization factors per sample
#   - dispersions(dds)     -> final (shrunk) per-gene dispersion estimates
#   - resultsNames(dds)    -> names of the coefficients you can pull with results()
#                             (e.g. "Day_T0_vs_C", "Day_TM_vs_C", etc.,
#                             depending on your reference level)
#   - results(dds, ...)    -> NOW usable; this is where you actually pull
#                             log2FoldChange, p-value, and padj per gene,
#                             specifying a contrast or coefficient name
dds <- DESeq(dds)

# save object for students that can't install DESeq2
saveRDS(dds, file = "DESeq_results_object.rds")



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
#
# results(dds, contrast = ...) is how you pull specific pairwise
# comparisons out of the fitted dds object. Internally it reads from
# the per-gene stats stored in rowRanges(dds)'s elementMetadata (the
# 26-column DFrame you saw earlier) — but you don't access that
# directly; results() handles the correct extraction and computes
# padj (BH-adjusted p-values) fresh for each specific contrast.
#
# contrast = c("factor_name", "level_of_interest", "reference_level")
#   - "Day"  -> the design variable being tested
#   - "D7"   -> numerator (the level of interest)
#   - "D0"   -> denominator (the reference/baseline level)
#   - log2FoldChange in the output is log2(D7 / D0) — positive values
#     mean higher expression at D7 relative to D0
#
# Each results() call returns its own DESeqResults object (one row
# per gene) with columns:
#   baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
#
# Note: this approach re-derives the contrast directly from the
# fitted model coefficients (via the design formula), so it works
# regardless of which level was set as the reference when the factor
# was created — you don't need "Day" releveled to D0 beforehand for
# this to give the correct D7-vs-D0 comparison.

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

