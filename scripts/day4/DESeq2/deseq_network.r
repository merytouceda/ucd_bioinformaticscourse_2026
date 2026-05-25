# UCD - bioinformatics summer course
# Week 2 - Day 2: RNA-seq
# Gene co-expression networks with WGCNA

# ==============================
# LIBRARY AND WORKING DIR SET UP
# ==============================

#BiocManager::install("DESeq2", force = TRUE)
library(DESeq2)
#BiocManager::install("preprocessCore")
#BiocManager::install("impute")
#BiocManager::install("WGCNA")
library(WGCNA)
library(tibble)
library(tidyverse)

setwd("~/Documents/GitHub/ucd_bioinformaticscourse_2026/scripts/day2/Sheep_RNA_seq_for_Mery/DESeq2/")



# ==============================
# READ IN DATA
# ==============================

Gene_matrix <- read.table("Gene_matrix.txt", sep="\t", header=TRUE)
dim(Gene_matrix)
#[1] 28524     9

pheno_data <- read.table("Metadata.txt", sep="\t", header=TRUE) %>%
  mutate(Day = factor(Day)) %>%
  column_to_rownames("Sample_ID") %>%
  mutate(Sample_ID = rownames(.))   # NOTE: fixed original bug (was rownames(pheno_data) before assignment)
dim(pheno_data)
#[1] 9 3

levels(pheno_data$Day)
#[1] "D0"  "D14" "D7"

all(rownames(pheno_data) == colnames(Gene_matrix))
#[1] TRUE

# ==============================
# DESeq2
# ==============================

dds <- DESeqDataSetFromMatrix(Gene_matrix, pheno_data, design = ~ Day)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds)
#[1] 15922     9

dds <- DESeq(dds)

# ==============================
# PCA
# ==============================

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup="Day", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaPlot <- ggplot(pcaData, aes(x=PC1, y=PC2, color=Day)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

pcaPlot
ggsave("PCA_plot.pdf", pcaPlot, device="pdf")

# ==============================
# WGCNA — CO-EXPRESSION NETWORK
# ==============================

# WGCNA requires multi-threading to be enabled
allowWGCNAThreads()

# --- Prepare expression matrix ---
# Reuse the VST-normalized matrix from DESeq2
# WGCNA expects samples as rows, genes as columns
datExpr <- t(assay(vsd))

# QC: check for genes or samples with too many missing values
# (unlikely with count data, but good practice to show students)
gsg <- goodSamplesGenes(datExpr, verbose=3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  message("Removed ", sum(!gsg$goodGenes), " genes and ",
          sum(!gsg$goodSamples), " samples after QC")
}

# --- Pick soft-thresholding power ---
# WGCNA transforms Pearson correlations to a weighted adjacency by raising
# them to a power β. We pick β so the network approximates scale-free topology.
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=5)

# Plot 1: Scale-free topology fit index
# — look for the lowest power where R² first plateaus above 0.80
sft_plot <- ggplot(data.frame(Power   = sft$fitIndices[,1],
                              SFT.R2 = -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]),
                   aes(x=Power, y=SFT.R2, label=Power)) +
  geom_point(color="steelblue") +
  geom_text(vjust=-0.5, size=3) +
  geom_hline(yintercept=0.80, linetype="dashed", color="red") +
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology R²") +
  ggtitle("Scale-free topology fit — choose power above red line") +
  theme_bw()

sft_plot
ggsave("WGCNA_soft_threshold.pdf", sft_plot, device="pdf")

# Plot 2: Mean connectivity
# — as a sanity check; too-high connectivity suggests power is too low
conn_plot <- ggplot(data.frame(Power       = sft$fitIndices[,1],
                               MeanK      = sft$fitIndices[,5]),
                    aes(x=Power, y=MeanK, label=Power)) +
  geom_point(color="darkorange") +
  geom_text(vjust=-0.5, size=3) +
  xlab("Soft Threshold (power)") +
  ylab("Mean Connectivity") +
  ggtitle("Mean connectivity vs power") +
  theme_bw()

conn_plot
ggsave("WGCNA_mean_connectivity.pdf", conn_plot, device="pdf")

# Set power — use the estimate from pickSoftThreshold, or override manually
# NOTE: with only 9 samples the automated estimate may be unreliable;
# students should look at the scale-free plot and choose the elbow.
soft_power <- sft$powerEstimate
message("Suggested soft-thresholding power: ", soft_power)

# --- Build network and detect modules ---
# blockwiseModules() computes the TOM (Topological Overlap Matrix),
# clusters genes by TOM-dissimilarity, and assigns each cluster a color label.
#
# Key parameters students should understand:
#   minModuleSize  — minimum number of genes per module (use small value for
#                    teaching; 30 is typical in real studies)
#   mergeCutHeight — threshold for merging similar modules (0 = no merging,
#                    1 = merge everything); 0.25 merges modules with >75%
#                    eigengene correlation
net <- blockwiseModules(
  datExpr,
  power            = 4,
  TOMType          = "unsigned",
  minModuleSize    = 5,          # low for teaching purposes
  mergeCutHeight   = 0.25,
  numericLabels    = FALSE,      # use color names, not numbers
  pamRespectsDendro = FALSE,
  saveTOMs         = FALSE,
  verbose          = 3
)

# How many modules, and how many genes per module?
table(net$colors)

# --- Visualize the gene dendrogram with module colors ---
pdf("WGCNA_dendrogram.pdf", width=10, height=5)
plotDendroAndColors(
  net$dendrograms[[1]],
  net$colors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang         = 0.03,
  addGuide     = TRUE,
  guideHang    = 0.05,
  main         = "Gene dendrogram and module assignment"
)
dev.off()

# --- Module eigengenes ---
# The eigengene is the first principal component of a module — it summarises
# the overall expression pattern of all genes in that module across samples.
MEs <- net$MEs

# Add sample metadata for plotting
MEs_plot <- MEs %>%
  rownames_to_column("Sample_ID") %>%
  left_join(pheno_data, by="Sample_ID") %>%          # pheno_data already has Sample_ID
  pivot_longer(cols=starts_with("ME"), names_to="Module", values_to="Eigengene") %>%
  mutate(Module = sub("^ME", "", Module))

# Plot eigengene expression across time points per module
eigengene_plot <- ggplot(MEs_plot, aes(x=Day, y=Eigengene, color=Day, group=Day)) +
  geom_jitter(width=0.1, size=2) +
  stat_summary(fun=mean, geom="crossbar", width=0.4, color="black") +
  facet_wrap(~Module, scales="free_y") +
  xlab("Day") +
  ylab("Module eigengene") +
  ggtitle("Module eigengene expression by time point") +
  theme_bw() +
  theme(legend.position="none")

eigengene_plot
ggsave("WGCNA_eigengenes.pdf", eigengene_plot, width=10, height=8, device="pdf")

# --- Module–trait correlation heatmap ---
# Correlate each module eigengene with numeric representations of the trait
# (here: days post-infection as a continuous variable)
trait_matrix <- data.frame(
  D0  = as.integer(pheno_data$Day == "D0"),
  D7  = as.integer(pheno_data$Day == "D7"),
  D14 = as.integer(pheno_data$Day == "D14")
)
rownames(trait_matrix) <- rownames(pheno_data)

# Pearson correlation + p-values
module_trait_cor <- cor(MEs, trait_matrix, use="pairwise.complete.obs")
module_trait_pval <- corPvalueStudent(module_trait_cor, nSamples=nrow(datExpr))

# Format labels as "r\n(p)" for the heatmap cells
text_matrix <- paste0(
  round(module_trait_cor, 2), "\n(",
  signif(module_trait_pval, 1), ")"
)
dim(text_matrix) <- dim(module_trait_cor)

pdf("WGCNA_module_trait_heatmap.pdf", width=15, height=12)
labeledHeatmap(
  Matrix     = module_trait_cor,
  xLabels    = colnames(module_trait_cor),
  yLabels    = rownames(module_trait_cor),
  ySymbols   = rownames(module_trait_cor),
  colorLabels = FALSE,
  colors     = blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = FALSE,
  cex.text   = 0.7,
  zlim       = c(-1, 1),
  main       = "Module–trait correlations\n(r value, p-value)"
)
dev.off()

# --- Export hub genes for a module of interest ---
# Hub genes = most highly connected genes within a module.
# Students can change "turquoise" to whichever module they want to explore.
module_of_interest <- "turquoise"

# Module membership (kME): correlation of each gene's expression with the eigengene
kME <- signedKME(datExpr, MEs)

hub_genes <- kME %>%
  rownames_to_column("Gene") %>%
  select(Gene, kME = paste0("kME", module_of_interest)) %>%
  filter(Gene %in% names(net$colors)[net$colors == module_of_interest]) %>%
  arrange(desc(kME)) %>%
  head(20)

hub_genes
write.csv(hub_genes, paste0("WGCNA_hub_genes_", module_of_interest, ".csv"),
          row.names=FALSE)
