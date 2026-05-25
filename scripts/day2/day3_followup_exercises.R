# =============================================================================
# Day 3 Follow-up Exercises — Microbial Taxonomic and Functional Ecology
# =============================================================================
# This script assumes the data loading, reordering, normalization, and alpha /
# beta diversity sections from the main Day 3 script have already been run.
# The following objects must already be in your environment:
#   md                   — metadata data frame (with urban, site columns)
#   taxonomy_counts      — raw taxonomy count table (samples as rows)
#   ko_counts            — raw KO count table (samples as rows)
#   rpm_taxonomy_counts  — RPM-normalised taxonomy table (samples as rows)
#   rpm_ko_counts        — RPM-normalised KO table (samples as rows)
#   taxonomy.bray        — Bray-Curtis distance matrix for taxonomy
#   ko.bray              — Bray-Curtis distance matrix for KOs
#   taxonomy.nmds        — metaMDS object for taxonomy
#   ko.nmds              — metaMDS object for KOs
# =============================================================================


# -----------------------------------------------------------------------------
# Install and load additional packages
# -----------------------------------------------------------------------------

# CRAN packages — only run once
install.packages(c("indicspecies", "cooccur", "igraph", "betapart", "ggrepel"))


# Load everything
library(tidyverse)
library(vegan)
library(indicspecies)   # indicator species analysis
library(betapart)       # nestedness vs. turnover decomposition
library(cooccur)        # probabilistic co-occurrence
library(igraph)         # network analysis
library(ggpubr)         # stat_cor() for correlation plots
library(ggrepel)        # non-overlapping labels on volcano plots




# =============================================================================
# BEGINNER
# =============================================================================


# -----------------------------------------------------------------------------
# Exercise 1 — Rarefaction curves
# -----------------------------------------------------------------------------
# Goal: visualise how richness accumulates with sequencing depth and decide
# whether the rarefaction cutoff used in the main script is justified.
# Try it: change step to a smaller number — how does curve resolution change?
#         change sample to a lower value — how does richness change?

taxonomy_counts_forrich <- taxonomy_counts %>%
  column_to_rownames(var = "sample")

ko_counts_forrich <- ko_counts %>%
  column_to_rownames(var = "sample")

# Colour lines by habitat
tax_colors <- ifelse(md$urban == "urban", "steelblue", "forestgreen")

# Taxonomic rarefaction
rarecurve(
  taxonomy_counts_forrich,
  step   = 100000,     # << try a smaller number, e.g. 50000
  sample = 7734272,    # << try a lower depth and see how richness estimates change
  col    = tax_colors,
  xlab   = "Sequencing depth (reads)",
  ylab   = "Number of taxa",
  main   = "Taxonomic rarefaction curves",
  label  = FALSE
)
legend("bottomright",
       legend = c("Urban", "Natural"),
       col    = c("steelblue", "forestgreen"),
       lty    = 1)

# Functional (KO) rarefaction
rarecurve(
  ko_counts_forrich,
  step   = 50000,      # << try a smaller number
  sample = 2352294,    # << try a lower depth
  col    = tax_colors,
  xlab   = "Sequencing depth (reads)",
  ylab   = "Number of KOs",
  main   = "Functional rarefaction curves",
  label  = FALSE
)
legend("bottomright",
       legend = c("Urban", "Natural"),
       col    = c("steelblue", "forestgreen"),
       lty    = 1)

# Are curves flattening near the chosen depth? If they are still rising steeply,
# richness estimates are sensitive to sequencing effort and the cutoff may be
# too conservative.




# -----------------------------------------------------------------------------
# Exercise 2 — Stacked bar chart of top N taxa
# -----------------------------------------------------------------------------
# Goal: compare which taxa dominate urban vs. natural communities.
# Try it: change N_TAXA from 20 to 10 or 50 — does the biological story change?

N_TAXA <- 20   # << change this

# Identify the N most abundant taxa (by mean RPM across all samples)
top_taxa <- rpm_taxonomy_counts %>%
  pivot_longer(-sample, names_to = "taxon", values_to = "RPM") %>%
  group_by(taxon) %>%
  summarise(mean_RPM = mean(RPM)) %>%
  slice_max(mean_RPM, n = N_TAXA) %>%
  pull(taxon)

# Collapse everything else into "Other"
stacked_data <- rpm_taxonomy_counts %>%
  pivot_longer(-sample, names_to = "taxon", values_to = "RPM") %>%
  left_join(md %>% select(sample, urban, site), by = "sample") %>%
  mutate(taxon = ifelse(taxon %in% top_taxa, taxon, "Other")) %>%
  group_by(sample, urban, site, taxon) %>%
  summarise(RPM = sum(RPM), .groups = "drop")

ggplot(stacked_data, aes(x = sample, y = RPM, fill = taxon)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ urban, scales = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x     = NULL,
    y     = "Relative abundance",
    fill  = "Taxon",
    title = paste("Top", N_TAXA, "taxa by mean RPM")
  ) +
  theme_bw() +
  theme(
    axis.text.x    = element_text(angle = 90, hjust = 1, size = 6),
    legend.key.size = unit(0.4, "cm")
  )

# Do the dominant taxa shift between urban and natural?
# Are there taxa that appear almost exclusively in one habitat type?




# -----------------------------------------------------------------------------
# Exercise 3 — Indicator species analysis
# -----------------------------------------------------------------------------
# Goal: statistically identify taxa (and KOs) that are most diagnostic of
# urban vs. natural sites using the indicator value (IndVal) metric.
# Try it: change alpha = 0.05 to 0.01 — how many indicators do you lose?

groups <- ifelse(md$urban == "urban", 1, 2)  # 1 = urban, 2 = natural

# Taxonomy indicators
taxonomy_matrix <- taxonomy_counts %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()

indval_taxonomy <- multipatt(
  taxonomy_matrix,
  groups,
  func    = "IndVal.g",    # corrected IndVal for unequal group sizes
  control = how(nperm = 999)
)

summary(indval_taxonomy, indvalcomp = TRUE, alpha = 0.05)  # << change alpha here

# Extract results as a tidy data frame for downstream use
indval_results <- indval_taxonomy$sign %>%
  as.data.frame() %>%
  rownames_to_column("taxon") %>%
  filter(p.value < 0.05) %>%   # << change this threshold to match alpha above
  arrange(p.value)

cat("Number of significant taxonomy indicators:", nrow(indval_results), "\n")

# KO indicators
ko_matrix <- ko_counts %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()

indval_ko <- multipatt(
  ko_matrix,
  groups,
  func    = "IndVal.g",
  control = how(nperm = 999)
)

summary(indval_ko, indvalcomp = TRUE, alpha = 0.05)




# =============================================================================
# INTERMEDIATE
# =============================================================================


# -----------------------------------------------------------------------------
# Exercise 4 — Differential abundance testing with DESeq2
# -----------------------------------------------------------------------------
# Goal: identify specific taxa/KOs that are significantly more or less abundant
# in urban vs. natural sites, controlling for multiple testing.
# Try it: change FC_THRESHOLD from 1 to 0.5 or 2 — how many taxa change?
#         change padj < 0.05 to padj < 0.01 — how many taxa do you lose?
#         the ggrepel labels are already included — try changing n = 10

# -- Taxonomy

tax_counts_mat <- taxonomy_counts %>%
  column_to_rownames(var = "sample") %>%
  t() %>%
  round() %>%         # DESeq2 requires integer counts
  as.matrix()

col_data <- md %>%
  select(sample, urban, site) %>%
  column_to_rownames("sample")

dds_tax <- DESeqDataSetFromMatrix(
  countData = tax_counts_mat,
  colData   = col_data,
  design    = ~ urban   # add + site to control for site effects
)
dds_tax <- dds_tax[rowSums(counts(dds_tax)) >= 10, ]  # remove very rare taxa
dds_tax <- DESeq(dds_tax)

res_tax <- results(
  dds_tax,
  contrast = c("urban", "urban", "natural"),  # positive FC = higher in urban
  alpha    = 0.05
)
summary(res_tax)

res_tax_df <- as.data.frame(res_tax) %>%
  rownames_to_column("taxon") %>%
  filter(!is.na(padj))

# Count how many taxa fall into each category
FC_THRESHOLD <- 1   # << change this: try 0.5 or 2

res_tax_df %>%
  mutate(direction = case_when(
    padj < 0.05 & log2FoldChange >  FC_THRESHOLD ~ "Urban-enriched",
    padj < 0.05 & log2FoldChange < -FC_THRESHOLD ~ "Natural-enriched",
    TRUE                                          ~ "Not significant"
  )) %>%
  count(direction)

# Labels for the 10 most significant taxa
top10_labels <- res_tax_df %>%
  filter(padj < 0.05) %>%
  slice_min(padj, n = 10)

# Volcano plot
res_tax_df %>%
  mutate(significant = padj < 0.05 & abs(log2FoldChange) > FC_THRESHOLD) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-FC_THRESHOLD, FC_THRESHOLD), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_label_repel(
    data         = top10_labels,
    aes(label    = taxon),
    color        = "black",
    size         = 3,
    max.overlaps = 20
  ) +
  scale_color_manual(values = c("grey60", "firebrick")) +
  labs(
    title = "Differential abundance: Urban vs. Natural (taxonomy)",
    x     = "log2 fold change",
    y     = "-log10 adjusted p-value"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# Save significant taxa to CSV
res_tax_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > FC_THRESHOLD) %>%
  arrange(padj) %>%
  write_csv("significant_taxa_deseq2.csv")

# -- KOs (same workflow)

ko_counts_mat <- ko_counts %>%
  column_to_rownames(var = "sample") %>%
  t() %>%
  round() %>%
  as.matrix()

dds_ko <- DESeqDataSetFromMatrix(
  countData = ko_counts_mat,
  colData   = col_data,
  design    = ~ urban
)
dds_ko <- dds_ko[rowSums(counts(dds_ko)) >= 10, ]
dds_ko <- DESeq(dds_ko)

res_ko <- results(dds_ko, contrast = c("urban", "urban", "natural"), alpha = 0.05)
summary(res_ko)

res_ko_df <- as.data.frame(res_ko) %>%
  rownames_to_column("KO") %>%
  filter(!is.na(padj))

res_ko_df %>%
  mutate(direction = case_when(
    padj < 0.05 & log2FoldChange >  FC_THRESHOLD ~ "Urban-enriched",
    padj < 0.05 & log2FoldChange < -FC_THRESHOLD ~ "Natural-enriched",
    TRUE                                          ~ "Not significant"
  )) %>%
  count(direction)

res_ko_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > FC_THRESHOLD) %>%
  arrange(padj) %>%
  write_csv("significant_KOs_deseq2.csv")




# -----------------------------------------------------------------------------
# Exercise 5 — Environmental gradient modeling with envfit
# -----------------------------------------------------------------------------
# Goal: test which continuous environmental variables (e.g. pH, temperature,
# % impervious surface) are significantly associated with the ordination axes,
# then overlay them as vectors on the NMDS plot.
# Note: replace the column names below with whatever continuous variables
#       you actually have in your metadata.
# Try it: colour the NMDS by a continuous variable instead of urban/natural —
#         does any gradient visually explain the ordination axes?

# Select only the continuous environmental columns from md
# (adjust this selection to match your actual metadata columns)
env_vars <- md %>%
  select(where(is.numeric)) %>%
  select(-ends_with("Axis01"), -ends_with("Axis02"),
         -taxonomic_rich, -taxonomic_shannon,
         -functional_rich, -functional_shannon,
         -beta_turnover, -beta_nestedness, -beta_total)  # exclude columns we added

# -- Taxonomy NMDS
envfit_tax <- envfit(
  taxonomy.nmds,
  env_vars,
  permutations = 999,
  na.rm        = TRUE
)
print(envfit_tax)   # r² and p-value for each variable

# Plot with significant vectors overlaid (p < 0.05)
plot(taxonomy.nmds, display = "sites", type = "n",
     main = "Taxonomy NMDS + environmental vectors")
points(taxonomy.nmds, display = "sites",
       col = ifelse(md$urban == "urban", "steelblue", "forestgreen"),
       pch = 16, cex = 1.5)
plot(envfit_tax, p.max = 0.05, col = "black", cex = 0.8)
legend("topright",
       legend = c("Urban", "Natural"),
       col    = c("steelblue", "forestgreen"),
       pch    = 16)

# Smooth surface for a single variable (uncomment and replace column name)
# ordisurf(taxonomy.nmds, md$your_variable, main = "Gradient: your_variable")

# -- KO NMDS
envfit_ko <- envfit(ko.nmds, env_vars, permutations = 999, na.rm = TRUE)
print(envfit_ko)

plot(ko.nmds, display = "sites", type = "n",
     main = "KO NMDS + environmental vectors")
points(ko.nmds, display = "sites",
       col = ifelse(md$urban == "urban", "steelblue", "forestgreen"),
       pch = 16, cex = 1.5)
plot(envfit_ko, p.max = 0.05, col = "black", cex = 0.8)
legend("topright",
       legend = c("Urban", "Natural"),
       col    = c("steelblue", "forestgreen"),
       pch    = 16)

# Does the same variable drive both taxonomic and functional composition, or
# do different gradients explain each?




# -----------------------------------------------------------------------------
# Exercise 6 — Variance partitioning
# -----------------------------------------------------------------------------
# Goal: quantify how much compositional variation is explained by urbanisation
# alone, site identity alone, and their shared (confounded) effect.
# Try it: add a third continuous variable from env_vars and re-run — does the
#         unexplained fraction shrink?

# Hellinger transformation is required for varpart (distance-based equivalent)
taxonomy_hell <- decostand(
  rpm_taxonomy_counts %>% column_to_rownames("sample"),
  method = "hellinger"
)

ko_hell <- decostand(
  rpm_ko_counts %>% column_to_rownames("sample"),
  method = "hellinger"
)

# Build explanatory matrices (dummy-code factors, drop intercept)
X1_urban <- model.matrix(~ urban, data = md)[, -1, drop = FALSE]
X2_site  <- model.matrix(~ site,  data = md)[, -1, drop = FALSE]

# -- Taxonomy
vp_tax <- varpart(taxonomy_hell, X1_urban, X2_site)
print(vp_tax)

plot(vp_tax,
     Xnames = c("Urbanisation", "Site"),
     main   = "Variance partitioning — Taxonomy",
     bg     = c("steelblue", "forestgreen"))

# Test significance of each fraction using partial RDA
# [a] urbanisation controlling for site
anova(rda(taxonomy_hell ~ urban + Condition(site), data = md),
      permutations = 999)

# [b] site controlling for urbanisation
anova(rda(taxonomy_hell ~ site + Condition(urban), data = md),
      permutations = 999)

# -- KOs
vp_ko <- varpart(ko_hell, X1_urban, X2_site)
print(vp_ko)

plot(vp_ko,
     Xnames = c("Urbanisation", "Site"),
     main   = "Variance partitioning — KOs",
     bg     = c("steelblue", "forestgreen"))

anova(rda(ko_hell ~ urban + Condition(site), data = md), permutations = 999)
anova(rda(ko_hell ~ site  + Condition(urban), data = md), permutations = 999)

# If the shared fraction [c] is large, urban and site effects are confounded:
# urban communities differ partly because of geography, not just urbanisation.




# =============================================================================
# ADVANCED
# =============================================================================


# -----------------------------------------------------------------------------
# Exercise 7 — Nestedness vs. turnover decomposition
# -----------------------------------------------------------------------------
# Goal: decompose beta diversity into two components:
#   Turnover   — species are replaced by different species across sites
#   Nestedness — poorer sites are subsets of richer sites (species loss)
# This tells you whether urban/natural communities differ because they have
# completely different taxa (turnover) or because urban sites lost taxa that
# natural sites still have (nestedness).

# Convert RPM to presence/absence
taxonomy_pa <- rpm_taxonomy_counts %>%
  column_to_rownames("sample") %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))

# Compute the betapart core object (only needs to be done once)
beta_core <- betapart.core(taxonomy_pa)

# Decompose into pairwise turnover, nestedness, and total Sorensen dissimilarity
beta_pair <- beta.pair(beta_core, index.family = "sorensen")
# beta_pair$beta.sim  = turnover (Simpson dissimilarity)
# beta_pair$beta.sne  = nestedness-resultant dissimilarity
# beta_pair$beta.sor  = total Sorensen beta diversity (sim + sne)

# Attach mean values per sample to md for plotting
md$beta_turnover   <- rowMeans(as.matrix(beta_pair$beta.sim))
md$beta_nestedness <- rowMeans(as.matrix(beta_pair$beta.sne))
md$beta_total      <- rowMeans(as.matrix(beta_pair$beta.sor))

# Which component dominates overall?
cat("Overall mean turnover:  ",
    mean(as.matrix(beta_pair$beta.sim)[upper.tri(as.matrix(beta_pair$beta.sim))]),
    "\n")
cat("Overall mean nestedness:",
    mean(as.matrix(beta_pair$beta.sne)[upper.tri(as.matrix(beta_pair$beta.sne))]),
    "\n")

# Plot turnover vs. nestedness by habitat
md %>%
  select(sample, urban, beta_turnover, beta_nestedness) %>%
  pivot_longer(c(beta_turnover, beta_nestedness),
               names_to  = "component",
               values_to = "dissimilarity") %>%
  mutate(component = recode(component,
                            beta_turnover   = "Turnover",
                            beta_nestedness = "Nestedness")) %>%
  ggplot(aes(x = urban, y = dissimilarity, color = urban)) +
  geom_jitter(width = 0.15) +
  geom_boxplot(outliers = FALSE, alpha = 0.4) +
  facet_wrap(~ component) +
  scale_color_manual(values = c("urban" = "steelblue", "natural" = "forestgreen")) +
  labs(x = NULL, y = "Mean pairwise dissimilarity",
       title = "Beta diversity: nestedness vs. turnover") +
  theme_bw() +
  theme(legend.position = "none")

# Test whether the components differ between habitats
anova(lm(beta_turnover   ~ urban, data = md))
anova(lm(beta_nestedness ~ urban, data = md))

# Multi-site decomposition per habitat (single summary value per group)
urban_pa_beta   <- taxonomy_pa[md$urban == "urban",   ]
natural_pa_beta <- taxonomy_pa[md$urban == "natural", ]

cat("\nUrban multi-site beta diversity:\n")
print(beta.multi(betapart.core(urban_pa_beta),   index.family = "sorensen"))
cat("\nNatural multi-site beta diversity:\n")
print(beta.multi(betapart.core(natural_pa_beta), index.family = "sorensen"))

# If turnover dominates: urban and natural sites have fundamentally different
# microbial assemblages. If nestedness dominates: urbanisation is filtering out
# taxa that natural sites still retain.




# -----------------------------------------------------------------------------
# Exercise 8 — Mantel test
# -----------------------------------------------------------------------------
# Goal: test whether communities that are taxonomically similar are also
# functionally similar, by directly correlating two distance matrices.

mantel_result <- mantel(
  taxonomy.bray,
  ko.bray,
  method       = "pearson",   # use "spearman" for a rank-based version
  permutations = 999
)
print(mantel_result)
# r  = Mantel statistic (strength of the matrix correlation)
# p  = significance estimated from permutations

# Scatter plot of pairwise distances
tax_dist_vec <- as.vector(as.dist(taxonomy.bray))
ko_dist_vec  <- as.vector(as.dist(ko.bray))

ggplot(data.frame(tax = tax_dist_vec, ko = ko_dist_vec),
       aes(x = tax, y = ko)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", color = "firebrick") +
  stat_cor(method = "pearson") +
  labs(
    x     = "Taxonomic Bray-Curtis dissimilarity",
    y     = "Functional Bray-Curtis dissimilarity",
    title = "Mantel scatter plot — all samples"
  ) +
  theme_bw()

# A high r means taxonomy predicts function well (low functional redundancy).
# A low r means different taxa perform the same functions across sites
# (high functional redundancy) — common in diverse soil microbiomes.




# -----------------------------------------------------------------------------
# Exercise 9 — Taxonomic-functional decoupling
# -----------------------------------------------------------------------------
# Goal: extend the Mantel test by asking whether the taxonomy-function
# relationship differs between urban and natural samples.
# If r is lower in urban sites, fewer taxa may be performing more of the
# functional work there (reduced functional redundancy).

urban_idx   <- which(md$urban == "urban")
natural_idx <- which(md$urban == "natural")

tax_urban   <- as.dist(as.matrix(taxonomy.bray)[urban_idx,   urban_idx])
ko_urban    <- as.dist(as.matrix(ko.bray)      [urban_idx,   urban_idx])
tax_natural <- as.dist(as.matrix(taxonomy.bray)[natural_idx, natural_idx])
ko_natural  <- as.dist(as.matrix(ko.bray)      [natural_idx, natural_idx])

mantel_urban   <- mantel(tax_urban,   ko_urban,   method = "pearson", permutations = 999)
mantel_natural <- mantel(tax_natural, ko_natural, method = "pearson", permutations = 999)

cat("Mantel r — Urban:  ", round(mantel_urban$statistic,   3),
    "  p =", mantel_urban$signif,   "\n")
cat("Mantel r — Natural:", round(mantel_natural$statistic, 3),
    "  p =", mantel_natural$signif, "\n")

bind_rows(
  data.frame(tax = as.vector(tax_urban),   ko = as.vector(ko_urban),   group = "Urban"),
  data.frame(tax = as.vector(tax_natural), ko = as.vector(ko_natural), group = "Natural")
) %>%
  ggplot(aes(x = tax, y = ko, color = group)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  scale_color_manual(values = c("Urban" = "steelblue", "Natural" = "forestgreen")) +
  facet_wrap(~ group) +
  labs(
    x     = "Taxonomic Bray-Curtis dissimilarity",
    y     = "Functional Bray-Curtis dissimilarity",
    title = "Taxonomy-function relationship by habitat"
  ) +
  theme_bw() +
  theme(legend.position = "none")




# -----------------------------------------------------------------------------
# Exercise 10 — Co-occurrence networks
# -----------------------------------------------------------------------------
# Goal: build co-occurrence networks for urban and natural communities and
# compare their structure — are urban networks more or less connected?
# Note: computation scales as N², so we limit to the top 20 most abundant
#       taxa by default. Increase N_NETWORK cautiously (30 is still fast,
#       50+ may take a few minutes).
# Try it: change N_NETWORK to 30 and see how the network changes.

N_NETWORK <- 20   # << change this cautiously — 30 is reasonable, 50+ is slow

top_network_taxa <- rpm_taxonomy_counts %>%
  pivot_longer(-sample, names_to = "taxon", values_to = "RPM") %>%
  group_by(taxon) %>%
  summarise(mean_RPM = mean(RPM)) %>%
  slice_max(mean_RPM, n = N_NETWORK) %>%
  pull(taxon)

# Helper: build a presence/absence matrix for one habitat (taxa as rows)
make_pa_matrix <- function(habitat) {
  rpm_taxonomy_counts %>%
    filter(sample %in% md$sample[md$urban == habitat]) %>%
    select(all_of(c("sample", top_network_taxa))) %>%
    column_to_rownames("sample") %>%
    mutate(across(everything(), ~ ifelse(. > 0, 1, 0))) %>%
    t() %>%
    as.data.frame()
}

urban_pa_net   <- make_pa_matrix("urban")
natural_pa_net <- make_pa_matrix("natural")

# Run probabilistic co-occurrence analysis
cooccur_urban   <- cooccur(urban_pa_net,   type = "spp_site", spp_names = TRUE)
cooccur_natural <- cooccur(natural_pa_net, type = "spp_site", spp_names = TRUE)

summary(cooccur_urban)
summary(cooccur_natural)

# Helper: convert significant co-occurrences to an igraph object
make_cooccur_graph <- function(cooccur_obj) {
  pairs <- prob.table(cooccur_obj) %>%
    filter(p_lt < 0.05 | p_gt < 0.05) %>%
    mutate(weight = ifelse(p_gt < 0.05, 1, -1))  # +1 positive, -1 negative
  if (nrow(pairs) == 0) {
    warning("No significant co-occurrences found — try increasing N_NETWORK")
    return(make_empty_graph())
  }
  graph_from_data_frame(
    pairs[, c("sp1_name", "sp2_name", "weight")],
    directed = FALSE
  )
}

g_urban   <- make_cooccur_graph(cooccur_urban)
g_natural <- make_cooccur_graph(cooccur_natural)

# Compare network metrics
cat("Urban network:   nodes =", vcount(g_urban),   " edges =", ecount(g_urban),
    " connectance =", round(edge_density(g_urban),   3), "\n")
cat("Natural network: nodes =", vcount(g_natural), " edges =", ecount(g_natural),
    " connectance =", round(edge_density(g_natural), 3), "\n")

# Plot both networks side by side
# Blue edges = positive co-occurrence, red = mutual exclusion
par(mfrow = c(1, 2))

plot(g_urban,
     vertex.size      = degree(g_urban) * 3,
     vertex.label     = V(g_urban)$name,
     vertex.label.cex = 0.6,
     edge.color       = ifelse(E(g_urban)$weight > 0, "steelblue", "firebrick"),
     main             = paste0("Urban (top ", N_NETWORK, " taxa)"))

plot(g_natural,
     vertex.size      = degree(g_natural) * 3,
     vertex.label     = V(g_natural)$name,
     vertex.label.cex = 0.6,
     edge.color       = ifelse(E(g_natural)$weight > 0, "steelblue", "firebrick"),
     main             = paste0("Natural (top ", N_NETWORK, " taxa)"))

par(mfrow = c(1, 1))

# Hub taxa (highest number of significant co-occurrences)
cat("\nTop 5 hub taxa — Urban:\n")
print(sort(degree(g_urban),   decreasing = TRUE)[1:5])
cat("\nTop 5 hub taxa — Natural:\n")
print(sort(degree(g_natural), decreasing = TRUE)[1:5])

# Are urban networks more or less connected than natural ones?
# More connected networks can indicate tighter species interactions, or
# strong environmental filtering forcing similar taxa to always co-occur.


# =============================================================================
# END OF EXERCISES
# =============================================================================
