# Extra questions: solutions
# Day 3 Exercise 1 follow-up — UCD Bioinformatics Course
# Author: Mery Touceda-Suarez
# Date: June 2026

library(tidyverse)

phage_inferences <- read_csv("~/Documents/GitHub/ucd_bioinformaticscourse_2026/data/phages_summary.csv") %>%
  mutate(type = case_when(str_detect(votu, "provirus") ~ "provirus",
                           TRUE ~ "virus"))

# Earlier simple filter, used for comparison later on
phage_inferences_filtered <- phage_inferences %>%
  filter(contig_length > 5000)


# ==========================================================================
# Q1. Complete the code: contamination and completeness
# ==========================================================================

# --------------- Completeness
hist(phage_inferences$completeness)
summary(phage_inferences$completeness)

ggplot(phage_inferences, aes(x = completeness)) +
  geom_histogram() +
  facet_wrap(~type) +
  theme_classic() +
  ggtitle("Completeness distribution by type")

# --------------- Contamination
hist(phage_inferences$contamination)
summary(phage_inferences$contamination)

ggplot(phage_inferences, aes(x = contamination)) +
  geom_histogram() +
  facet_wrap(~type) +
  theme_classic() +
  ggtitle("Contamination distribution by type")

# --------------- checkv_quality categories (CheckV's own summary classification)
phage_inferences %>%
  count(checkv_quality) %>%
  arrange(desc(n))

ggplot(phage_inferences, aes(x = checkv_quality)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Number of contigs per CheckV quality tier")


# ==========================================================================
# Q2. Trade-offs of filtering — when to be more/less strict?
# ==========================================================================

# This is a discussion question rather than a coding one, but we can use the
# data to illustrate the trade-off concretely: how much does a stricter
# completeness/contamination filter cost us in terms of retained contigs?

quality_filters <- list(
  "no filter"             = function(df) df,
  "lenient (comp>50, contam<10)" = function(df) filter(df, completeness > 50, contamination < 10),
  "moderate (comp>70, contam<5)" = function(df) filter(df, completeness > 70, contamination < 5),
  "strict (comp>90, contam<1)"   = function(df) filter(df, completeness > 90, contamination < 1)
)

quality_filter_summary <- map_dfr(names(quality_filters), function(filter_name) {
  filtered <- quality_filters[[filter_name]](phage_inferences)
  tibble(
    filter = filter_name,
    n_retained = nrow(filtered),
    pct_retained = round(100 * nrow(filtered) / nrow(phage_inferences), 1)
  )
})

quality_filter_summary

# Discussion notes for students:
# - Being MORE strict reduces false positives (contaminated or incomplete
#   sequences that aren't really "good" viral genomes) but increases false
#   negatives (you throw away real, low-coverage, or partially-assembled
#   viruses, shrinking your dataset and potentially biasing it toward only
#   the most abundant/easy-to-assemble taxa).
# - Being LESS strict keeps more sequences (useful for diversity surveys,
#   exploratory analyses, or rare-taxa discovery) but increases the risk of
#   including artifacts, chimeras, or host contamination in downstream
#   analyses like phylogenetics or functional annotation, where a single bad
#   sequence can distort results.
# - In general: be MORE strict when the downstream analysis is sensitive to
#   individual sequence quality (e.g. genome annotation, phylogenomics,
#   building reference databases). Be LESS strict when you care more about
#   overall community-level patterns (e.g. diversity, presence/absence
#   across samples) where a few noisy sequences matter less in aggregate.


# ==========================================================================
# Q3. What is Read_GC about? What is GC? How variable is it? Why?
# ==========================================================================

# Read_GC is the GC content (the percentage of a sequence's bases that are
# guanine or cytosine, as opposed to adenine or thymine) calculated from the
# reads mapped back to each contig. GC content is a basic compositional
# property of DNA and varies by organism/lineage, so it can act as a rough
# fingerprint of what's being sequenced.

summary(phage_inferences$Read_GC)
sd(phage_inferences$Read_GC, na.rm = TRUE)

ggplot(phage_inferences, aes(x = Read_GC)) +
  geom_histogram() +
  theme_classic() +
  ggtitle("Distribution of Read_GC across phage contigs")

# Discussion notes for students:
# - GC content varies because different viral lineages (and their bacterial
#   hosts, in the case of proviruses) have evolved different baseline GC
#   contents — this is a known, heritable genomic trait, not noise.
# - A wide spread in Read_GC across contigs suggests you're sampling a
#   taxonomically diverse set of phages/hosts in this wastewater community,
#   rather than one dominant lineage.
# - It's also worth checking whether Read_GC clusters by type (virus vs
#   provirus), since proviruses inherit some compositional signal from their
#   host genome:
ggplot(phage_inferences, aes(x = type, y = Read_GC)) +
  geom_boxplot(aes(color = type)) +
  theme_classic()


# ==========================================================================
# Q4. Avg_fold vs contig_length, and vs viral_genes
# ==========================================================================

ggplot(phage_inferences, aes(x = contig_length, y = Avg_fold)) +
  geom_point(alpha = 0.4) +
  theme_classic() +
  ggtitle("Coverage (Avg_fold) vs contig length")

cor.test(phage_inferences$contig_length, phage_inferences$Avg_fold, method = "spearman")

ggplot(phage_inferences, aes(x = viral_genes, y = Avg_fold)) +
  geom_point(alpha = 0.4) +
  theme_classic() +
  ggtitle("Coverage (Avg_fold) vs viral gene count")

cor.test(phage_inferences$viral_genes, phage_inferences$Avg_fold, method = "spearman")

# Discussion notes for students:
# - We use Spearman correlation rather than Pearson because Avg_fold is
#   often heavily right-skewed (a few very highly-covered contigs), and
#   Spearman is robust to that skew since it works on ranks rather than raw
#   values.
# - In general, contig length and coverage are NOT expected to correlate
#   strongly: coverage reflects how ABUNDANT an organism is in the sample,
#   while length reflects how well that organism's genome assembled, which
#   depends more on assembly quality, repeat content, and strain
#   heterogeneity than on abundance.
# - Viral gene count may show a weak positive relationship with coverage,
#   since well-covered (likely longer, less fragmented) contigs naturally
#   have more room to contain identifiable genes — but this is partly
#   confounded with contig length itself.


# ==========================================================================
# Q5. Disagreement between completeness, contamination, checkv_quality
# ==========================================================================

# Group A: "complete" but highly contaminated
complete_contaminated <- phage_inferences %>%
  filter(completeness > 90, contamination > 10)

# Group B: incomplete but very clean
incomplete_clean <- phage_inferences %>%
  filter(completeness <= 90, contamination < 1)

nrow(complete_contaminated)
nrow(incomplete_clean)

# Combine for a side-by-side comparison
comparison_groups <- bind_rows(
  complete_contaminated %>% mutate(group = "Complete + contaminated"),
  incomplete_clean %>% mutate(group = "Incomplete + clean")
)

# Contig length comparison
ggplot(comparison_groups, aes(x = group, y = contig_length)) +
  geom_boxplot(aes(color = group)) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Contig length: complete+contaminated vs incomplete+clean")

# Viral gene count comparison
ggplot(comparison_groups, aes(x = group, y = viral_genes)) +
  geom_boxplot(aes(color = group)) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Viral gene count: complete+contaminated vs incomplete+clean")

# Quick numeric summary alongside the plots
comparison_groups %>%
  group_by(group) %>%
  summarise(
    n = n(),
    median_length = median(contig_length),
    median_viral_genes = median(viral_genes)
  )

# Discussion notes for students:
# - If "complete + contaminated" contigs tend to be LONGER with similar or
#   higher viral gene counts, that's consistent with these being real viral
#   genomes that happen to carry some extra (possibly host-derived or
#   misassembled) sequence — i.e. contamination doesn't necessarily mean
#   "not a virus," just "not a pure single-genome assembly."
# - If "incomplete + clean" contigs are SHORTER with fewer viral genes,
#   that's consistent with these being genuine fragments of real viruses
#   that just didn't assemble fully, rather than artifacts.
# - The broader point: completeness alone is NOT a reliable single filter,
#   since a contig can be incomplete (because of assembly issues) yet
#   perfectly clean and trustworthy as a partial genome, while a "complete"
#   contig can still be impure. This is the rationale for using checkv
#   quality tiers (which combine multiple metrics) rather than any single
#   metric in isolation.


# ==========================================================================
# Q6. viral_gene_ratio = viral_genes / gene_count
# ==========================================================================

phage_inferences <- phage_inferences %>%
  mutate(viral_gene_ratio = viral_genes / gene_count)

# Plot for viruses vs proviruses
ggplot(phage_inferences, aes(x = type, y = viral_gene_ratio)) +
  geom_boxplot(aes(color = type)) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Viral gene ratio: virus vs provirus")

# Discussion notes for students:
# - Expect proviruses to generally show LOWER viral_gene_ratio than free
#   viruses, since proviruses are embedded in host DNA and the gene-calling
#   window often includes some flanking host genes alongside viral ones.

# Filter: at least half viral genes AND length > 5000bp
phage_inferences_ratio_filtered <- phage_inferences %>%
  filter(viral_gene_ratio > 0.5, contig_length > 5000)

nrow(phage_inferences_ratio_filtered)
nrow(phage_inferences_filtered)  # the simple length-only filter from earlier

# How many were lost by adding the ratio criterion?
n_lost <- nrow(phage_inferences_filtered) - nrow(phage_inferences_ratio_filtered)
pct_lost <- round(100 * n_lost / nrow(phage_inferences_filtered), 1)

n_lost
pct_lost

# Which specific contigs were lost?
lost_contigs <- phage_inferences_filtered %>%
  anti_join(phage_inferences_ratio_filtered, by = "contig_id")

summary(lost_contigs$viral_gene_ratio)
summary(lost_contigs$contig_length)

# Discussion notes for students:
# - Whether it was "worth it" depends on what you saw above: if the lost
#   contigs have a viral_gene_ratio just barely under 0.5, you're discarding
#   borderline cases that are arguably still fine. If they cluster at very
#   low ratios (e.g. < 0.2), the ratio filter is doing real work in removing
#   sequences that are mostly non-viral despite being long enough to pass
#   the simple length filter — i.e. it's catching hybrid/misclassified
#   contigs that length alone would have let through.


# ==========================================================================
# Q7. Coverage bins on the filtered phage contigs (NEEDS CORRECTION)
# ==========================================================================

phage_inferences_filtered_binned <- phage_inferences_filtered %>%
  left_join(phage_inferences %>% select(contig_id, Avg_fold), by = "contig_id") %>%
  mutate(coverage_bin = cut(
    Avg_fold,
    breaks = c(0, 10, 50, 100, Inf),
    labels = c("0-10x", "10-50x", "50-100x", ">100x"),
    right = FALSE
  ))

coverage_bin_counts <- phage_inferences_filtered_binned %>%
  count(coverage_bin)

coverage_bin_counts

ggplot(coverage_bin_counts, aes(x = coverage_bin, y = n)) +
  geom_col() +
  xlab("Coverage bin") +
  ylab("Number of contigs") +
  theme_classic() +
  ggtitle("Distribution of filtered phage contigs across coverage bins")

# Discussion notes for students:
# - A distribution skewed heavily toward the lowest bin (0-10x), with
#   progressively fewer contigs in each higher bin, suggests the viral
#   community is dominated by many low-abundance lineages with only a
#   handful of highly abundant ones — a common pattern in diverse natural
#   communities (sometimes summarized as "few abundant, many rare").
# - A more uniform spread across bins would instead suggest a community
#   with several lineages at comparable abundance, without one or two
#   dominant taxa.
# - This distribution shape connects directly back to the alpha diversity
#   work from Exercise 2: a skewed coverage distribution is consistent with
#   the kind of unevenness that Shannon diversity (as opposed to simple
#   richness) is specifically designed to capture.
