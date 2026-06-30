# Filter phage and plasmid output from geNomad
# Day 3 Exercise 1, UCD 
# Author: Mery Touceda-Suarez
# Date: June 2026



# Load libraries
library(tidyverse)



# ---------------------------------------- viruses  ----------------------------------------

# --------------- 1. Load data
phage_inferences <- read_csv("~/Documents/GitHub/ucd_bioinformaticscourse_2026/data/phages_summary.csv") 

colnames(phage_inferences)

# --------------- 2. Inspect data
summary(phage_inferences)


# --------------- 3. Explore
# We use this exploration to choose filter parameters
# We have viruses and proviruses, we should explore those separately maybe

phage_inferences <- phage_inferences %>%
  mutate(type = case_when(str_detect(votu, "provirus") ~ "provirus", 
                   TRUE ~ "virus"))

# How many of each type?
dim(phage_inferences %>% filter(type == "provirus"))
dim(phage_inferences %>% filter(type == "virus"))

# How do we know which is provirus and which one is virus (host genes)
summary(phage_inferences %>% filter(type == "provirus") %>% select(host_genes))
summary(phage_inferences %>% filter(type == "virus") %>% select(host_genes))


# 3.a. contig length
# all together
hist(phage_inferences$contig_length)
summary(phage_inferences$contig_length)

# separated by type
ggplot(phage_inferences, aes(y = contig_length))+
  geom_histogram()+
  facet_wrap(~type)

summary(phage_inferences %>% filter(type == "provirus") %>% select(contig_length))
summary(phage_inferences %>% filter(type == "virus") %>% select(contig_length))


# 3.b contamination 

# 3.c completeness




# --------------- 4. Filter
# I will filter out viruses smaller than 5000bp

phage_inferences_filtered <- phage_inferences %>%
  filter(contig_length > 5000)
# let's see how many we have left
dim(phage_inferences_filtered )

# We can be more strict, I want them also to have more than two viral genes
phage_inferences_filtered_strict <- phage_inferences %>%
  filter(contig_length > 5000, viral_genes > 2)

# Or less strict, I want them to be >5000bp OR to have more than 2 viral genes
phage_inferences_filtered_strict <- phage_inferences %>%
  filter(contig_length > 5000 | viral_genes > 2)




# ---------------------------------------- plasmids  ----------------------------------------
# Let's now look at the plasmids!

# --------------- 1. Load data
plasmid_inferences <- read_csv("~/Documents/GitHub/ucd_bioinformaticscourse_2026/data/plasmids_summary.csv") 


# --------------- 2. Inspect data
summary(plasmid_inferences)


# --------------- 3. Explore
# We use this exploration to choose filter parameters

# 3.a. length
hist(plasmid_inferences$length)
summary(plasmid_inferences$length)

# 3.b. plasmid_score (geNomad's confidence score for plasmid classification)
hist(plasmid_inferences$plasmid_score)
summary(plasmid_inferences$plasmid_score)

# 3.c. topology, if present (circular vs linear)
# How many of each topology?
plasmid_inferences %>% count(topology)

# length distribution separated by topology
ggplot(plasmid_inferences, aes(x = length)) +
  geom_histogram() +
  facet_wrap(~topology)

# 3.d. number of genes per plasmid
hist(plasmid_inferences$n_genes)
summary(plasmid_inferences$n_genes)

# relationship between score and length
ggplot(plasmid_inferences, aes(x = length, y = plasmid_score)) +
  geom_point(alpha = 0.4) +
  theme_classic()

# --------------- 4. Filter
# I will filter out plasmids with a low plasmid_score
plasmid_inferences_filtered <- plasmid_inferences %>%
  filter(plasmid_score > 0.8)

# let's see how many we have left
dim(plasmid_inferences_filtered)

# We can be more strict: score AND a minimum length
plasmid_inferences_filtered_strict <- plasmid_inferences %>%
  filter(plasmid_score > 0.8, length > 1000)

# Or less strict: score OR a minimum length
plasmid_inferences_filtered_lenient <- plasmid_inferences %>%
  filter(plasmid_score > 0.8 | length > 1000)


































# Sample code: comparing the effect of filtering permissiveness
# For both phages (viruses) and plasmids


# ==========================================================================
# PHAGES (viruses)
# ==========================================================================
# --------------- A. Compare a few named filters side by side
# Build a small function that takes a filter expression (as a function) and reports counts

phage_filters <- list(
  "no filter"          = function(df) df,
  "length > 5000"      = function(df) filter(df, contig_length > 5000),
  "strict (AND)"       = function(df) filter(df, contig_length > 5000, viral_genes > 2),
  "lenient (OR)"       = function(df) filter(df, contig_length > 5000 | viral_genes > 2)
)

phage_filter_summary <- map_dfr(names(phage_filters), function(filter_name) {
  filtered <- phage_filters[[filter_name]](phage_inferences)
  tibble(
    filter = filter_name,
    n_retained = nrow(filtered),
    pct_retained = round(100 * nrow(filtered) / nrow(phage_inferences), 1)
  )
})

phage_filter_summary

# --------------- B. Threshold sweep: how does contig_length cutoff affect retained count?

length_thresholds <- seq(1000, 10000, by = 500)

phage_length_sweep <- map_dfr(length_thresholds, function(thresh) {
  tibble(
    threshold = thresh,
    n_retained = sum(phage_inferences$contig_length > thresh)
  )
})

ggplot(phage_length_sweep, aes(x = threshold, y = n_retained)) +
  geom_line() +
  geom_point() +
  xlab("Minimum contig length (bp)") +
  ylab("Number of viral contigs retained") +
  ggtitle("Effect of contig length cutoff on retained phages") +
  theme_classic()

# --------------- C. Threshold sweep: viral_genes cutoff

gene_thresholds <- seq(0, 10, by = 1)

phage_gene_sweep <- map_dfr(gene_thresholds, function(thresh) {
  tibble(
    threshold = thresh,
    n_retained = sum(phage_inferences$viral_genes > thresh)
  )
})

ggplot(phage_gene_sweep, aes(x = threshold, y = n_retained)) +
  geom_line() +
  geom_point() +
  xlab("Minimum number of viral genes") +
  ylab("Number of viral contigs retained") +
  ggtitle("Effect of viral gene count cutoff on retained phages") +
  theme_classic()

# --------------- D. Combined sweep: AND vs OR across length thresholds (genes fixed at 2)

phage_combined_sweep <- map_dfr(length_thresholds, function(thresh) {
  tibble(
    threshold = thresh,
    n_AND = sum(phage_inferences$contig_length > thresh & phage_inferences$viral_genes > 2),
    n_OR  = sum(phage_inferences$contig_length > thresh | phage_inferences$viral_genes > 2)
  )
}) %>%
  pivot_longer(cols = c(n_AND, n_OR), names_to = "logic", values_to = "n_retained")

ggplot(phage_combined_sweep, aes(x = threshold, y = n_retained, color = logic)) +
  geom_line() +
  geom_point() +
  xlab("Minimum contig length (bp)") +
  ylab("Number of viral contigs retained") +
  ggtitle("AND vs OR filtering across length thresholds") +
  theme_classic()


# ==========================================================================
# PLASMIDS
# ==========================================================================
# --------------- A. Compare a few named filters side by side
colnames(plasmid_inferences)
plasmid_filters <- list(
  "no filter"             = function(df) df,
  "score > 0.8"           = function(df) filter(df, plasmid_score > 0.8),
  "strict (AND)"          = function(df) filter(df, plasmid_score > 0.8, length > 1000),
  "lenient (OR)"          = function(df) filter(df, plasmid_score > 0.8 | length > 1000)
)

plasmid_filter_summary <- map_dfr(names(plasmid_filters), function(filter_name) {
  filtered <- plasmid_filters[[filter_name]](plasmid_inferences)
  tibble(
    filter = filter_name,
    n_retained = nrow(filtered),
    pct_retained = round(100 * nrow(filtered) / nrow(plasmid_inferences), 1)
  )
})

plasmid_filter_summary

# --------------- B. Threshold sweep: plasmid_score cutoff

score_thresholds <- seq(0.5, 0.95, by = 0.05)

plasmid_score_sweep <- map_dfr(score_thresholds, function(thresh) {
  tibble(
    threshold = thresh,
    n_retained = sum(plasmid_inferences$plasmid_score > thresh, na.rm = TRUE)
  )
})

ggplot(plasmid_score_sweep, aes(x = threshold, y = n_retained)) +
  geom_line() +
  geom_point() +
  xlab("Minimum plasmid score") +
  ylab("Number of plasmids retained") +
  ggtitle("Effect of plasmid score cutoff on retained plasmids") +
  theme_classic()

# --------------- C. Threshold sweep: length cutoff

plasmid_length_thresholds <- seq(500, 10000, by = 500)

plasmid_length_sweep <- map_dfr(plasmid_length_thresholds, function(thresh) {
  tibble(
    threshold = thresh,
    n_retained = sum(plasmid_inferences$length > thresh, na.rm = TRUE)
  )
})

ggplot(plasmid_length_sweep, aes(x = threshold, y = n_retained)) +
  geom_line() +
  geom_point() +
  xlab("Minimum plasmid length (bp)") +
  ylab("Number of plasmids retained") +
  ggtitle("Effect of length cutoff on retained plasmids") +
  theme_classic()

# --------------- D. Combined sweep: AND vs OR across score thresholds (length fixed at 1000)

plasmid_combined_sweep <- map_dfr(score_thresholds, function(thresh) {
  tibble(
    threshold = thresh,
    n_AND = sum(plasmid_inferences$plasmid_score > thresh & plasmid_inferences$length > 1000, na.rm = TRUE),
    n_OR  = sum(plasmid_inferences$plasmid_score > thresh | plasmid_inferences$length > 1000, na.rm = TRUE)
  )
}) %>%
  pivot_longer(cols = c(n_AND, n_OR), names_to = "logic", values_to = "n_retained")

ggplot(plasmid_combined_sweep, aes(x = threshold, y = n_retained, color = logic)) +
  geom_line() +
  geom_point() +
  xlab("Minimum plasmid score") +
  ylab("Number of plasmids retained") +
  ggtitle("AND vs OR filtering across score thresholds") +
  theme_classic()
