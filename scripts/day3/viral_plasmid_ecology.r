# Phage and plasmid ecology analyses dww
# Day 3 Exercise 2, UCD 
# Author: Mery Touceda-Suarez
# Date: June 2026

# Load libraries
library(tidyverse)
library(vegan)



# ---------------------------------------- Load and prep data ----------------------------------------

# ---------------- Metadata
md <- read.csv("~/Documents/GitHub/ucd_bioinformaticscourse_2026/data/metadata_day3.csv")
md <- md %>%
  arrange(desc(sample)) %>%
  # make site variable
  mutate(site = paste(country, city, plant)) %>%
  # make date character into date class
  mutate(collection_date_n = as.Date(collection_date))


# ---------------- Phages
votus_count <- read.table("~/Documents/GitHub/ucd_bioinformaticscourse_2026/data/votus_count_table.txt", sep = "\t", header = T)
# prepare data
votus_clean <- votus_count %>%
  pivot_longer(!Contig, names_to = "sample", values_to = "counts") %>%
  mutate(sample = str_remove(sample, ".Read.Count")) %>%
  filter(sample %in% md$sample) %>%
  pivot_wider(names_from = "Contig", values_from = "counts") %>%
  arrange(desc(sample)) 

  
# ---------------- Plasmids
potus_count <- read.table("~/Documents/GitHub/ucd_bioinformaticscourse_2026/data/potus_count_table.txt", sep = "\t", header = T)
# prepare data
potus_clean <- potus_count %>%
  pivot_longer(!Contig, names_to = "sample", values_to = "counts") %>%
  mutate(sample = str_remove(sample, ".Read.Count")) %>%
  filter(sample %in% md$sample) %>%
  pivot_wider(names_from = "Contig", values_from = "counts") %>%
  arrange(desc(sample)) 


# confirm that samples are in same order
md$sample == votus_clean$sample
md$sample == potus_clean$sample





# ---------------------------------------- Alppha diversity ----------------------------------------
# prepare count table
votus_for_rich <- votus_clean %>%
  column_to_rownames(var = "sample") %>%
  mutate(across(where(is.character), as.numeric))

potus_for_rich <- potus_clean %>%
  column_to_rownames(var = "sample") %>%
  mutate(across(where(is.character), as.numeric))

# check counts distribution
summary(rowSums(votus_for_rich))
hist(rowSums(votus_for_rich))

summary(rowSums(potus_for_rich))
hist(rowSums(potus_for_rich))

# calculate diversity measures
md$viral_rich <- specnumber(rrarefy(votus_for_rich, sample = 71456))
md$viral_shannon <- diversity(rrarefy(votus_for_rich, sample = 71456), index = "shannon")
md$plasmid_rich <- specnumber(rrarefy(potus_for_rich, sample = 63508))
md$plasmid_shannon <- diversity(rrarefy(potus_for_rich, sample = 63508), index = "shannon")

# differences per site
ggplot(md, aes(x = site, y = viral_shannon))+
  geom_boxplot(aes(color = site)) +
  xlab("")+
  ylab("Number of vOTUs")+
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("CHANGE THIS PATH! /figures/viral_shannon_site.pdf", width = 5, height = 7, units = "in")


# differences per site
ggplot(md, aes(x = site, y = plasmid_shannon))+
  geom_boxplot(aes(color = site)) +
  xlab("")+
  ylab("Number of potus")+
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))





# ---------------------------------------- Community Composition ----------------------------------------


# --------------------------------- RPKM abundance normalization

# PHAGES
total_reads <- md %>%
  select(sample, total_reads) 

# load votus length
votus_length <- read_csv("~/Documents/GitHub/ucd_bioinformaticscourse_2026/data/phages_summary.csv") %>%
  select(c("votu", "contig_length"))

# convert count table to long format
votus_count_long <- votus_clean %>%
  pivot_longer(cols = -c("sample"), values_to = "count", names_to= "votu") 

votus_count_long_plusreads <- left_join(votus_count_long , total_reads, by = "sample")
votus_count_long_plusreads_pluslength <- left_join(votus_count_long_plusreads , votus_length, by = "votu")

# compute normalization
# RPKM (for viral load)
rpkm_votu_count <- votus_count_long_plusreads_pluslength %>%
  mutate(RPKM = (count * 1e6) / (contig_length * total_reads)) %>%
  dplyr::select(-c("count", "total_reads", "contig_length")) %>%
  pivot_wider(names_from = "sample", values_from = "RPKM") %>%
  column_to_rownames(var = "votu") %>%
  drop_na()

# TPM (for composition comparisons, beta)
tpm_votu_count <- votus_count_long_plusreads_pluslength %>%
  mutate(RPK = (count * 1e6) / contig_length) %>%
  group_by(sample) %>%
  mutate(TPM = RPK / sum(RPK, na.rm = TRUE) * 1e6) %>%
  ungroup() %>%
  dplyr::select(-c("count", "total_reads", "contig_length", "RPK")) %>%
  pivot_wider(names_from = "sample", values_from = "TPM") %>%
  column_to_rownames(var = "votu") %>%
  drop_na()




# PLASMIDS
potus_length <- read_csv("~/Documents/GitHub/ucd_bioinformaticscourse_2026/data/plasmids_summary.csv") %>%   # change for your selection
  select(c("seq_name", "length")) %>%
  rename(potu = seq_name)

potus_count_long <- potus_clean %>%
  pivot_longer(cols = -c("sample"), values_to = "count", names_to= "potu") 

potus_count_long_plusreads <- left_join(potus_count_long , total_reads, by = "sample")
potus_count_long_plusreads_pluslength <- left_join(potus_count_long_plusreads , potus_length, by = "potu")

# compute normalization
rpkm_potu_count <- potus_count_long_plusreads_pluslength %>%
  mutate(RPKM = (count * 1e9) / (length * total_reads)) %>%
  dplyr::select(-c("count", "total_reads", "length")) %>%
  pivot_wider(names_from = "sample", values_from = "RPKM") %>%
  column_to_rownames(var = "potu") %>%
  drop_na()





# --------------------------------- Community composition

# PHAGES
viral.bray <- vegdist(t(rpkm_votu_count), method="bray")

#ordination (non-multidimensional scaling)
viral.nmds <- metaMDS(viral.bray)
md$bracken.Axis01 = viral.nmds$points[,1]
md$bracken.Axis02 = viral.nmds$points[,2]
viral.nmds$stress # 0.1588582

# stats
adonis2(viral.bray ~  site, data = md, permutations = 999, method = "bray")
adonis2(viral.bray ~  collection_date_n, data = md, permutations = 999, method = "bray")


# Ordination plot
ggplot(md, aes(bracken.Axis01, bracken.Axis02))+
  geom_point(aes(alpha=collection_date_n, color = site), size=3.5)+
  stat_ellipse(aes(color=site)) +
  theme_classic()+
  theme(legend.position="right", text = element_text(size=12))
#ggsave("CHANGE THIS PATH! /viral_composition_bydate.pdf", device = "pdf", width = 9, height = 5 , units = "in")


# PLASMIDS
plasmid.bray <- vegdist(t(rpkm_potu_count), method="bray")

#ordination (non-multidimensional scaling)
plasmid.nmds <- metaMDS(plasmid.bray)
md$bracken.Axis01 = plasmid.nmds$points[,1]
md$bracken.Axis02 = plasmid.nmds$points[,2]
plasmid.nmds$stress # 0.1588582

# stats
adonis2(plasmid.bray ~  site, data = md, permutations = 999, method = "bray")
adonis2(plasmid.bray ~  collection_date_n, data = md, permutations = 999, method = "bray")

# Ordination plot
ggplot(md, aes(bracken.Axis01, bracken.Axis02))+
  geom_point(aes(alpha=collection_date_n, color = site), size=3.5)+
  stat_ellipse(aes(color=site)) +
  theme_classic()+
  theme(legend.position="right", text = element_text(size=12))
#ggsave("CHANGE THIS PATH!/plasmid_composition_bydate.pdf", device = "pdf", width = 9, height = 5 , units = "in")



# Other questions
# What variables are determining viral and phage diversity and community composition, apart from the location and the time?
# Do plasmid and viral diversity correlate? 
# Do their community compositions vary together? 
# What are the most abundant phages/plasmids? What is their length distribution? 
# Do we see changes in phage/plasmid abundance associated with their length? (i.e. shorter viruses are more abundant)







