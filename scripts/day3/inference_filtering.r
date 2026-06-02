# Filter phage and plasmid output from geNomad
# Day 3 Exercise 1, UCD 
# Author: Mery Touceda-Suarez
# Date: June 2026



# Load libraries
library(tidyverse)



# ---------------------------------------- viruses  ----------------------------------------

# --------------- 1. Load data
phage_inferences <- read_csv("~/Documents/GitHub/ucd_bioinformaticscourse_2026/data/phages_summary.csv") 


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
