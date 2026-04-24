# Day 3 Activity. Microbial taxonomic and functional ecology analysis.

# Load libraries
library(tidyverse)
library(vegan)
library(car)
library(lme4)
library(ggpubr)


#------------------------ Load and inspect data ----------------------------------------------
setwd("~/Documents/GitHub/ucd_bioinformaticscourse_2026/")

# ---- metadata
md <- read_csv("./data/metadata.csv")
# take a quick look at it
glimpse(md)
# get some summary numbers
summary(md)
# check for missing values per column
md %>% summarise(across(everything(), ~ sum(is.na(.))))


# ---- taxonomic count table
taxonomy_counts <- read_csv("./data/taxonomy_counts_table.csv")
glimpse(taxonomy_counts)
summary(taxonomy_counts)
taxonomy_counts %>% summarise(across(everything(), ~ sum(is.na(.))))

# we need to substitute the NAs with 0s
taxonomy_counts <- taxonomy_counts %>% mutate(across(where(is.numeric), ~ replace_na(., 0)))

# ---- gene (KO) count table
ko_counts <- read_table("./data/ko_table.txt")
glimpse(ko_counts)
summary(ko_counts)
ko_counts %>% summarise(across(everything(), ~ sum(is.na(.))))



#------------------------ Reorder data so they all agree ----------------------------------------------
# this is important later on for the alpha diversity calculation
md <- md %>%
  arrange(sample)

taxonomy_counts <- taxonomy_counts %>%
  column_to_rownames(var = "taxon") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  arrange(sample)

ko_counts <- ko_counts %>%
  column_to_rownames(var = "KO") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  arrange(sample)

# check that they are all in the same order
all(md$sample == taxonomy_counts$sample)
all(md$sample == ko_counts$sample)



#------------------------ Normalization ----------------------------------------------

# First we extract the total_reads information from the md
total_reads_table <- md %>%
  select(sample, total_reads)

# Then we include this information with the counts for each taxon
taxonomy_counts_plustotalreads <- taxonomy_counts %>%
  pivot_longer(-sample, names_to = "taxon", values_to = "count") %>%
  left_join(., total_reads_table , by = "sample")

# Finally we use this information to do a RPM normalization (Reads per million)
rpm_taxonomy_counts <- taxonomy_counts_plustotalreads %>%
  mutate(RPM = (count / total_reads) * 1e6) %>%
  select(-count, -total_reads) %>%
  pivot_wider(names_from = "taxon", values_from = "RPM")

# Now we do the same for the ko table
ko_counts_plustotalreads <- ko_counts %>%
  pivot_longer(-sample, names_to = "KO", values_to = "count") %>%
  left_join(., total_reads_table , by = "sample")

# Finally we use this information to do a RPM normalization (Reads per million)
rpm_ko_counts <- ko_counts_plustotalreads %>%
  mutate(RPM = (count / total_reads) * 1e6) %>%
  select(-count, -total_reads) %>%
  pivot_wider(names_from = "KO", values_from = "RPM")



#------------------------ Alpha diveristy ----------------------------------------------

# First we need to prepare the table
taxonomy_counts_forrich <- taxonomy_counts %>%
  column_to_rownames(var = "sample")

# Then we calculate the rarefaction number
summary(rowSums(taxonomy_counts_forrich))
# We can also visualize it
hist(rowSums(taxonomy_counts_forrich))

# calculate diversity measures
md$taxonomic_rich <- specnumber(rrarefy(taxonomy_counts_forrich, sample = 7734272))
md$taxonomic_shannon <- diversity(rrarefy(taxonomy_counts_forrich, sample = 7734272), index = "shannon")


# Now we do the same for the kos
ko_counts_forrich <- ko_counts %>%
  column_to_rownames(var = "sample")
summary(rowSums(ko_counts_forrich))
hist(rowSums(ko_counts_forrich))

md$functional_rich <- specnumber(rrarefy(ko_counts_forrich, sample = 2352294))
md$functional_shannon <- diversity(rrarefy(ko_counts_forrich, sample = 2352294), index = "shannon")



# ----------- Alpha diversity plots

# First we need to convert the explanatory variables to factors: 
md$urban <- as.factor(md$urban)
md$site <- as.factor(md$site)

# Taxonomic richness
ggplot(md, aes(x = urban, y = taxonomic_rich))+
  geom_jitter(aes(color = urban))+
  geom_boxplot(aes(color = urban), outliers = F) +
  xlab("")+
  ylab("Number of bacterial taxa")+
  theme_bw() +
  theme(legend.position = "none")
ggsave("./figures/taxonomic_rich.pdf", width = 5, height = 7, units = "in")

# Taxonomic shannon
ggplot(md, aes(x = urban, y = taxonomic_shannon))+
  geom_jitter(aes(color = urban))+
  geom_boxplot(aes(color = urban), outliers = F) +
  xlab("")+
  ylab("Taxonomic Shannon Diversity")+
  theme_bw() +
  theme(legend.position = "none")
ggsave("./figures/taxonomic_shannon.pdf", width = 5, height = 7, units = "in")


# KO richness
ggplot(md, aes(x = urban, y = functional_rich))+
  geom_jitter(aes(color = urban))+
  geom_boxplot(aes(color = urban), outliers = F) +
  xlab("")+
  ylab("Number of KOs")+
  theme_bw() +
  theme(legend.position = "none")
ggsave("./figures/functional_rich.pdf", width = 5, height = 7, units = "in")


# KO richness
ggplot(md, aes(x = urban, y = functional_shannon))+
  geom_jitter(aes(color = urban))+
  geom_boxplot(aes(color = urban), outliers = F) +
  xlab("")+
  ylab("Functional Shannon Diversity")+
  theme_bw() +
  theme(legend.position = "none")
ggsave("./figures/functional_shannon.pdf", width = 5, height = 7, units = "in")




# ---------- Statistics

# simple linear model
anova(lm(taxonomic_rich ~ urban, data = md))
anova(lm(taxonomic_shannon ~ urban, data = md))
anova(lm(functional_rich ~ urban, data = md))
anova(lm(functional_shannon ~ urban, data = md))

# mixed effects model
Anova(lmer(taxonomic_rich ~ urban + (1|site), data = md))
Anova(lmer(taxonomic_shannon ~ urban + (1|site), data = md))
Anova(lmer(functional_rich ~ urban + (1|site), data = md))
Anova(lmer(functional_shannon ~ urban + (1|site), data = md))



# ------------ Are they correlated?
ggplot(md, aes(x = taxonomic_rich, y = functional_rich))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor()+
  xlab("Taxonomic Richness")+
  ylab("Functional Richness")+
  theme_bw()

# by urban vs. natural
ggplot(md, aes(x = taxonomic_rich, y = functional_rich, color = urban, group = urban))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor()+
  xlab("Taxonomic Richness")+
  ylab("Functional Richness")+
  theme_bw()

ggplot(md, aes(x = taxonomic_shannon, y = functional_shannon))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor()+
  xlab("Taxonomic Shannon Diversity")+
  ylab("Functional Shannon Diversity")+
  theme_bw()

ggplot(md, aes(x = taxonomic_shannon, y = functional_shannon, color = urban, group = urban))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor()+
  xlab("Taxonomic Shannon Diversity")+
  ylab("Functional Shannon Diversity")+
  theme_bw()





#------------------------ Beta diveristy ----------------------------------------------
# First we prepare the table
taxonomy_counts_forbray <- rpm_taxonomy_counts %>%
  column_to_rownames(var = "sample")

# Then we calculate a dissimilarity matrix
taxonomy.bray<- vegdist(taxonomy_counts_forbray, method="bray")

# Then we calculate the ordination (non-multidimensional scaling) and add it to md for plotting
taxonomy.nmds <- metaMDS(taxonomy.bray)
md$taxonomy.Axis01 = taxonomy.nmds$points[,1]
md$taxonomy.Axis02 = taxonomy.nmds$points[,2]
taxonomy.nmds$stress #0.07968216

# Now we do the same with the ko table
# First we prepare the table
ko_counts_forbray <- rpm_ko_counts %>%
  column_to_rownames(var = "sample")

ko.bray<- vegdist(ko_counts_forbray, method="bray")

ko.nmds <- metaMDS(ko.bray)
md$ko.Axis01 = ko.nmds$points[,1]
md$ko.Axis02 = ko.nmds$points[,2]
ko.nmds$stress #0.02760141



# ------------- Plot
# Ordination plot
ggplot(md, aes(taxonomy.Axis01, taxonomy.Axis02))+
  geom_point(aes(color = urban), size=3.5)+
  stat_ellipse(aes(color= urban)) +
  #scale_color_see()+
  theme_classic()+
  #scale_color_manual(values = md$coloring)+
  theme(legend.position="right", text = element_text(size=12))

ggplot(md, aes(taxonomy.Axis01, taxonomy.Axis02))+
  geom_point(aes(color = site), size=3.5)+
  stat_ellipse(aes(color= site)) +
  #scale_color_see()+
  theme_classic()+
  #scale_color_manual(values = md$coloring)+
  theme(legend.position="right", text = element_text(size=12))

ggplot(md, aes(ko.Axis01, ko.Axis02))+
  geom_point(aes(color = urban), size=3.5)+
  stat_ellipse(aes(color= urban)) +
  #scale_color_see()+
  theme_classic()+
  #scale_color_manual(values = md$coloring)+
  theme(legend.position="right", text = element_text(size=12))

ggplot(md, aes(ko.Axis01, ko.Axis02))+
  geom_point(aes(color = site), size=3.5)+
  stat_ellipse(aes(color= site)) +
  #scale_color_see()+
  theme_classic()+
  #scale_color_manual(values = md$coloring)+
  theme(legend.position="right", text = element_text(size=12))



# ------------- Stats
adonis2(taxonomy.bray ~  urban, data = md, permutations = 999, method = "bray")
adonis2(ko.bray ~  urban, data = md, permutations = 999, method = "bray")
adonis2(taxonomy.bray ~  site, data = md, permutations = 999, method = "bray")
adonis2(ko.bray ~  site, data = md, permutations = 999, method = "bray")
