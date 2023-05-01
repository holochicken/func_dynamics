# Temporal development of microbiome composition - updated 22/03/2023

# Packages
library(ape)
library(vegan)
library(nlme)
library(sjPlot)
library(hillR)
library(tidyverse)


# 1- Load data ---------------------
# MAG counts and metadata
mag_counts <-
  read_tsv("data/metagenomics/mag_counts.tsv") %>%
  column_to_rownames(var = 'mag_id')

metadata <-
  read_tsv("data/metadata/metadata.tsv") %>%
  mutate(sampling_time = factor(sampling_time, levels = c('7',
                                                          '21',
                                                          '35')))

metadata_2 <-
  metadata %>%
  dplyr::rename(site2 = 'animal_code',
                sampling_time_2 = 'sampling_time',
                pen_2 = 'pen',
                trial_2 = 'trial')

# MAG characteristics
stats <-
  read_tsv("data/metagenomics/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)

taxonomy <-
  read_tsv("data/metagenomics/taxonomy.tsv")

tree <-
  read.tree("data/metagenomics/tree.nwk")

kegg_paths <-
  read_tsv("data/metagenomics/dram.tsv") %>%
  column_to_rownames(var = 'mag_id') %>%
  select(1:399)


# 2- Standardisation ---------------------
# By MAG genome length
mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  data.frame()

# By sequencing depth
hel <- decostand(mag_weighted, 'hellinger')


rm(mag_counts, mag_weighted, stats)
# 3- Beta diversity ---------------------
# 3.1- Neutral
betadiv_n <- hill_taxa_parti_pairwise(hel, q = 1, pairs = "full")

# Calculate local dissimilarity
betadiv_n_dis <-
  betadiv_n %>%
  mutate(local_dissimilarity = 1 - local_similarity) %>%
  select(site1, site2, local_dissimilarity) %>%
  write_tsv(file = "results/tables/betadiv_n_dis.tsv")

# Compare individuals from same trial and pen between days 7-21 and 21-35
betadiv_n_dis <-
  read_tsv(file = "results/tables/betadiv_n_dis.tsv") %>%
  rename(animal_code = 'site1') %>%
  left_join(metadata, by = 'animal_code') %>%
  left_join(metadata_2, by = 'site2') %>%
  mutate(diff = case_when(sampling_time == '7' & sampling_time_2 == '21' ~ '7_21',
                          sampling_time == '21' & sampling_time_2 == '35' ~ '21_35',
                          sampling_time == '7' & sampling_time_2 == '35' ~ '7_35')) %>%
  mutate(diff_trial = case_when(trial == 'CA' & trial_2 == 'CA' ~ 'CA',
                                trial == 'CB' & trial_2 == 'CB' ~ 'CB')) %>%
  mutate(diff_pen = case_when(pen == pen_2  ~ 'same')) %>%
  drop_na(diff, diff_trial, diff_pen)

# Stats
betadiv_n_dis %>%
  group_by(diff) %>%
  summarise(mean = mean(local_dissimilarity), sd = sd(local_dissimilarity))

# Plot
pdf("results/figures/figS2a_beta_p_div.pdf", width = 5, height = 4)

plot_n <-
  betadiv_n_dis %>%
  filter(!diff == '7_35') %>%
  mutate(diff = factor(diff, levels = c('7_21', '21_35'))) %>%
  ggplot(aes(x = diff, y = local_dissimilarity, fill = diff)) +
  geom_boxplot() +
  ylim(0,0.5) +
  theme_minimal() +
  theme(legend.position = 'none') +
  ggtitle('Neutral beta div')

dev.off()


# 3.2- Phylogenetic
betadiv_p <- hill_phylo_parti_pairwise(hel, q = 1, tree, pairs = "full")

# Calculate local dissimilarity
betadiv_p_dis <-
  betadiv_p %>%
  mutate(local_dissimilarity = 1 - local_similarity) %>%
  select(site1, site2, local_dissimilarity) %>%
  write_tsv("results/tables/betadiv_p_dis.tsv")

# Compare individuals from same trial and pen between days 7-21 and 21-35
betadiv_p_dis <-
  read_tsv("results/tables/betadiv_p_dis.tsv") %>%
  rename(animal_code = 'site1') %>%
  left_join(metadata, by = 'animal_code') %>%
  left_join(metadata_2, by = 'site2') %>%
  mutate(diff = case_when(sampling_time == '7' & sampling_time_2 == '21' ~ '7_21',
                          sampling_time == '21' & sampling_time_2 == '35' ~ '21_35',
                          sampling_time == '7' & sampling_time_2 == '35' ~ '7_35')) %>%
  mutate(diff_trial = case_when(trial == 'CA' & trial_2 == 'CA' ~ 'CA',
                                trial == 'CB' & trial_2 == 'CB' ~ 'CB')) %>%
  mutate(diff_pen = case_when(pen == pen_2  ~ 'same')) %>%
  drop_na(diff, diff_trial, diff_pen)


# Stats
betadiv_p_dis %>%
  group_by(diff) %>%
  summarise(mean = mean(local_dissimilarity), sd = sd(local_dissimilarity))

# Plot
pdf("results/figures/figS2b_beta_p_div.pdf", width = 5, height = 4)

plot_p <-
  betadiv_p_dis %>%
  filter(!diff == '7_35') %>%
  mutate(diff = factor(diff, levels = c('7_21', '21_35'))) %>%
  ggplot(aes(x = diff, y = local_dissimilarity, fill = diff)) +
  geom_boxplot() +
  ylim(0,0.5) +
  theme_minimal() +
  theme(legend.position = 'none') +
  ggtitle('Phylogenetic beta div')

dev.off()


# 3.3- Functional
betadiv_f <- hill_func_parti_pairwise(hel, q = 1, kegg_paths, pairs = "full")

# Calculate local dissimilarity
betadiv_f_dis <-
  betadiv_f %>%
  mutate(local_dissimilarity = 1 - local_similarity) %>%
  select(site1, site2, local_dissimilarity) %>%
  write_tsv("results/tables/betadiv_f_dis.tsv")

# Compare individuals from same trial and pen between days 7-21 and 21-35
betadiv_f_dis <-
  read_tsv("results/tables/betadiv_f_dis.tsv") %>%
  rename(animal_code = 'site1') %>%
  left_join(metadata, by = 'animal_code') %>%
  left_join(metadata_2, by = 'site2') %>%
  mutate(diff = case_when(sampling_time == '7' & sampling_time_2 == '21' ~ '7_21',
                          sampling_time == '21' & sampling_time_2 == '35' ~ '21_35',
                          sampling_time == '7' & sampling_time_2 == '35' ~ '7_35')) %>%
  mutate(diff_trial = case_when(trial == 'CA' & trial_2 == 'CA' ~ 'CA',
                                trial == 'CB' & trial_2 == 'CB' ~ 'CB')) %>%
  mutate(diff_pen = case_when(pen == pen_2  ~ 'same')) %>%
  drop_na(diff, diff_trial, diff_pen)

# Stats
betadiv_f_dis %>%
  group_by(diff) %>%
  summarise(mean = mean(local_dissimilarity), sd = sd(local_dissimilarity))

# Plot
pdf("results/figures/figS2c_beta_f_div.pdf", width = 5, height = 4)

plot_f <-
  betadiv_f_dis %>%
  filter(!diff == '7_35') %>%
  mutate(diff = factor(diff, levels = c('7_21', '21_35'))) %>%
  ggplot(aes(x = diff, y = local_dissimilarity, fill = diff)) +
  geom_boxplot() +
  ylim(0,0.5) +
  theme_minimal() +
  theme(legend.position = 'none') +
  ggtitle('Functional beta div')

dev.off()


# 4- Dissimilarity tables ---------------------
# Neutral
b_n_dis_table <- read_tsv("results/tables/betadiv_n_dis.tsv") %>%
  pivot_wider(names_from = site2, values_from = local_dissimilarity) %>%
  column_to_rownames('site1') %>%
  as.dist()

# Phylogenetic
b_p_sim_table <- read_tsv("results/tables/betadiv_p_dis.tsv") %>%
  pivot_wider(names_from = site2, values_from = local_dissimilarity) %>%
  column_to_rownames('site1') %>%
  as.dist()

# Functional
b_f_sim_table <- read_tsv("results/tables/betadiv_f_dis.tsv") %>%
  pivot_wider(names_from = site2, values_from = local_dissimilarity) %>%
  column_to_rownames('site1') %>%
  as.dist()


# 5- Effect of design in microbiome composition ---------------------
perm <- how(nperm = 999)
setBlocks(perm) <- with(metadata, pen)

adonis2(betadiv_n_dis ~ trial + sex + breed + treatment * age,
        permutations = perm,
        data = metadata)
adonis2(betadiv_p_dis ~ trial + sex + breed + treatment * age,
        permutations = perm,
        data = metadata)
adonis2(betadiv_f_dis ~ trial + sex + breed + treatment * age,
        permutations = perm,
        data = metadata)


rm(list = ls())
