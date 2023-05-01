# Functional landscape - updated 22/03/2023

# Packages
library(Rtsne)
library(RColorBrewer)
library(tidyverse)


# 1- Load data ---------------------
# MAG characteristics
gifts <-
  read_tsv("results/tables/gifts_elements.tsv") %>%
  column_to_rownames(var = 'mag_id')

taxonomy <-
  read_tsv("data/metagenomics/taxonomy.tsv")

trends <-
  read_tsv("results/tables/hmsc_mag_trend.tsv")

taxonomy_colors <-
  read_tsv("data/metadata/taxonomy_colours.tsv")


# 2- Tsne ---------------------
set.seed(100)
tsne_func <- Rtsne(X = gifts, dims = 2, check_duplicates = FALSE)
avg_mci <- rowMeans(gifts)
trends$avg_mci <- avg_mci

tsne_df <-
  tsne_func$Y %>%
  as.data.frame() %>%
  mutate(mag_id = rownames(gifts)) %>%
  left_join(trends, by = 'mag_id') %>%
  left_join(taxonomy, by = 'mag_id') %>%
  mutate_at(vars(phylum, class, order, hmsc_trend), factor) %>%
  mutate(order = factor(order, levels = taxonomy_colors$order)) %>%
  mutate(phylum = factor(phylum, levels = unique(taxonomy_colors$phylum))) %>%
  rename(tsne1 = 'V1', tsne2 = 'V2')


# 3- Plots ---------------------
# Order
pdf("results/figures/fig1a_tsine_order.pdf", width = 6, height = 4)

tsne_df %>%
  ggplot(aes(x = tsne1, y = tsne2, color = order)) +
  geom_point(size = 2, shape = 16, alpha = 0.8) +
  scale_color_manual(values = taxonomy_colors$color_order) +
  theme_minimal() +
  theme(legend.position = 'none')

dev.off()


# Phylum
pdf("results/figures/fig1a_tsine_phylum.pdf", width = 6, height = 4)

tsne_df %>%
  ggplot(aes(x = tsne1, y = tsne2, color = phylum)) +
  geom_point(size = 2, shape = 16, alpha = 0.8) +
  scale_color_manual(values = unique(taxonomy_colors$color_phylum)) +
  theme_minimal() +
  theme(legend.position = 'none')

dev.off()


# DistillR MCI value
pdf("results/figures/fig1b_tsine_mci.pdf", width = 6, height = 4)

tsne_df %>%
  ggplot(aes(x = tsne1, y = tsne2, color = avg_mci)) +
  geom_point(size = 2, shape = 16, alpha = 0.8) +
  scale_colour_gradientn(colours = rev(terrain.colors(10))) +
  theme_minimal() +
  theme(legend.position = 'none')

dev.off()


#Trend
pdf("results/figures/fig1c_tsine_trend.pdf", width = 6, height = 4)

tsne_df %>%
  ggplot(aes(x = tsne1, y = tsne2, color = hmsc_trend)) +
  geom_point(size = 2, shape = 16, alpha = 0.8) +
  scale_color_manual(values = c('#b85454','#9bb854','#cecece')) +
  theme_minimal() +
  theme(legend.position = 'none')

dev.off()


rm(list = ls())
