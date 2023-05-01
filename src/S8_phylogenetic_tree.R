# Radial plot - updated 22/03/2023

# Packages
library(vegan)
library(phyloseq)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(tidytree)
library(ggnewscale)
library(RColorBrewer)
library(plotly)
library(tidyverse)


# 1- Load and tidy data ---------------------
# 1.1- Relative abundance table
mag_counts <-
  read_tsv("data/metagenomics/mag_counts.tsv") %>%
  arrange(., by_group = mag_id) %>%
  column_to_rownames(var = 'mag_id')

# 1.2- Metadata
metadata <-
  read_tsv(file = "data/metadata/metadata.tsv") %>%
  select(animal_code, sampling_time)

# 1.3- MAG stats
stats <-
  read_tsv("data/metagenomics/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)

# 1.4- Standardisation
# By MAG genome length
mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  data.frame()

# By sequencing depth
total <- decostand(mag_weighted, 'total')
rowSums(total) # check if all values are 1

# 1.5- Phylogenetic tree
tree <-
  read.tree("data/metagenomics/tree.nwk")

# 1.6- Taxonomy
taxonomy <-
  read_tsv("data/metagenomics/taxonomy.tsv") %>%
  rename(label = 'mag_id') %>%
  left_join(taxonomy_colors %>% select(!phylum), by = 'order')

# 1.7- Joining taxonomy to tree
tax_tree <-
  tree %>%
  as_tibble() %>%
  left_join(taxonomy, by = 'label') %>%
  mutate(order = factor(order, levels = taxonomy_colors$order)) %>%
  mutate(phylum = factor(phylum, levels = unique(taxonomy_colors$phylum))) %>%
  as.treedata()

phylum_ids <-
  taxonomy %>%
  group_by(phylum) %>%
  summarise(phylum_node_id = MRCA(tree, label)) %>%
  rename(phylum_id = 'phylum')

# 1.8- Trends
trends <-
  read_tsv("results/tables/hmsc_mag_trend.tsv")

# 1.9- Compounds table - obtain mean values
mci_com <-
  read_tsv(file = "results/tables/gifts_elements.tsv") %>%
  mutate(mci = rowMeans(across(where(is.numeric)))) %>%
  select(mag_id, mci) %>%
  left_join(stats) %>%
  left_join(trends) %>%
  select(mag_id, mci, completeness_score, hmsc_trend) %>%
  mutate(hmsc_trend = factor(hmsc_trend, levels = c('increaser',
                                                    'stable',
                                                    'decreaser')))

# 1.10- Relative abundance table - separated by day
mag_counts_day <-
  total %>%
  rownames_to_column('animal_code') %>%
  left_join(metadata) %>%
  pivot_longer(!c(animal_code, sampling_time),
               names_to = 'mag_id', values_to = 'values') %>%
  group_by(sampling_time, mag_id) %>%
  summarise(total = colMeans(across(where(is.numeric)))) %>%
  arrange(., by_group = mag_id) %>%
  mutate(
    sampling_time = factor(sampling_time, levels = c('7',
                                                     '21',
                                                     '35')))

# clean environmnet
rm(stats, tree, mag_counts, total, metadata, decreasers, increasers)

# 2- Radial plot ---------------------
# Define colors
phylum_colors <- c(Halobacteriota = "#f2f2f2",
                   Methanobacteriota = "#bfbfbf",
                   Patescibacteria = "#dacce3",
                   Campylobacterota = "#cdadca",
                   Synergistota = "#c08eb3",
                   Thermoplasmatota = "#777dae",
                   Verrucomicrobiota = "#0066ae",
                   Desulfobacterota = "#68b0dc",
                   Proteobacteria = "#60cfde",
                   Bacteroidota = "#4dc87c",
                   Cyanobacteria = "#92e09f",
                   Actinobacteriota = "#c9e0af",
                   Deferribacterota = "#dff77e",
                   Firmicutes = "#ffd366",
                   Firmicutes_A = "#fd8854",
                   Firmicutes_B = "#8b1222",
                   Firmicutes_C = "#5a0c37")

order_colors <- c(Methanomicrobiales = "#f2f2f2",
                  Methanobacteriales = "#bfbfbf",
                  Saccharimonadales = "#dacce3",
                  Campylobacterales = "#cdadca",
                  Synergistales = "#c08eb3",
                  Methanomassiliicoccales = "#777dae",
                  Victivallales = "#0066ae",
                  Opitutales = "#1c7ebc",
                  Verrucomicrobiales = "#4a96cc",
                  Desulfovibrionales = "#68b0dc",
                  Enterobacterales = "#60cfde",
                  RF32 = "#60dfd2",
                  Burkholderiales = "#5cdfb5",
                  'Rs-D84' = "#40df91",
                  Bacteroidales = "#4dc87c",
                  Flavobacteriales = "#88c88b",
                  Gastranaerophilales = "#92e09f",
                  Coriobacteriales = "#c9e0af",
                  Actinomycetales = "#d8e093",
                  Deferribacterales = "#dff77e",
                  RF39 = "#ecf76d",
                  Lactobacillales = "#feef68",
                  'ML615J-28' = "#fde671",
                  Erysipelotrichales = "#ffd366",
                  Acholeplasmatales = "#fdc151",
                  Bacillales = "#fc953d",
                  RFN20 = "#fd8035",
                  Oscillospirales = "#fd8854",
                  TANB77 = "#f4814d",
                  Christensenellales = "#c26340",
                  Lachnospirales = "#c28c5c",
                  Clostridiales = "#ce9360",
                  Monoglobales_A = "#ce734f",
                  Peptostreptococcales = "#c65631",
                  Monoglobales = "#b93725",
                  UBA1212 = "#b10f19",
                  UBA4068 = "#8b1222",
                  Selenomonadales = "#7a1a12",
                  Acidaminococcales = "#64121f",
                  Veillonellales = "#5a0c37"
)

# 2.1- Basic circular tree
p <-
  ggtree(tax_tree, layout = 'fan', open.angle = 15, size = 0.2)

# 2.2- Rotate, open, highlight and label tree
p1 <-
  rotate_tree(p, 105) +
  geom_highlight(data = phylum_ids,
                 mapping = aes(node = phylum_node_id, fill = phylum_id),
                 show.legend = FALSE,
                 align = 'right',
                 alpha = 0.5,
                 to.bottom = TRUE
                ) +
  scale_fill_manual(values = phylum_colors)

# 2.3- Colour orders
p2 <-
  p1 +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = mag_id, fill = order),
    offset = 0.06,
    width = 0.1,
    grid.params = list(color = NA)) +
  scale_fill_manual(values = order_colors)

# 2.4- MCI value
p3 <-
  p2 +
  new_scale_fill() +
  geom_fruit(data = mci_com,
             geom = geom_tile,
             mapping = aes(y = mag_id, fill = mci),
             offset = 0.06,
             width = 0.1,
             grid.params = list(color = NA)) +
  scale_fill_gradientn(colours = brewer.pal(7, 'YlGnBu'), name = 'MCI')

# 2.5- Relative abundance (by sampling day)
p4 <-
  p3 +
  new_scale_fill() +
  geom_fruit(data = mag_counts_day,
             geom = geom_tile,
             mapping = aes(y = mag_id, x = sampling_time, fill = total),
             color = NA,
             offset = 0.03,
             width = 0.1,
             pwidth = 0.1
             ) +
  scale_fill_gradientn(colours = brewer.pal(7, 'Greys'), name = 'Abundance')

# 2.6- Temporal trend
p5 <-
  p4 +
  new_scale_fill() +
  geom_fruit(data = mci_com,
             geom = geom_tile,
             mapping = aes(y = mag_id, fill = hmsc_trend),
             offset = 0.06,
             width = 0.1,
             grid.params = list(color = NA)) +
  scale_fill_manual(values = c('#a8d48c','#737373', '#d06a6b'),name = 'Trend')

# 2.7- Theme options
pdf("results/figures/fig2a_radial_plot.pdf", width = 8, height = 6)

p5  +
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.position = 'none',
        legend.direction = 'horizontal')

dev.off()


rm(list = ls())
