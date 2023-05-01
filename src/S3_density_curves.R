# Phylum densities - updated 22/03/2023

# Packages
library(vegan)
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

# MAG characteristics
stats <-
  read_tsv("data/metagenomics/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)

taxonomy <-
  read_tsv("data/metagenomics/taxonomy.tsv")

taxonomy_colors <-
  read_tsv("data/metadata/taxonomy_colours.tsv")


# 2- Standardisation ---------------------
# By MAG genome length
mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  data.frame()

# By sequencing depth
total <- decostand(mag_weighted, 'total')
rowSums(total)


# 3- Group by phylum ---------------------
sum_ind <-
  total %>%
  t() %>%
  data.frame() %>%
  rownames_to_column(var = 'mag_id') %>%
  left_join(taxonomy %>% select(mag_id, phylum), by = 'mag_id') %>%
  mutate(phylum = factor(phylum, levels = unique(taxonomy_colors$phylum))) %>%
  pivot_longer(-c(mag_id, phylum), values_to = 'value', names_to = 'animal_code') %>%
  group_by(phylum, animal_code) %>%
  summarise(total_ind = sum(value), mean_ind = mean(value))

sum_phylum <-
  sum_ind %>%
  group_by(phylum) %>%
  summarise(total_phylum = sum(total_ind), mean_phylum = mean(mean_ind))

# Filter rare phylums
count_taxa <-
  sum_ind %>%
  filter(phylum %in% c('Firmicutes_A',
                       'Firmicutes',
                       'Bacteroidota',
                       'Cyanobacteria',
                       'Proteobacteria',
                       'Actinobacteriota',
                       'Verrucomicrobiota',
                       'Campylobacterota'))

# 4- Plot density figure ---------------------
phylum_colors <- c(Campylobacterota = "#cdadca",
                   Verrucomicrobiota = "#0066ae",
                   Proteobacteria = "#60cfde",
                   Bacteroidota = "#4dc87c",
                   Cyanobacteria = "#92e09f",
                   Actinobacteriota = "#c9e0af",
                   Firmicutes = "#ffd366",
                   Firmicutes_A = "#fd8854",
                   Firmicutes_B = "#8b1222",
                   Firmicutes_C = "#5a0c37")

pdf("results/figures/fig1d_taxonomy_density.pdf", width = 7, height = 8)

ggplot(count_taxa, aes(x = total_ind, colour = phylum, fill = phylum)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = phylum_colors) +
  scale_color_manual(values = phylum_colors) +
  scale_x_log10() +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none')

dev.off()


rm(list = ls())
