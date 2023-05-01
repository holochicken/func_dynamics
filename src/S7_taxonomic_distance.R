# Tree depth - taxonomic distance - updated 22/03/2023

# Packages
library(ape)
library(tidyverse)


# 1- Load data ---------------------
# MAG characteristics
tree <-
  read.tree("data/metagenomics/tree.nwk")

taxonomy <-
  read_tsv("data/metagenomics/taxonomy.tsv")

pw_phylo_dist <- cophenetic.phylo(tree)
xy <- t(combn(colnames(pw_phylo_dist), 2))
pw_phylo_dist <- data.frame(xy, dist = pw_phylo_dist[xy])
colnames(pw_phylo_dist) <- c('mag_id.x', 'mag_id.y', 'distance')

pw_phylo_dist_taxa <-
  pw_phylo_dist %>%
  inner_join(taxonomy, by = c('mag_id.x' = 'mag_id')) %>%
  inner_join(taxonomy, by = c('mag_id.y' = 'mag_id'))


# 2- Create distance table ---------------------
tax_table <- c()
for (i in c(1:nrow(pw_phylo_dist_taxa))) {
  pair <- pw_phylo_dist_taxa[i,]
  #Phylum level
  if (!is.na(pair$phylum.x != pair$phylum.y) & pair$phylum.x != pair$phylum.y) {
    row <- c('Phylum', pair$distance)
  } else {
    #Class level
    if (!is.na(pair$class.x != pair$class.y) & pair$class.x != pair$class.y) {
      row <- c('Class', pair$distance)
    } else {
      #Order level
      if (!is.na(pair$order.x != pair$order.y) & pair$order.x != pair$order.y) {
        row <- c('Order', pair$distance)
      } else {
        #Family level
        if (!is.na(pair$family.x != pair$family.y) & pair$family.x != pair$family.y) {
          row <- c('Family', pair$distance)
        } else {
          # Genus
          if (!is.na(pair$genus.x != pair$genus.y) & pair$genus.x != pair$genus.y) {
            row <- c('Genus', pair$distance)
          }
        }
      }
    }
  }
  tax_table <- rbind(tax_table,row)
}

tax_table <-
  tax_table %>%
  data.frame() %>%
  dplyr::rename(taxonomy = 'X1', distance = 'X2') %>%
  mutate(distance = as.numeric(distance)) %>%
  write_tsv("results/tables/taxonomy_distance.tsv")


# 3- Plot ---------------------
pdf("results/figures/fig1g_taxonomy_distance.pdf", width = 6, height = 2)

tax_table %>%
ggplot(aes(x = distance,
           fill = taxonomy,
           color = taxonomy,
           alpha = 0.1
           )) +
  geom_density(adjust = 5) +
  scale_fill_manual(values = c('#E69F00', '#e28cff', '#999999',
                               '#56B4E9','#99cc00')) +
  scale_color_manual(values = c('#E69F00', '#e28cff', '#999999',
                                '#56B4E9','#99cc00')) +
  theme_void()

dev.off()


rm(list = ls())
