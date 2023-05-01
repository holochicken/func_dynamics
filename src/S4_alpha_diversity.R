# Alpha diversity metrics - updated 22/03/2023

# Packages
library(vegan)
library(ape)
library(hillR)
library(nlme)
library(sjPlot)
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

log_seq_depth <-
  as.data.frame(log(colSums(mag_counts))) %>%
  rename(seq_depth = 'log(colSums(mag_counts))') %>%
  cbind(animal_code = metadata$animal_code, .)


# 3- Diversity indices ---------------------
# Neutral
alphadiv_n <-
  hill_taxa(hel, q = 1) %>%
  as.data.frame() %>%
  rename('n' = '.') %>%
  rownames_to_column('animal_code')

# Functional
alphadiv_f <-
  hill_func(hel, kegg_paths, q = 1 , fdis = FALSE) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('animal_code') %>%
  select(animal_code, q = Q)

# Phylogenetic
hel_filtered <-
  hel %>%
  select(!c(cmag_274, cmag_430, cmag_554))

alphadiv_p <-
  hill_phylo(hel_filtered, tree, q = 1) %>%
  as.data.frame() %>%
  rename('p' = '.') %>%
  rownames_to_column('animal_code')


div_all <-
  alphadiv_n %>%
  left_join(alphadiv_p, by = 'animal_code') %>%
  left_join(alphadiv_f, by = 'animal_code') %>%
  left_join(metadata, by = 'animal_code') %>%
  left_join(log_seq_depth, by = 'animal_code') %>%
  write_tsv("results/tables/alphadiv_all.tsv")

# 4- Diversity plots ---------------------
# Neutral
pdf("results/figures/fig1e_alpha_n_div.pdf", width = 5, height = 4)

p_n <- ggplot(div_all,
       aes(x = sampling_time,
           y = n,
           group = sampling_time,
           colour = sampling_time)
       ) +
  geom_boxplot(aes_string(colour = 'sampling_time', fill = 'sampling_time'),
               width = 0.3,
               lwd = 3,
               outlier.color = NA,
               position = position_nudge(x = -.4)
               ) +
  ylim(200,600) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  stat_summary(geom = 'crossbar',
               width = 0.3,
               fatten = 0,
               color = 'white',
               position = position_nudge(x = -.4),
               fun.data = function(x){ return(c(y = median(x),
                                                ymin = median(x),
                                                ymax = median(x))
               ) }) +
  scale_fill_manual(values = c('#E69F00', '#CC6677', '#56B4E9')) +
  scale_color_manual(values = c('#E69F00', '#CC6677', '#56B4E9')) +
  theme_minimal() +
  ggtitle("Neutral diversity") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

dev.off()

# Phylogenetic
pdf("results/figures/fig1e_alpha_p_div.pdf", width = 5, height = 4)

p_p <- ggplot(div_all,
       aes(x = sampling_time,
           y = p,
           group = sampling_time,
           colour = sampling_time)
       ) +
  geom_boxplot(aes_string(colour = 'sampling_time', fill = 'sampling_time'),
               width = 0.3,
               lwd = 3,
               outlier.color = NA,
               position = position_nudge(x = -.4)
               ) +
  ylim(9,21) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  stat_summary(geom = "crossbar",
               width = 0.3,
               fatten = 0,
               color = "white",
               position = position_nudge(x = -.4),
               fun.data = function(x){ return(c(y = median(x),
                                                ymin = median(x),
                                                ymax = median(x))
               ) }) +
  scale_fill_manual(values = c('#E69F00', '#CC6677', '#56B4E9')) +
  scale_color_manual(values = c('#E69F00', '#CC6677', '#56B4E9')) +
  theme_minimal() + ylab("Effective number of strains") +
  ggtitle("Phylogenetic diversity") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = 'none')

dev.off()

# Functional
pdf("results/figures/fig1e_alpha_f_div.pdf", width = 5, height = 4)

p_f <- ggplot(div_all,
       aes(x = sampling_time,
           y = q,
           group = sampling_time,
           colour = sampling_time)
       ) +
  geom_boxplot(aes_string(colour = 'sampling_time', fill = 'sampling_time'),
               width = 0.3,
               lwd = 3,
               outlier.color = NA,
               position = position_nudge(x = -.4)
  ) +
  ylim(22,27) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  stat_summary(geom = "crossbar",
               width = 0.3,
               fatten = 0,
               color = "white",
               position = position_nudge(x = -.4),
               fun.data = function(x){ return(c(y = median(x),
                                                ymin = median(x),
                                                ymax = median(x))
               ) }) +
  scale_fill_manual(values = c("#E69F00", "#CC6677", "#56B4E9")) +
  scale_color_manual(values = c("#E69F00", "#CC6677", "#56B4E9")) +
  theme_minimal() + ylab("Effective number of strains") +
  ggtitle("Functional diversity") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

dev.off()

grid.arrange(p_n, p_p, p_f, ncol = 3)
# 5- Linear mixed models - temporal development ---------------------
# Neutral
m_div_n <-
  lme(n ~ seq_depth + sex + breed + treatment + trial * age,
      random = ~1|pen,
      data = alpha_diversity)

summary(m_div_n)
anova(m_div_n)
plot(m_div_n)

# Phylogenetic
m_div_p <-
  lme(p ~ seq_depth + sex + breed + treatment + trial * age,
      random = ~1|pen,
      data = alpha_diversity)

summary(m_div_p)
anova(m_div_p)
plot(m_div_p)

# Functional
m_div_f <-
  lme(q ~ seq_depth + sex + breed + treatment + trial * age,
      random = ~1|pen,
      data = alpha_diversity)

summary(m_div_f)
anova(m_div_f)
plot(m_div_f)

plot_model(m_div_n,type = 'eff',terms = c('age'))
plot_model(m_div_p,type = 'eff',terms = c('age'))
plot_model(m_div_f,type = 'eff',terms = c('age'))


rm(list = ls())
