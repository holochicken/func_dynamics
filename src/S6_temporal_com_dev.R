# Temporal changes and sources of variation on the microbiota - updated 22/03/2023


# packages
library(vegan)
library(tidyverse)


# 1- Load data ---------------------
# MAG counts and metadata
mag_counts <-
  read_tsv(file = "data/metagenomics/mag_counts.tsv") %>%
  column_to_rownames(var = 'mag_id')

metadata <-
  read_tsv(file = "data/metadata/metadata.tsv") %>%
  mutate(sampling_time = factor(sampling_time, levels = c('7',
                                                          '21',
                                                          '35'))) %>%
  column_to_rownames(var = 'animal_code')

# MAG characteristics
stats <-
  read_tsv("data/metagenomics/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)

taxonomy <-
  read_tsv("data/metagenomics/taxonomy.tsv")


# 2- Standardisation ---------------------
# By MAG genome length
mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  data.frame()

# By sequencing depth
hel <- decostand(mag_weighted, 'hellinger')

# 3- Permanova  ---------------------
perm <- how(nperm = 999)
setBlocks(perm) <- with(metadata, pen)

adonis2(hel ~ trial * sampling_time + breed * sampling_time +
          sex * sampling_time+ treatment * sampling_time,
        permutations = perm,
        data = metadata)

adonis2(hel ~ trial * sampling_time + breed + sex + treatment,
        permutations = perm,
        data = metadata)


# 4- Distance-based RDA ---------------------
set.seed(4)
hel_bray <- vegdist(hel, method = 'bray')
hel_bray_sqrt <- sqrt(hel_bray)

pcoa <- cmdscale(hel_bray_sqrt, k = nrow(hel) - 1, eig = TRUE)
pcoa_scores <- pcoa$points

db_rda <- rda(pcoa_scores ~ sampling_time, data = metadata)

# Stats
anova(db_rda, by = 'term')
summary(db_rda)
RsquareAdj(db_rda)

db_rda_scores <-
  data.frame(
    scores(db_rda, display = 'wa'),
    pen = metadata$pen,
    time = metadata$sampling_time
    )

# Plot
pdf(file = "results/figures/fig1f_db_rda.pdf", width = 7, height = 4)

db_rda_scores %>%
  ggplot(aes(x = RDA1, y = RDA2, colour = time)) +
  geom_point() +
  scale_color_manual(values = c('#e6a024', '#cc6777', '#5bb4e5')) +
  theme_bw() +
  theme(legend.position = 'none')

dev.off()


rm(list = ls())
