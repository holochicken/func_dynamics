# correlating increasers and decreasers with BW

# Packages
library(tidyverse)
library(vegan)
library(nlme)
library(sjPlot)
library(gridExtra)


# 1- Load  - ENA id-s
stats <-
  read_tsv("data/metatranscriptomics/stats.tsv") %>%
  mutate(correction_factor = median(gene_length) / gene_length)

trends <-
  read_tsv("results/tables/hmsc_mag_trend.tsv") %>%
  select(mag_id, hmsc_trend)


# Load abundance table
mag_counts <-
  read_tsv("data/metatranscriptomics/mag_expressions.tsv") %>%
  column_to_rownames("mag_id") %>%
  select(-CA13.10, -CA17.09)


# 1- Standardization and correction  -----------------
# Correction
mag_weighted <- round(
  sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  rownames_to_column("mag_id")

# increasers
increasers <-
  mag_weighted %>%
  left_join(trends) %>%
  filter(hmsc_trend == "increaser") %>%
  select(!hmsc_trend) %>%
  column_to_rownames("mag_id") %>%
  t() %>%
  as.data.frame() %>%
  mutate(summed_increasers = rowSums(.)) %>%
  select(summed_increasers)

# decreasers
decreasers <-
  mag_weighted %>%
  left_join(trends) %>%
  filter(hmsc_trend == "decreaser") %>%
  select(!hmsc_trend) %>%
  column_to_rownames("mag_id") %>%
  t() %>%
  as.data.frame() %>%
  mutate(summed_decreasers = rowSums(.)) %>%
  select(summed_decreasers)


# stable
stables <-
  mag_weighted %>%
  left_join(trends) %>%
  filter(hmsc_trend == "stable") %>%
  select(!hmsc_trend) %>%
  column_to_rownames("mag_id") %>%
  t() %>%
  as.data.frame() %>%
  mutate(summed_stables = rowSums(.)) %>%
  select(summed_stables)


abundances <- data.frame(increasers, decreasers, stables) %>%
  compositions::clr() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")


# Load metadata and split samples by sampling day
metadata <-
  read_tsv(file = 'data/metadata/metadata.tsv') %>%
  filter(animal_code %in% colnames(mag_counts)) %>%
  mutate(
    sampling_time = factor(sampling_time, levels = c('7', '21', '35'))
  ) %>%
  column_to_rownames(var = 'animal_code')

# day 7
metadata_d7 <-
  as.data.frame(metadata %>%
                  rownames_to_column('sample_id') %>%
                  left_join(abundances) %>%
                  column_to_rownames('sample_id') %>%
                  filter(sampling_time == '7') %>%
                  mutate(
                    age = as_factor(age),
                    pen = as_factor(pen),
                    trial = as_factor(trial)
                  )) %>%
  mutate(chicken_body_weight = log(chicken_body_weight))

# day 21
metadata_d21 <-
  as.data.frame(metadata %>%
                  rownames_to_column('sample_id') %>%
                  left_join(abundances) %>%
                  column_to_rownames('sample_id') %>%
                  filter(sampling_time == '21') %>%
                  mutate(
                    age = as_factor(age),
                    pen = as_factor(pen),
                    trial = as_factor(trial)
                  )) %>%
  mutate(chicken_body_weight = log(chicken_body_weight))


metadata_d35 <-
  as.data.frame(metadata %>%
                  rownames_to_column('sample_id') %>%
                  left_join(abundances) %>%
                  column_to_rownames('sample_id') %>%
                  filter(sampling_time == '35') %>%
                  mutate(
                    age = as_factor(age),
                    pen = as_factor(pen),
                    trial = as_factor(trial)
                  )) %>%
  mutate(chicken_body_weight = log(chicken_body_weight))

# Linear Mixed Effects Models
model_d7 <- lme(chicken_body_weight ~  age + trial + sex + breed + treatment + summed_increasers,
                random = ~ 1 | pen,
                data = metadata_d7)
summary(model_d7)
# anova(model_d7)
a <- plot_model(model_d7,
                type = 'eff',
                title = 'increasers - d7',
                terms = 'summed_increasers',
                show.data = TRUE)

model_d21 <- lme(chicken_body_weight ~  age + trial + sex + breed + treatment + summed_increasers,
                 random = ~ 1 | pen,
                 data = metadata_d21)
summary(model_d21)
# anova(model_d21)
b <- plot_model(model_d21,
                type = 'eff',
                title = 'increasers - d21',
                terms = 'summed_increasers',
                show.data = TRUE)

model_d35 <- lme(chicken_body_weight ~ age + trial + sex + breed + treatment + summed_increasers,
                 random = ~ 1 | pen,
                 data = metadata_d35)
summary(model_d35)

# anova(model_d35)
c <- plot_model(model_d35,
                type = 'eff',
                title = 'increasers - d35',
                terms = 'summed_increasers',
                show.data = TRUE)

# Plot
# step 6- Theme options
# pdf(
#   file = 'C:/Users/smarcos007/Desktop/ratio_bw_v4.pdf',
#   width = 8,
#   height = 2.5)

grid.arrange(a, b, c, ncol = 3)

# dev.off()
