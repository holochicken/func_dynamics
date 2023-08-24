# Correlating increasers and increasers with BW - updated 22/03/2023

# Packages
library(vegan)
library(nlme)
library(sjPlot)
library(gridExtra)
library(tidyverse)


# 1- Load data ---------------------
mag_counts <-
  read_tsv(file = "data/metagenomics/mag_counts.tsv") %>%
  column_to_rownames(var = 'mag_id')

metadata <-
  read_tsv(file = "data/metadata/metadata.tsv") %>%
  mutate(sampling_time = factor(sampling_time, levels = c('7',
                                                          '21',
                                                          '35'))) %>%
  column_to_rownames(var = 'animal_code')

# MAG stats
stats <-
  read_tsv("data/metagenomics/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)

trends <-
  read_tsv("results/tables/hmsc_mag_trend.tsv") %>%
  select(mag_id, hmsc_trend)

# 2- Standardisation ---------------------
# By MAG genome length
mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  rownames_to_column('mag_id')


# 3- Group MAGs by trend and samples by sampling day ---------------------
# Increasers
increasers <-
  mag_weighted %>%
  left_join(trends) %>%
  filter(hmsc_trend == 'increaser') %>%
  select(!hmsc_trend) %>%
  column_to_rownames('mag_id') %>%
  t() %>%
  as.data.frame() %>%
  mutate(summed_increasers = rowSums(.)) %>%
  select(summed_increasers)

# Decreasers
decreasers <-
  mag_weighted %>%
  left_join(trends) %>%
  filter(hmsc_trend == 'decreaser') %>%
  select(!hmsc_trend) %>%
  column_to_rownames('mag_id') %>%
  t() %>%
  as.data.frame() %>%
  mutate(summed_decreasers = rowSums(.)) %>%
  select(summed_decreasers)

# Stables
stables <-
  mag_weighted %>%
  left_join(trends) %>%
  filter(hmsc_trend == 'stable') %>%
  select(!hmsc_trend) %>%
  column_to_rownames('mag_id') %>%
  t() %>%
  as.data.frame() %>%
  mutate(summed_stables = rowSums(.)) %>%
  select(summed_stables)

# Centred log ratio
abundances <- data.frame(increasers, decreasers, stables) %>%
  compositions::clr() %>%
  as.data.frame() %>%
  rownames_to_column('sample_id')


# Split samples by sampling day
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
                    ))

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
                    ))

# day 35
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
                    ))


# 4- Linear Mixed Effects Models ---------------------
# 4.1- Increasers
# Day 7 model
model_d7 <-
  lme(chicken_body_weight ~ trial + age + breed + sex + treatment + summed_increasers,
      random = ~ 1 | pen,
      data = metadata_d7)

summary(model_d7)
anova(model_d7)

plot_d7 <-
  plot_model(model_d7,
             type = 'eff',
             title = 'increasers - d7',
             terms = 'summed_increasers',
             show.data = TRUE)

# Day 21 model
model_d21 <-
  lme(chicken_body_weight ~ trial + age + breed + sex + treatment + summed_increasers,
      random = ~ 1 | pen,
      data = metadata_d21)

summary(model_d21)
anova(model_d21)

plot_d21 <-
  plot_model(model_d21,
             type = 'eff',
             title = 'increasers - d21',
             terms = 'summed_increasers',
             show.data = TRUE)

# Day 35 model
model_d35 <-
  lme(chicken_body_weight ~ trial + age + breed + sex + treatment + summed_increasers,
      random = ~ 1 | pen,
      data = metadata_d35)

summary(model_d35)
anova(model_d35)

plot_d35  <-
  plot_model(model_d35,
             type = 'eff',
             title = 'increasers - d35',
             terms = 'summed_increasers',
             show.data = TRUE)

# Refit the model with the double interactions between breed * summed_increasers + sex* summed_increasers
model_d35_breed_sex <-
  lme(chicken_body_weight ~ trial + age + treatment + breed * summed_increasers + sex* summed_increasers ,
      random = ~ 1 | pen,
      data = metadata_d35)

# Hypothesis tests
summary(model_d35_breed_sex)
anova(model_d35_breed_sex)

# Plots with sex- and breed-specific lines
plot_model(model_d35_breed_sex,
           type = 'eff',
           title = 'increasers - d35',
           terms = c('summed_increasers','sex'),
           show.data = TRUE)
plot_model(model_d35_breed_sex,
           type = 'eff',
           title = 'increasers - d35',
           terms = c('summed_increasers','breed'),
           show.data = TRUE)

## Difference in BW between lowest and highest percentiles of clr(summed_increasers)
increasers05q <- metadata_d35[metadata_d35$summed_increasers < quantile(metadata_d35$summed_increasers, probs = 0.05),]
increasers95q <- metadata_d35[metadata_d35$summed_increasers > quantile(metadata_d35$summed_increasers, probs = 0.95),]
boxplot(increasers05q$chicken_body_weight, increasers95q$chicken_body_weight)

# Difference in body weight between groups
mean(BW_increasers95q) - mean(BW_increasers05q)

t.test(increasers05q$chicken_body_weight, increasers95q$chicken_body_weight, var.equal = FALSE)

# Representation of different kinds of chicken in both groups.
data.frame(increasers05q$breed, increasers95q$breed)
data.frame(increasers05q$sex, increasers95q$sex)
data.frame(increasers05q$trial, increasers95q$trial)
data.frame(increasers05q$treatment, increasers95q$treatment)
data.frame(increasers05q$age, increasers95q$age)


#### decreasers
# Linear Mixed Effects Models
# Day 7 model
model_d7_d <-
  lme(chicken_body_weight ~ trial + age + breed + sex + treatment + summed_decreasers,
      random = ~ 1 | pen,
      data = metadata_d7)

summary(model_d7_d)
anova(model_d7_d)

plot_d7_d <-
  plot_model(model_d7,
             type = 'eff',
             title = 'decreasers - d7',
             terms = 'summed_decreasers',
             show.data = TRUE,
             dot.size = 1.5)

# Day 21 model
model_d21_d <-
  lme(chicken_body_weight ~ trial + age + breed + sex + treatment + summed_decreasers,
      random = ~ 1 | pen,
      data = metadata_d21)

summary(model_d21_d)
anova(model_d21_d)

plot_d21_d <-
  plot_model(model_d21,
             type = 'eff',
             title = 'decreasers - d21',
             terms = 'summed_decreasers',
             show.data = TRUE,
             dot.size = 1.5)

# Day 35 model
model_d35_d <-
  lme(chicken_body_weight ~ trial + age + breed + sex + treatment + summed_decreasers,
      random = ~ 1 | pen,
      data = metadata_d35)

summary(model_d35_d)
anova(model_d35_d)

plot_d35_d <-
  plot_model(model_d35,
             type = 'eff',
             title = 'decreasers - d35',
             terms = 'summed_decreasers',
             show.data = TRUE,
             dot.size = 1.5)


# 5- Plot
# Increasers
pdf(file = "results/figures/fig3_corr_bw_increasers.pdf", width = 8, height = 2.5)

grid.arrange(plot_d7, plot_d21, plot_d35, ncol = 3)

dev.off()


# Decreasers
pdf(file = "results/figures/figs3_corr_bw_decreasers.pdf", width = 8, height = 2.5)

grid.arrange(plot_d7_d, plot_d21_d, plot_d35_d, ncol = 3)

dev.off()


rm(list = ls())

# 6- Trends taxonomy ----
# Load data
trend <-
  read_tsv("results/tables/hmsc_mag_trend.tsv")

taxonomy <-
  read_tsv('data/metagenomics/taxonomy.tsv')

# Group by trend
# Decreasers
decreasers <-
  trend %>%
  filter(grepl('decreaser', hmsc_trend)) %>%
  left_join(taxonomy)

# phylum
sort(prop.table(table(decreasers$phylum)), decreasing = TRUE)
# order
sort(prop.table(table(decreasers$order)), decreasing = TRUE)


# Increasers
increasers <-
  trend %>%
  filter(grepl('increaser', hmsc_trend)) %>%
  left_join(taxonomy)

# phylum
sort(prop.table(table(increasers$phylum)), decreasing = TRUE)
# order
sort(prop.table(table(increasers$order)), decreasing = TRUE)


rm(list = ls())
