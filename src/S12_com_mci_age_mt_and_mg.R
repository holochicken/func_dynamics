# Comparing Community MCI values for MG and MT

# Packages
library(nlme)
library(tidyverse)


# 1- Load data ---------------------
# Community mci - MG
elements_com_mg <-
  read_tsv("results/tables/elements_com.tsv") %>%
  t() %>%
  as.data.frame()

metadata <-
  read_tsv("data/metadata/metadata.tsv") %>%
  mutate(sampling_time = factor(sampling_time, levels = c('7',
                                                          '21',
                                                          '35'))) %>%
  mutate(trial_b = case_when(trial == 'CB' ~ 1,
                             trial == 'CA' ~ 0))

# Community mci - MT
elements_com_mt <-
  read_tsv("results/tables/elements_com_mt.tsv") %>%
  t() %>%
  as.data.frame()

metadata_mt <-
  metadata %>%
  filter(animal_code %in% rownames(elements_com_mt))


# 2- Modelling ---------------------
# MG

p_val <- vector(mode = "numeric", length = ncol(elements_com_mg))
parameter <- vector(mode = "numeric", length = ncol(elements_com_mg))
p_val_adj <- vector(mode = "numeric", length = ncol(elements_com_mg))

results_mg <- data.frame(parameter, p_val, p_val_adj)
rownames(results_mg) <- colnames(elements_com_mg)

for (i in 1:ncol(elements_com_mg)) {
  if ((sum(elements_com_mg[,i] > 0) > 4) & (var(elements_com_mg[,i]) > 0)) {
    comp <- elements_com_mg[,i]
    # glmm
    m <- glmmTMB(comp ~ trial + sex + breed + treatment + age + (1 | pen) + (1 | Sample_ID),
                 family = binomial, data = metadata)
    # coeff
    results_mg[i,1] <- summary(m)$coefficients$cond[7,1]
    results_mg[i,2] <- summary(m)$coefficients$cond[7,4]
  }else{
    results_mg[i,1] <- NA
    results_mg[i,2] <- NA
    results_mg[i,3] <- NA
  }
}

results_mg$p_val_adj <- p.adjust(results_mg$p_val, method = "fdr")
results_mg <- results_mg[!is.na(results_mg$parameter),]

results_mg[results_mg$p_val < 0.05,]
results_mg[results_mg$p_val_adj < 0.05,]

dim(results_mg)


# MT
p_val <- vector(mode = "numeric", length = ncol(elements_com_mt))
parameter <- vector(mode = "numeric", length = ncol(elements_com_mt))
p_val_adj <- vector(mode = "numeric", length = ncol(elements_com_mt))

results_mt <- data.frame(parameter, p_val, p_val_adj)
rownames(results_mt) <- colnames(elements_com_mt)

for (i in 1:ncol(elements_com_mt)) {
  if ((sum(elements_com_mt[,i] > 0) > 4) & (var(elements_com_mt[,i]) > 0)) {
    comp <- elements_com_mt[,i]
    # lme
    m <- lme(comp ~ trial + sex + breed + treatment + age, random = ~1|pen, data = metadata_mt)
    # summ
    results_mt[i,1] <- summary(m)$tTable[7,1]
    results_mt[i,2] <- summary(m)$tTable[7,5]
  }else{
    results_mt[i,1] <- NA
    results_mt[i,2] <- NA
    results_mt[i,3] <- NA
  }
}

results_mt$p_val_adj <- p.adjust(results_mt$p_val, method = "fdr")
results_mt <- results_mt[!is.na(results_mt$parameter),]

results_mt[results_mt$p_val < 0.05,]
results_mt[results_mt$p_val_adj < 0.05,]

dim(results_mt)

# 3- Compare MG and MT based results ---------------------
results_mg_comp <- results_mg[rownames(results_mg) %in% rownames(results_mt),]
results_mt_comp <- results_mt[rownames(results_mt) %in% rownames(results_mg),]
results_mt_comp <- results_mt_comp[match(rownames(results_mg_comp), rownames(results_mt_comp)),]
mean(rownames(results_mg_comp) == rownames(results_mt_comp))

plot(results_mg_comp$parameter,results_mt_comp$parameter)
abline(lm(results_MT_comp$parameter ~ results_Mmgcomp$parameter))

abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
abline(0, 1, lty = 2)

beta_mci_df <- data.frame(mci = rownames(results_mg_comp),
                          mg_age_beta = results_mg_comp$parameter,
                          mt_age_beta = results_mt_comp$parameter)

# Plot
pdf("results/figures/fig2e_corr_mci_mg_mt.pdf", width = 6, height = 4)

ggplot(data = beta_mci_df, mapping = aes(x = mg_age_beta, y = mt_age_beta)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  geom_text_repel(aes(label = mci), size = 2)

dev.off()
