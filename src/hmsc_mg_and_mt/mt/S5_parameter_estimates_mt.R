# Exploring parameter estimates - updated 22/03/2023

# Packages
library(Hmsc)
library(colorspace)
library(corrplot)
library(gridExtra)
library(grid)
library(factoextra)
library(tidyverse)
library(ggrepel)
library(ape)
library(ade4)
library(phylosignal)
library(phylobase)
library(abind)


# 1- Open fitted model ---------------------
# Metatranscriptomics directory
localDir = "."
DataDir = file.path(localDir, "results/clean/mags/caecum")
MTModelDir = file.path(localDir, "src/hmsc_mg_and_mt/mt/models")
MTPanelDir = file.path(localDir, "src/hmsc_mg_and_mt/mt/panels")

thin = 10
samples = 250
nChains = 4
filename = paste("model","_thin_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_", as.character(nChains),
                 ".Rdata", sep = "")
load(file = file.path(MTModelDir,filename))


# 2- Functional structure ---------------------
beta_post <- getPostEstimate(m, parName = "Beta")
plotBeta(m, beta_post, supportLevel = 0.95, plotTree = TRUE)

Gradient_age = constructGradient(m,
                                 focalVariable = "age",
                                 non.focalVariables = 1)
predY_age = predict(m, Gradient = Gradient_age, expected = TRUE)

# Example using cmag_376
plotGradient(m,
             Gradient_age,
             pred = predY_age,
             yshow = 0,
             measure = "Y",
             index = 376,
             showData = TRUE)


# MAG trends table
# Increasing MAGs
increasers <- colnames(beta_post$support)[beta_post$support[8,] > 0.95]
increasers_df <-
  data.frame(mag_id = increasers,
             parameter = beta_post$mean[8, colnames(m$Y) %in% increasers])
increasers_df$hmsc_trend_mt <- "increaser"

# Decreasing MAGs
decreasers <- colnames(beta_post$support)[beta_post$support[8,] < 0.05]
decreasers_df <-
  data.frame(mag_id = decreasers,
             parameter = beta_post$mean[8,colnames(m$Y) %in% decreasers])
decreasers_df$hmsc_trend_mt <- "decreaser"

trends <- colnames(beta_post$support) %>%
  as.data.frame() %>%
  dplyr::rename(mag_id = '.') %>%
  left_join(rbind(increasers_df, decreasers_df), by = 'mag_id') %>%
  replace_na(list(hmsc_trend_mt = 'stable')) %>%
  write_tsv("results/tables/hmsc_mag_trend_mt.tsv")


# Variance partitioning
VP <- computeVariancePartitioning(m)
plotVariancePartitioning(hM = m, VP = VP)

mpost <- convertToCodaObject(m)
quantile(unlist(mpost$Rho), probs = c(.05,.5,.95))

age_parameter <- beta_post$mean[8,]
age_parameter_phyloSorted <-
  data.frame(age_parameter = age_parameter[
    match(m$phyloTree$tip.label,
          names(age_parameter))
    ])
mean(rownames(age_parameter_phyloSorted) == m$phyloTree$tip.label)
new_tree <- m$phyloTree
new_tree$node.label <- NULL

obj <- phylo4d(new_tree, tip.data = age_parameter_phyloSorted)
barplot.phylo4d(obj)
# age.cg_MT <- phyloCorrelogram(obj, trait = "age_parameter")
# saveRDS(age.cg_MT,file="age.cg_MT.rds")
# age.cg_MT<-readRDS(file="age.cg.rds")
plot(age.cg_MT)


# 3- Compare MG and MT trends ---------------------
m_MT <- m

# Open fitted MT model
MGModelDir = file.path(localDir, "src/mags/hmsc_mg_and_mt/mg/models")
MGPanelDir = file.path(localDir, "src/mags/hmsc_mg_and_mt/mg/panels")

thin = 10
samples = 250
nChains = 4
filename = paste("model","_thin_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_",as.character(nChains),".Rdata",sep = "")
load(file = file.path(MGModelDir,filename))

m_MG <- m


# Posterior beta
beta_post_MG <- getPostEstimate(m_MG,parName = "Beta")
beta_post_MT <- getPostEstimate(m_MT,parName = "Beta")

beta_post_df <- data.frame(beta_age_MG = beta_post_MG$mean[8,],
                           beta_age_MT = beta_post_MT$mean[8,])

# Plot
pdf("results/figures/fig2d_corr_mg_mt.pdf", width = 6, height = 4)

ggplot(data = beta_post_df, mapping = aes(x = beta_age_MG, y = beta_age_MT)) +
  geom_point() +
  geom_smooth(method = "loess") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic()

dev.off()
