# Exploring parameter estimates - updated 22/03/2023

# Packages
library(Hmsc)
library(colorspace)
library(corrplot)
library(gridExtra)
library(grid)
library(factoextra)
library(tidyverse)
library(ape)
library(ade4)
library(phylosignal)
library(phylobase)


# 1- Open fitted model ---------------------
# Directory
localDir = "."
DataDir = file.path(localDir, "data/metagenomics")
ModelDir = file.path(localDir, "src/hmsc_mg_and_mt/mg/models")
PanelDir = file.path(localDir, "src/hmsc_mg_and_mt/mg/panels")

# Chain characteristics
thin = 10
samples = 250
nChains = 4
filename = paste("model","_thin_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_", as.character(nChains),
                 ".Rdata", sep = "")
load(file = file.path(ModelDir,filename))


# 2- Functional structure ---------------------
beta_post <- getPostEstimate(m, parName = "Beta")
plotBeta(m, beta_post,supportLevel = 0.95, plotTree = TRUE)

Gradient_age = constructGradient(m,
                                 focalVariable = "age",
                                 non.focalVariables = 1)
predY_age = predict(m, Gradient = Gradient_age, expected = TRUE)

# Example using cmag_376
plotGradient(m, Gradient_age,
             pred = predY_age,
             yshow = 0,
             measure = "Y",
             index = 376,
             showData = TRUE)

# MAG trends table
# Increasing MAGs
increasers <- colnames(beta_post$support)[beta_post$support[8,] > 0.95]
increasers_df <- data.frame(mag_id = increasers,
                            parameter = beta_post$mean[4,colnames(m$Y) %in% increasers])
increasers_df$hmsc_trend <- "increaser"
# Decreasing MAGs
decreasers <- colnames(beta_post$support)[beta_post$support[8,] < 0.05]
decreasers_df <- data.frame(mag_id = decreasers,
                            parameter = beta_post$mean[4,colnames(m$Y) %in% decreasers])
decreasers_df$hmsc_trend <- "decreaser"

trends <- colnames(beta_post$support) %>%
  as.data.frame() %>%
  dplyr::rename(mag_id = '.') %>%
  left_join(rbind(increasers_df, decreasers_df), by = 'mag_id') %>%
  replace_na(list(hmsc_trend = 'stable')) %>%
  write_tsv("results/tables/hmsc_mag_trend.tsv")

# Variance partitioning
VP <- computeVariancePartitioning(m)
plotVariancePartitioning(hM = m, VP = VP)

# Phylogenetic signal
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

age.cg <- phyloCorrelogram(obj, trait = "age_parameter")
saveRDS(age.cg,file = "age.cg.rds")
age.cg <- readRDS(file = "age.cg.rds")

plot(age.cg)

