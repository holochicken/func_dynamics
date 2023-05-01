# Setting model structure - updated 22/03/2023

# Packages
library(Hmsc)
library(ape)
library(tidyverse)


# 1- Load data ---------------------
# Create directories
dir.create(path = "src/hmsc_mg_and_mt/mg/models", showWarnings = FALSE, recursive = TRUE)
dir.create(path = "src/hmsc_mg_and_mt/mg/model_fit", showWarnings = FALSE, recursive = TRUE)
dir.create(path = "src/hmsc_mg_and_mt/mg/panels", showWarnings = FALSE, recursive = TRUE)

# Directory
localDir = "."
DataDir = file.path(localDir, "data/metagenomics")
ModelDir = file.path(localDir, "src/hmsc_mg_and_mt/mg/models")

# Abundance table
mag_counts <-
  read_tsv(file.path(DataDir, "mag_counts.tsv")) %>%
  arrange(., by_group = mag_id) %>%
  column_to_rownames(var = 'mag_id')

# Metadata
metadata <-
  read_tsv(file = "data/metadata/metadata.tsv") %>%
  mutate(
    sampling_time = factor(sampling_time, levels = c('Day 7',
                                                     'Day 21',
                                                     'Day 35'))
    ) %>%
  column_to_rownames(var = 'animal_code')

# MAG stats
stats <- read_tsv(file.path(DataDir, "stats.tsv"))


# 2- Standardization and correction ---------------------
# Correction
stats$correction_factor <- median(stats$mag_length) / stats$mag_length

mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  as.data.frame()

# metadata summary
design <- metadata[, c(
  'trial', 'pen', 'age','breed',
  'sex', 'treatment', 'chicken_body_weight'
)]

design$log_seq_depth <- log(colSums(mag_counts))
design$chicken_ID <- rownames(design)

for (i in 1:ncol(design)) {
  if (is.character(design[,i])) {
    design[,i] <- factor(design[,i])
  }
}

head(design)


# 3- Define tables for hmsc model ---------------------
# YData
YData <- log1p(mag_weighted)

# XData
XData <- design
mean(rownames(YData) == rownames(XData))  # Ydata and XData in the same order

# Define XFormula
XFormula_Time <- ~log_seq_depth + trial + sex + breed + treatment + age

# Phylogenetic tree
tree <- read.tree(file.path(DataDir,"tree.nwk"))
YData <- YData[,colnames(YData) %in% tree$tip.label]  # Filter missing MAGs in tree

# Study Design
StudyDesign <- design[,c(2,9)]
rL.Pen <- HmscRandomLevel(units = levels(StudyDesign$pen))


# 4- Model ---------------------
m <- Hmsc(Y = YData,
          XData = XData,
          XFormula = XFormula_Time,
          studyDesign = StudyDesign,
          ranLevels = list("pen" = rL.Pen),
          phyloTree = tree,
          distr = "normal",
          YScale = TRUE)

save(m, file = file.path(ModelDir,"Unfitted_model.Rdata"))
