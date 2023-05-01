# Setting model structure - updated 22/03/2023

# Packages
library(Hmsc)
library(ape)
library(tidyverse)


# 1- Load data ---------------------
# Create directories
dir.create(path = "src/hmsc_mg_and_mt/mt/models", showWarnings = FALSE, recursive = TRUE)
dir.create(path = "src/hmsc_mg_and_mt/mt/model_fit", showWarnings = FALSE, recursive = TRUE)
dir.create(path = "src/hmsc_mg_and_mt/mt/panels", showWarnings = FALSE, recursive = TRUE)

# Directory
localDir = "."
MGDataDir = file.path(localDir, "data/metagenomics")
MTDataDir = file.path(localDir, "data/metatranscriptomics")
ModelDir = file.path(localDir, "src/hmsc_mg_and_mt/mt/models")

# Abundance table
mag_counts <-
  read_tsv(file = file.path(MTDataDir, "mag_expressions.tsv")) %>%
  arrange(., by_group = mag_id) %>%
  column_to_rownames(var = 'mag_id')

# Metadata
metadata <-
  read_tsv(file = "data/metadata/metadata.tsv") %>%
  filter(animal_code %in% colnames(mag_counts)) %>%
  mutate(
    sampling_time = factor(sampling_time, levels = c('Day 7',
                                                     'Day 21',
                                                     'Day 35'))
    ) %>%
  column_to_rownames(var = 'animal_code')

# Metadata summary
design <-
  metadata[, c(
    'trial', 'pen', 'age', 'breed',
    'sex', 'treatment','chicken_body_weight'
    )]

design$log_seq_depth <- log(colSums(mag_counts))

for (i in 1:ncol(design)) {
  if (is.character(design[,i])) {
    design[,i] <- factor(design[,i])
  }
}

head(design)

# MAG stats
stats <- read.delim(file.path(MTDataDir,"stats.tsv"))


# 2- Standardization and correction ---------------------
# Correction
stats$correction_factor <- median(stats$gene_length) / stats$gene_length

mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  as.data.frame()


# 3- Define tables for hmsc model ---------------------
# YData
YData <- log1p(mag_weighted)

# XData
XData <- design
mean(rownames(YData) == rownames(XData))  # Ydata and XData in the same order

# Define XFormula
XFormula_Time <- ~log_seq_depth + trial + sex + breed + treatment + age

# Phylogenetic tree
tree <- read.tree(file.path(MGDataDir, "tree.nwk"))
YData <- YData[, colnames(YData) %in% tree$tip.label]  # Filter missing MAGs in tree

# Study Design
StudyDesign <- data.frame(pen = design[,c(2)])
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
