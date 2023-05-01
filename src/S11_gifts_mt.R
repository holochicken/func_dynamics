# Quantitative MAG gifts - updated 19/01/2023

# Packages
library(vegan)
library(distillR)
library(data.table)
library(RColorBrewer)
library(tidyverse)


# 1- Load data ---------------------
# MAG expression
gene_expr <-
  read_tsv("data/metatranscriptomics/gene_expressions.tsv.xz") %>%
  select(1:129) %>%
  rename(mag_name = 'MAG') %>%
  relocate(mag_name)

# Chunk analysis to sets of 100 MAGs
mags <-  sort(unique(gene_expr$mag_name))

gene_expr_1 <- gene_expr[gene_expr$mag_name %in% mags[c(1:100)],c(2:129)]
gene_expr_2 <- gene_expr[gene_expr$mag_name %in% mags[c(101:200)],c(2:129)]
gene_expr_3 <- gene_expr[gene_expr$mag_name %in% mags[c(201:300)],c(2:129)]
gene_expr_4 <- gene_expr[gene_expr$mag_name %in% mags[c(301:400)],c(2:129)]
gene_expr_5 <- gene_expr[gene_expr$mag_name %in% mags[c(401:500)],c(2:129)]
gene_expr_6 <- gene_expr[gene_expr$mag_name %in% mags[c(501:600)],c(2:129)]
gene_expr_7 <- gene_expr[gene_expr$mag_name %in% mags[c(601:700)],c(2:129)]
gene_expr_8 <- gene_expr[gene_expr$mag_name %in% mags[c(701:825)],c(2:129)]

save(gene_expr_1, file = "results/tables/gene_expr_1.Rdata")
save(gene_expr_2, file = "results/tables/gene_expr_2.Rdata")
save(gene_expr_3, file = "results/tables/gene_expr_3.Rdata")
save(gene_expr_4, file = "results/tables/gene_expr_4.Rdata")
save(gene_expr_5, file = "results/tables/gene_expr_5.Rdata")
save(gene_expr_6, file = "results/tables/gene_expr_6.Rdata")
save(gene_expr_7, file = "results/tables/gene_expr_7.Rdata")
save(gene_expr_8, file = "results/tables/gene_expr_8.Rdata")

rm(list = ls())
gc()


# 2- Distillation ---------------------
# MAG gene annotations
mag_ann <-
  read_tsv("data/metagenomics/gene_annotations.tsv.xz") %>%
  mutate(gene_length = end_position - start_position)

# Chunk analysis
load("results/tables/gene_expr_1.Rdata")
distq_1 <- distillq(gene_expr_1, mag_ann, GIFT_db1, genecol = 1, genomecol = 2, annotcol = c(9, 10, 19))
save(distq_1, file = "results/tables/distq_1.Rdata")
rm(gene_expr_1, distq_1)

load("results/tables/gene_expr_2.Rdata")
distq_2 <- distillq(gene_expr_2, mag_ann, GIFT_db1, genecol = 1, genomecol = 2, annotcol = c(9, 10, 19))
save(distq_2, file = "results/tables/distq_2.Rdata")
rm(gene_expr_2, distq_2)

load("results/tables/gene_expr_3.Rdata")
distq_3 <- distillq(gene_expr_3, mag_ann, GIFT_db1, genecol = 1, genomecol = 2, annotcol = c(9, 10, 19))
save(distq_3, file = "results/tables/distq_2.Rdata")
rm(gene_expr_3, distq_3)

load("results/tables/gene_expr_4.Rdata")
distq_4 <- distillq(gene_expr_4, mag_ann, GIFT_db1, genecol = 1, genomecol = 2, annotcol = c(9, 10, 19))
save(distq_4, file = "results/tables/distq_4.Rdata")
rm(gene_expr_4, distq_4)

load("results/tables/gene_expr_5.Rdata")
distq_5 <- distillq(gene_expr_5, mag_ann, GIFT_db1, genecol = 1, genomecol = 2, annotcol = c(9, 10, 19))
save(distq_5, file = "results/tables/distq_5.Rdata")
rm(gene_expr_5, distq_5)

load("results/tables/gene_expr_6.Rdata")
distq_6 <- distillq(gene_expr_6, mag_ann, GIFT_db1, genecol = 1, genomecol = 2, annotcol = c(9, 10, 19))
save(distq_6, file = "results/tables/distq_6.Rdata")
rm(gene_expr_6, distq_6)

load("results/tables/gene_expr_7.Rdata")
distq_7 <- distillq(gene_expr_7, mag_ann, GIFT_db1, genecol = 1, genomecol = 2, annotcol = c(9, 10, 12))
save(distq_7, file = "results/tables/distq_7.Rdata")
rm(gene_expr_7, distq_7)

load("results/tables/gene_expr_8.Rdata")
distq_8 <- distillq(gene_expr_8, mag_ann, GIFT_db1, genecol = 1, genomecol = 2, annotcol = c(9, 10, 12))
save(distq_8, file = "results/tables/distq_8.Rdata")
rm(gene_expr_8, distq_8)


# 3- Tidy data ---------------------
# Load chunks
load("results/tables/distq_1.Rdata")
load("results/tables/distq_2.Rdata")
load("results/tables/distq_3.Rdata")
load("results/tables/distq_4.Rdata")
load("results/tables/distq_5.Rdata")
load("results/tables/distq_6.Rdata")
load("results/tables/distq_7.Rdata")
load("results/tables/distq_8.Rdata")

# Join chunks
distq_raw <- c(distq_1, distq_2, distq_3, distq_4,
               distq_5, distq_6, distq_7, distq_8)

rm(distq_1, distq_2, distq_3, distq_4, distq_5, distq_6, distq_7, distq_8)

# Tidy table
distq <-
 lapply(distq_raw, function(x) {
  rownames(x) <- gsub("F1a","",rownames(x))
  return(x)
})

#Filter relevant samples
mag_counts <-
  read_tsv("data/metagenomics/mag_counts.tsv")

distq <- lapply(distq, function(x) {
  x <- x[intersect(colnames(mag_counts),rownames(x)),]
  return(x)
})


# 4- Aggregated community -----
#Merge all MAGs
distq_merged <- Reduce("+", distq)

#Aggregate to compounds
distq_merged_compounds <- to.elements(distq_merged, GIFT_db)
write.table(distq_merged_compounds,"results/tables/community_mci_expression.tsv", col.names = T, row.names = T, sep = "\t")

# Transformation
giftqs_elements_clr <- compositions::clr(t(distq_merged_compounds))
write.table(giftqs_elements_clr,"results/tables/community_MCI_expression_clr.tsv", row.names = T, col.names = T, sep = "\t", quote = F)


# 5- Community level average stats -----
distq_merged_compounds <- distq_merged_compounds[sort(rownames(distq_merged_compounds)),]

compounds_table_df <- reshape2::melt(distq_merged_compounds)
colnames(compounds_table_df) <- c("Samples", "Compounds", "MCI")

compounds_table_df2 <- merge(compounds_table_df, pathway_db, by.x = "Compounds", by.y = "Code_compound")

metadata <-
  read_tsv("data/metadata/metadata.tsv") %>%
  column_to_rownames('animal_code')

compounds_table_df2_merged <- merge(compounds_table_df2, metadata, by.x = "Samples", by.y = "row.names")

#Average
avg_caecum <- compounds_table_df2_merged %>%
  summarise(mean = mean(MCI),
            std = sd(MCI))

#Average per time point
compounds_table_df2_merged %>%
  filter(trial != "CC") %>%
  group_by(sampling_time) %>%
  summarise(mean = mean(MCI), std = sd(MCI))

compounds_table_df2_merged %>%
  group_by(sampling_time, trial) %>%
  summarise(mean = mean(MCI), std = sd(MCI))

#Total community MCI vs. bodyweight at day 35
compounds_table_df2_merged %>%
  filter(sampling_time == 35) %>%
  select(Samples, MCI, chicken_body_weight, sex, breed, trial) %>%
  group_by(Samples) %>%
  mutate(MCI_total = sum(MCI)) %>%
  select(Samples, chicken_body_weight, MCI_total, sex, breed, trial) %>%
  distinct() %>%
  ggplot(., aes(x = chicken_body_weight, y = log(MCI_total), color = sex)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y~x)

#Function grpup vs. bodyweight at day 35
compounds_table_df2_merged %>%
  filter(sampling_time == 35) %>%
  #filter(Type == "Biosynthesis") %>%
  filter(Type == "Degradation") %>%
  select(Samples, MCI, chicken_body_weight, sex, breed, trial) %>%
  group_by(Samples) %>%
  mutate(MCI_total = sum(MCI)) %>%
  select(Samples, chicken_body_weight, MCI_total, sex, breed, trial) %>%
  distinct() %>%
  ggplot(., aes(x = chicken_body_weight, y = log(MCI_total), color = breed)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y~x)

#Specific function MCI vs. bodyweight at day 35
compounds_table_df2_merged %>%
  filter(sampling_time == 35) %>%
  filter(Code_function == "B08") %>%
  select(Samples, MCI, chicken_body_weight, sex, breed, trial) %>%
  group_by(Samples) %>%
  mutate(MCI_total = sum(MCI)) %>%
  select(Samples, chicken_body_weight, MCI_total, sex, breed, trial) %>%
  distinct() %>%
  ggplot(., aes(x = chicken_body_weight, y = log(MCI_total), color = trial)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y~x)

#Heatmap
ggplot(compounds_table_df2_merged,
       aes(x = Compounds, y = Samples, fill = log(MCI), group = Code_function)
       ) +
  geom_tile() +
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_fill_gradientn(colours = brewer.pal(7, "YlGnBu")) +
  facet_grid(. ~ Code_function, scales = "free", space = "free") +
  theme_grey(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

rm(list = ls())
gc()


# 6- Individual bacteria -----
load("results/tables/distq_1.RData")
load("results/tables/distq_2.RData")
load("results/tables/distq_3.RData")
load("results/tables/distq_4.RData")
load("results/tables/distq_5.RData")
load("results/tables/distq_6.RData")
load("results/tables/distq_7.RData")
load("results/tables/distq_8.RData")

distq_raw <- c(distq_1, distq_2, distq_3, distq_4,
               distq_5, distq_6, distq_7, distq_8)

rm(distq_1, distq_2, distq_3, distq_4, distq_5, distq_6, distq_7, distq_8)

# Tidy table
distq <-
  lapply(distq_raw, function(x) {
    rownames(x) <- gsub("F1a","",rownames(x))
    return(x)
  })

metadata <-
  read_tsv("data/metadata/metadata.tsv") %>%
  mutate(sampling_time = factor(sampling_time, levels = c('7',
                                                          '21',
                                                          '35')))
metadata$Sample <- metadata$animal_code

mag_counts <-
  read_tsv("data/metagenomics/mag_counts.tsv")

colnames(mag_counts)[1] <- "Sample"

#Subset expression data
distq_mag <- distq$ERR4836887_bin.54 #cmag_379
distq_mag <- distq$ERR4836949_bin.36 #cmag_463
distq_mag <- distq$ERR4836887_bin.35 #cmag_377
distq_mag <- distq$ERR4968588_bin.12 #cmag_570
distq_mag <- distq$ERR4835994_bin.32 #cmag_348
distq_mag <- distq$ERR4303161bin.49  #cmag_062
distq_mag <- distq$ERR4968597_bin.38 # cmag_590
distq_mag <- distq$ERR4836964_bin.21 #cmag_510

#Subset relative abundance data
#Normalise by average relative abundance
mag_counts <-
  read_tsv(file = "data/metagenomics/mag_counts.tsv") %>%
  arrange(., by_group = mag_id) %>%
  column_to_rownames(var = 'mag_id')

stats <-
  read_tsv("data/metagenomics/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)

mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  data.frame()

# By sequencing depth
hel <-
  decostand(mag_weighted, 'hellinger') %>%
  t() %>%
  as.data.frame()

# select mag
# cmag_379, cmag_463, cmag_377, cmag_570, cmag_348, cmag_062, cmag_590, cmag_510
hel_mag <- as.data.frame(t(hel["cmag_379",]))
hel_mag$Sample <- rownames(hel_mag)
colnames(hel_mag)[1] <- "RelAbun"

#Distill compounds
distq_mag_compounds <- distill_compounds(distq_mag, pathway_db)

#Plot mt curve
pdf("results/figures/mg_curve.pdf", width = 8, height = 4)

distq_mag_compounds %>%
  reshape2::melt() %>%
  rename(Sample = Var1, Code_compound = Var2, MCI = value) %>%
  inner_join(pathway_db, by = "Code_compound") %>%
  mutate(Sample = gsub("F1a", "", Sample)) %>%
  inner_join(metadata, by = "Sample") %>%
  inner_join(hel_mag, by = "Sample") %>%
  group_by(Sample) %>%
  mutate(MCI_total = sum(MCI)) %>%
  select(Sample, chicken_body_weight, MCI_total, sampling_time, sex, breed, trial, RelAbun) %>%
  filter(trial != "CC") %>%
  #mutate(sampling_time = factor(sampling_time, levels = c(7,21,35))) %>%
  distinct() %>%
  ggplot(., aes(x = sampling_time, y = RelAbun)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  theme_void()
dev.off()

#Plot mt curve
pdf("results/figures/mt_curve.pdf", width = 8, height = 4)

distq_mag_compounds %>%
  reshape2::melt() %>%
  rename(Sample = Var1, Code_compound = Var2, MCI = value) %>%
  inner_join(pathway_db, by = "Code_compound") %>%
  mutate(Sample = gsub("F1a", "", Sample)) %>%
  inner_join(metadata, by = "Sample") %>%
  inner_join(hel_mag, by = "Sample") %>%
  group_by(Sample) %>%
  mutate(MCI_total = sum(MCI)) %>%
  select(Sample, chicken_body_weight, MCI_total, sampling_time, sex, breed, trial, RelAbun) %>%
  filter(trial != "CC") %>%
  #mutate(sampling_time = factor(sampling_time, levels = c(7,21,35))) %>%
  distinct() %>%
  ggplot(., aes(x = sampling_time, y = log(MCI_total))) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  theme_void()
dev.off()


# 7- Pretty heatmap selected taxa ------
load("results/tables/distq_1.RData")
load("results/tables/distq_2.RData")
load("results/tables/distq_3.RData")
load("results/tables/distq_4.RData")
load("results/tables/distq_5.RData")
load("results/tables/distq_6.RData")
load("results/tables/distq_7.RData")
load("results/tables/distq_8.RData")

distq_raw <- c(distq_1, distq_2, distq_3, distq_4,
               distq_5, distq_6, distq_7, distq_8)

rm(distq_1, distq_2, distq_3, distq_4, distq_5, distq_6, distq_7, distq_8)

distq <-
  lapply(distq_raw, function(x) {
    rownames(x) <- gsub("F1a","",rownames(x))
    return(x)
  })

distq_compounds <- lapply(distq_raw, function(x) to.elements(x, GIFT_db))

distq_compounds2 <- lapply(distq, function(x) x[-grep('^CC', row.names(x)),])
distq_compounds3 <- lapply(distq_compounds, function(x) colMeans(x))
distq_compounds_merged <- as.data.frame(do.call(rbind, distq_compounds3))

selected_genomes <- c("ERR4836887_bin.54",
                      "ERR4836949_bin.36",
                      "ERR4836887_bin.35",
                      "ERR4968588_bin.12",
                      "ERR4835994_bin.32",
                      "ERR4303161bin.49",
                      "ERR4836964_bin.21",
                      "ERR7167063_bin.29")

selected_gifts <- c(
  # Nucleic acid biosynthesis
  "B0101","B0102","B0103","B0104","B0105","B0106",
  # Amino acid biosynthesis
  "B0204","B0205","B0206","B0207","B0208","B0209","B0210","B0211","B0212",
  "B0213","B0214","B0215","B0216","B0217","B0218","B0220","B0221",
  # Amino acid derivative biosynthesis
  "B0303","B0307","B0309","B0310",
  # SCFA biosynthesis
  "B0401","B0402","B0403",
  # Organic anion biosynthesis
  "B0601","B0602","B0603","B0604","B0605",
  # Vitamin biosynthesis
  "B0701","B0702","B0703","B0704","B0705","B0706","B0707","B0708","B0709",
  "B0710","B0711","B0712",
  # Aromatic compound biosynthesis
  "B0802","B0803","B0804",
  # Antibiotic biosynthesis
  "B1004",
  # Lipid degradation
  "D0101","D0102","D0103","D0104",
  # Polysaccharide degradation
  "D0201","D0202","D0203","D0204","D0205","D0206","D0207","D0208","D0209",
  "D0210","D0211","D0212","D0213",
  # Sugar degradation
  "D0301","D0302","D0304","D0305","D0306","D0307","D0308","D0309","D0310",
  # Amino acid degradation
  "D0501","D0502","D0504","D0505","D0506","D0507","D0508","D0509","D0510",
  "D0511","D0512","D0513","D0516","D0518",
  # Nitrogen compound degradation
  "D0601","D0602","D0604","D0606","D0607","D0609","D0610","D0612","D0613",
  # Alcohol degradation
  "D0701","D0702","D0704","D0705","D0706","D0708",
  # Xenobiotic degradation
  "D0807","D0816","D0817",
  # Antibiotic degradation
  "D0901","D0902","D0904","D0905","D0907","D0908","D0910","D0911"
)

distq_compounds_merged_selected <- t(distq_compounds_merged[selected_genomes, selected_gifts])

#Normalise by average relative abundance
mag_counts <-
  read_tsv(file = "data/metagenomics/mag_counts.tsv") %>%
  arrange(., by_group = mag_id) %>%
  column_to_rownames(var = 'mag_id')

stats <-
  read_tsv("data/metagenomics/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)

mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  data.frame()

# By sequencing depth
hel <-
  decostand(mag_weighted, 'hellinger') %>%
  t() %>%
  as.data.frame()

hel_mean <- rowMeans(hel)
hel_mean_selected <- hel_mean[selected_genomes]
distq_compounds_merged_selected_corrected <- sweep(distq_compounds_merged_selected, 2, hel_mean_selected, "/")

# Plot
pdf("distillR/heatmap_selection_expression.pdf", width = 14, height = 2)

distq_compounds_merged_selected_corrected %>%
  reshape2::melt() %>%
  rename(Code_element = Var1, Genome = Var2, MCI = value) %>%
  inner_join(GIFT_db, by = "Code_element") %>%
  mutate(Genome = factor(Genome, levels = rev(c("ERR4836887_bin.54","ERR4836949_bin.36",
                                                "ERR4836887_bin.35","ERR4968588_bin.12",
                                                "ERR4835994_bin.32","ERR4303161bin.49",
                                                "ERR4836964_bin.21","ERR7167063_bin.29"))
                         )) %>%
  ggplot(., aes(x = Code_element, y = Genome, fill = Genome, group = Code_function, alpha = log(MCI,2)))+
  geom_tile(color = "#ffffff") +
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_fill_manual(values = c("#66cada",
                               "#66cada",
                               "#99d19c",
                               "#e9ea74",
                               "#c28c5c",
                               "#f4814c",
                               "#f4814c",
                               "#f4814c"),
                    na.value = "#f4f4f4") +
  facet_grid(. ~ Code_function, scales = "free", space = "free") +
  theme_void(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")

dev.off()
