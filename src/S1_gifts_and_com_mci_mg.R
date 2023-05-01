# MAG gifts and community level mci values - updated 22/03/2023

# Packages
library(distillR)
library(vegan)
library(tidyverse)


# 1- Load data ---------------------
# Create directories
dir.create(path = "results/figures/", showWarnings = FALSE, recursive = TRUE)
dir.create(path = "results/tables/", showWarnings = FALSE, recursive = TRUE)

# MAG counts and metadata
mag_counts <-
  read_tsv(file = "data/metagenomics/mag_counts.tsv") %>%
  arrange(., by_group = mag_id) %>%
  column_to_rownames(var = 'mag_id')

metadata <-
  read_tsv("data/metadata/metadata.tsv") %>%
  mutate(sampling_time = factor(sampling_time, levels = c('7',
                                                          '21',
                                                          '35')))

# MAG characteristics
mag_annotations <-
  read_tsv("data/metagenomics/gene_annotations.tsv.xz")

ena_to_mag_id <-
  read_tsv("data/metadata/ena_to_mag_id.tsv")

stats <-
  read_tsv("data/metagenomics/stats.tsv") %>%
  mutate(correction_factor = median(mag_length) / mag_length)

trends <-
  read_tsv("results/tables/hmsc_mag_trend.tsv")

taxonomy_colors <-
  read_tsv("data/metadata/taxonomy_colours.tsv")

# 2- Standardisation ---------------------
# By MAG genome length
mag_weighted <-
  round(sweep(mag_counts, MARGIN = 1, stats$correction_factor, `*`), 0) %>%
  t() %>%
  data.frame()

# By sequencing depth
hel <-
  decostand(mag_weighted, 'hellinger') %>%
  t() %>%
  as.data.frame()


# 3- Distillate annotations ---------------------
gifts <- distill(mag_annotations, GIFT_db1, genomecol = 2, annotcol = c(9,10,19))

gifts_raw <-
  gifts %>%
  data.frame() %>%
  rownames_to_column(var = 'mag_name') %>%
  left_join(ena_to_mag_id %>%
              select(mag_name, mag_id),
            by = 'mag_name') %>%
  select(-mag_name) %>%
  relocate(mag_id) %>%
  write_tsv("results/tables/gifts_raw.tsv") %>%
  column_to_rownames('mag_id') %>%
  as.matrix()

genome_completeness <-
  stats %>%
  select(mag_id, completeness_score) %>%
  as.matrix()


# 4- Define completeness corretion ---------------------
# Correction function defined in
# https://github.com/anttonalberdi/DAMMA/blob/main/R/damma_correction.R
# Based on https://www.nature.com/articles/s43705-023-00221-z
gifts_correction <- function(MCI_table, genome_completeness, stats = TRUE) {
  #Sort Genomes and test matching
  MCI_table <- MCI_table[order(rownames(MCI_table)),]
  genome_completeness <- genome_completeness[order(genome_completeness[,1]),]
  if (all(as.character(rownames(MCI_table)) !=
          as.character(genome_completeness[,1])))
    stop("Pathway table is missing")
  #Get completeness values
  Genome_completeness <- as.numeric(genome_completeness[,2])
  #Create corrected fullness matrix
  MCI_table_corrected <- matrix(0,nrow = nrow(MCI_table),ncol = ncol(MCI_table))
  colnames(MCI_table_corrected) = colnames(MCI_table)
  rownames(MCI_table_corrected) = rownames(MCI_table)
  #Iterate modelling and correction for each function
  suppressWarnings(
    for (i in 1:ncol(MCI_table)) {
      Model <- glm(MCI_table[,i]~Genome_completeness,family = "binomial")
      slope_coef <- coef(Model)[2]
      if (slope_coef > 0) {
        for (j in 1:nrow(MCI_table)) {
          # Model prediction of fullness if completeness was 100%
          pred_100 <- round(predict(Model,
                                    newdata = data.frame(
                                      Genome_completeness = 100),
                                    type = "response"),1)
          # Model prediction of fullness for actual completeness of the focal MAG
          pred_focal <- round(predict(Model,
                                      newdata = data.frame(
                                        Genome_completeness = Genome_completeness[j]),
                                      type = "response"), 1)
          # The expected change in function fullness if focal MAG was 100% complete
          pred_diff <- pred_100 - pred_focal
          MCI_table_corrected[j,i] = MCI_table[j,i] + pred_diff
        }
      } else if (slope_coef <= 0) {
        MCI_table_corrected[,i] <- MCI_table[,i]
      }
    }
  )
  # If corrected fullness >1, convert it to 1.
  MCI_table_corrected[MCI_table_corrected > 1] <- 1
  # Outout overall correction statistics on screen
  total <- nrow(MCI_table)*ncol(MCI_table)
  changes <- c(MCI_table == MCI_table_corrected)
  changes <- length(changes[!(changes)])
  percentage <- round(changes / total * 100,1)
  cat(paste0(changes," out of ",total," (",percentage,"%) fullness values were corrected\n"))
  if (stats == TRUE) {
    #Outout overall correction statistics on screen
    total <- ncol(MCI_table_corrected)
    for (r in rownames(MCI_table_corrected)) {
      completeness <- round(as.numeric(genome_completeness[genome_completeness[,1] == r,2]),1)
      changes <- MCI_table[r,] == MCI_table_corrected[r,]
      changes <- length(changes[!(changes)])
      percentage <- round(changes / total * 100,1)
      cat(paste0("\t",r," (",completeness,"%): ",changes,
                 "/",total," (",percentage,"%) fullness values were corrected\n"))
    }
  }
  #Output corrected table
  return(MCI_table_corrected)
}


# 5- Apply completeness corretion ---------------------
gifts_raw <- read_tsv("results/tables/gifts_raw.tsv") %>%
  column_to_rownames('mag_id')

gifts_corrected <- gifts_correction(gifts_raw, genome_completeness)

# Aggregate into compound level
gifts_elements <- to.elements(gifts_corrected,GIFT_db1)

gifts_elements <-
  gifts_elements %>%
  as.data.frame() %>%
  rownames_to_column(var = 'mag_id') %>%
  write_tsv("results/tables/gifts_elements.tsv") %>%
  column_to_rownames(var = 'mag_id')

# Aggregate into function level
gifts_elements <- read_tsv("results/tables/gifts_elements.tsv") %>%
  column_to_rownames(var = 'mag_id')

gifts_functions <- to.functions(gifts_elements,GIFT_db1)

gifts_functions <-
  gifts_functions %>%
  as.data.frame() %>%
  rownames_to_column(var = 'mag_id') %>%
  write_tsv("results/tables/gifts_functions.tsv") %>%
  column_to_rownames(var = 'mag_id')

# Aggregate into overall domain level
gifts_domains <- to.domains(gifts_functions,GIFT_db1)

gifts_domains <-
  gifts_functions %>%
  as.data.frame() %>%
  rownames_to_column(var = 'mag_id') %>%
  write_tsv("results/tables/gifts_domains.tsv") %>%
  column_to_rownames(var = 'mag_id')


# 6- Plot all taxa ---------------------
gifts_elements <- read_tsv("results/tables/gifts_elements.tsv")

png("distillR/figs1_heatmap_all.png", width = 2480, height = 3208)

gifts_elements %>%
  pivot_longer(!mag_id, names_to = 'Code_element', values_to = 'mci') %>%
  inner_join(GIFT_db1, by = 'Code_element') %>%
  ggplot(aes(x = mag_id,
             y = Code_element,
             fill = mci)) +
  geom_tile(color = '#ffffff') +
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  facet_grid(Code_function ~ ., scales = 'free', space = 'free') +
  scale_fill_gradientn(colours = rev(terrain.colors(10))) +
  theme_void(base_size = 18) +
  theme(axis.text.y = element_text(),
        legend.position = 'top')

dev.off()


# 7- Plot selected taxa ---------------------
# gifts_elements <-  read_tsv("results/tables/gifts_elements.tsv") %>%
#   column_to_rownames(var = 'mag_id')
# Selected MAGs
selected_genomes <- c('cmag_379', 'cmag_463', 'cmag_377', 'cmag_570',
                      'cmag_348', 'cmag_062', 'cmag_510', 'cmag_815')

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

# Filter MAGs and Functions
gifts_elements_filt <-
  gifts_elements %>%
  select(mag_id, all_of(selected_gifts)) %>%
  filter(mag_id %in% selected_genomes)

# Plot
pdf("results/figures/fig2b_heatmap_selection.pdf", width = 14, height = 2)

gifts_elements_filt %>%
  pivot_longer(!mag_id, names_to = 'Code_element', values_to = 'mci') %>%
  inner_join(GIFT_db1, by = 'Code_element') %>%
  mutate(mag_id = factor(mag_id, levels = rev(c('cmag_379',
                                                'cmag_463',
                                                'cmag_377',
                                                'cmag_570',
                                                'cmag_348',
                                                'cmag_062',
                                                'cmag_510',
                                                'cmag_815'))
                         )) %>%
  ggplot(aes(x = Code_element,
             y = mag_id,
             fill = mag_id,
             group = Code_function,
             alpha = mci)) +
  geom_tile(color = '#ffffff') +
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_fill_manual(values = c('#66cada','#66cada','#99d19c','#e9ea74',
                               '#c28c5c','#f4814c','#f4814c','#f4814c'),
                    na.value = '#f4f4f4') +
  facet_grid(. ~ Code_function, scales = "free", space = "free") +
  theme_void(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'none')

dev.off()


# 8- Community level average MCIs ---------------------
elements_com <- to.community(gifts_elements,
                            sweep(hel, 2, colSums(hel), FUN = "/"),
                            GIFT_db1)

elements_com <-
  elements_com %>%
  as.data.frame() %>%
  rownames_to_column('animal_code') %>%
  write_tsv("results/tables/elements_com.tsv")

funcs_com <- to.community(gifts_functions,
                           sweep(hel, 2, colSums(hel), FUN = "/"),
                           GIFT_db1)

funcs_com <-
  funcs_com %>%
  as.data.frame() %>%
  rownames_to_column('animal_code') %>%
  write_tsv("results/tables/funcs_com.tsv")


domains_com <- to.community(gifts_domains,
                                  sweep(hel, 2, colSums(hel), FUN = "/"),
                                  GIFT_db1)

domains_com <-
  domains_com %>%
  as.data.frame() %>%
  rownames_to_column('animal_code') %>%
  write_tsv("results/tables/domains_com.tsv")



# Plot
domains_com <- read_tsv("results/tables/domains_com.tsv")

pdf("results/figures/fig1h_community_mci_time.pdf", width = 5, height = 4)

domains_com %>%
    # rownames_to_column('animal_code') %>%
    rowwise() %>%
    mutate(overall = mean(c_across(B01:D09))) %>%
    select(animal_code, overall) %>%
    left_join(metadata %>%
                select(animal_code, sampling_time),
              by = 'animal_code') %>%
    ggplot(aes(x = sampling_time,
               y = overall,
               group = sampling_time,
               colour = sampling_time)) +
    geom_boxplot(aes_string(colour = 'sampling_time', fill = 'sampling_time'),
                 width = 0.3,
                 lwd = 3,
                 outlier.color = NA,
                 position = position_nudge(x = -.4)) +
    ylim(0.26,0.38) +
    geom_jitter(width = 0.15, alpha = 0.6) +
    stat_summary(geom = 'crossbar',
                 width = 0.3,
                 fatten = 0,
                 color = 'white',
                 position = position_nudge(x = -.4),
                 fun.data = function(x){ return(c(y = median(x),
                                                  ymin = median(x),
                                                  ymax = median(x)))
                 }) +
    scale_fill_manual(values = c('#E69F00', '#CC6677', '#56B4E9')) +
    scale_color_manual(values = c('#E69F00', '#CC6677', '#56B4E9')) +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = 'none')

dev.off()


rm(list = ls())
