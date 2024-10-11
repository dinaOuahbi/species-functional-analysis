
---
title: "Immunolife : Analyse Metagenomique Fonctionnelle"
author: "DO"
output:
  pdf_document: default
keep_tex: true
html_document: default
date: "2024-08-27"
---
# 703 obs / 1990 msp

library(ggplot2)
library(ggpubr)
library(tidyr)
library(ComplexHeatmap)
library(vegan)
library(ape)
library(momr)
library(kableExtra)
library(wesanderson)
library(stringr)
library(dplyr)
library(readxl)
library(svglite)
library(data.table)
library(grid)
library(gridExtra)
library(cowplot) # or library(patchwork)

# FUNCTIONS  --------------------------------------------------------------------------
source('~/projects/mgp_tools/scripts/ord_to_df.R') # convert ordinations (pcoa, nmds...) to ggplottable df
source('~/projects/mgp_tools/scripts/plot_es.R') # plots the effect size, abundance and prevalence (v2)
source('~/projects/mgp_tools/scripts/m4_taxonomy_shrink.R') # shrink taxonomy
source('~/projects/mgp_tools/scripts/readRData.R') # read .RData format
source('~/projects/mgp_tools/scripts/merge_msp_profiles.R')# merge msp oral/gut profiles
source('~/projects/immunolife/scripts/utilities/legends_N.R') # add N to legends

kegg <- read.delim("/projects/programs/FAnToMet/all_modules_definition_GMM_GBM_KEGG_89.tsv")

# TAXONUMY --------------------------------------------------------------------------
msp_taxonomy_old <- readRData('/projects/biodatabank/catalogue/all_catalogue/human_MSP_set/latest/MSP_set_ref_oss_gss_status_612_ref_paper_20240627_gtdb_220.RData')
msp_taxonomy_gut <- read.delim('/projects/biodatabank/catalogue/hs_10_4_igc2/mgs/latest/taxonomy/latest/hs_10_4_1990_MSP.20240626.gtdb_r220_annot_long.tsv')
msp_taxonomy_oral <- read.delim('/projects/biodatabank/catalogue/hs_8_4_oral/mgs/latest/taxonomy/latest/hs_8_4_oral_853_MSP.20240627.gtdb_r220_annot_long.tsv')

msp_taxonomy_all <- rbind(
  mutate(msp_taxonomy_gut,ref_MSP_test='gut'),
  mutate(msp_taxonomy_oral[!msp_taxonomy_oral$msp_name %in%msp_taxonomy_gut$msp_name, ], ref_MSP_test='oral')
) 
table(msp_taxonomy_all$ref_MSP_test)
msp_taxonomy_all$ref_MSP <- msp_taxonomy_old$ref_MSP[match(msp_taxonomy_all$msp_name,msp_taxonomy_old$msp_id)]

# PROLFILES --------------------------------------------------------------------------
meteor_pro_list <- 
  list(
    gut =  read.delim("/projects/shrcan/immunolife/data_9/genes_1/mgs_1/immunolife_vs_hs_10_4_igc2_smart_shared_reads_29032024.downHQ.20000000.fpkm.mgs.100.tsv", header = TRUE, comment.char = '#', sep = '\t', row.names = 1),
    oral =  read.delim("/projects/shrcan/immunolife/data_10/genes_1/mgs_1/immunolife_vs_hs_8_4_oral_smart_shared_reads_05042024.1.downHQ.20000000.fpkm.mgs.100.tsv", header = TRUE, comment.char = '#', sep = '\t', row.names = 1)
  )

meteor_pro_list$oral <- meteor_pro_list$oral[,match(colnames(meteor_pro_list$gut), colnames(meteor_pro_list$oral))]
meteor_pro_list_raw <- 
  merge_msp_profiles(gut = meteor_pro_list$gut, oral = meteor_pro_list$oral, taxonomy = msp_taxonomy_all, save = FALSE, msp_id = 'msp_name')
meteor_pro_raw <- meteor_pro_list_raw$profile
rownames(meteor_pro_raw) <- meteor_pro_raw$id_msp
meteor_pro_raw$id_msp <- NULL
meteor_pro_raw <- as.data.frame(t(meteor_pro_raw))


# COMPLETUDE --------------------------------------------------------------------------
chemin_gut = "/shrcan/immunolife/data_9/genes_1/mgs_1/taxonomy_1/ftmt_1/completude_ind_mgs_mod_data_downHQ_20000000fpkm_ftmt.txt"
chemin_oral = "/shrcan/immunolife/data_10/genes_1/mgs_1/taxonomy_1/ftmt_1/completude_ind_mgs_mod_data_downHQ_20000000fpkm_ftmt.txt"
seuil_completude = 0.9

completude = list(gut = setNames(as.data.frame(fread(chemin_gut, header = F, sep = "\t"), stringsAsFactors = F)[,c(1:4)], c("sample", "module", "msp", "complet")),
                  oral = setNames(as.data.frame(fread(chemin_oral, header = F, sep = "\t"), stringsAsFactors = F)[,c(1:4)], c("sample", "module", "msp", "complet")))

completude = rbind(completude$gut, completude$oral)
completude = completude[completude$complet >= seuil_completude,] 
completude$concat = paste(completude$sample, completude$module, completude$msp, sep = "__")
completude = completude[!duplicated(completude$concat),]
completude = completude[completude$sample %in% rownames(meteor_pro_raw),]
completude$sample = factor(completude$sample)
completude$module = factor(completude$module)
completude$msp = factor(completude$msp)
completude$concat_sample_module = factor(paste(completude$sample, completude$module, sep = "__"))
completude$concat_sample_msp = factor(paste(completude$sample, completude$msp, sep = "__"))



# Fragmentation de la matrice en plusieurs matrice chacune centré sur un échantillons
modules = setNames(as.data.frame(matrix(0, nrow = length(unique(completude$sample)), ncol = length(unique(completude$module))), row.names = unique(completude$sample), stringAsFactors = F), unique(completude$module))

for (i in levels(completude$sample)){
  print(i)
  a = completude[completude$sample ==i,] # table des modules d'un echantillon donnée
  
  for (j in unique(a$module)) # pour chaque module de cet echantillon
  {
    modules[i,j] = sum(meteor_pro_raw[i, # echantillon en question
                                                           as.character(a[a$module == j,"msp"]) # les msp d'un module
                                                           ]
                       )
  }
}

# CLINICAL DATA
# CLMD PARAMETERS
clmd_retro <-'/projects/immunolife/analysis/data/metadata/rhu_lumiere/n_481_LUNG_TOPO_LONGITUDINAL_20240820.xlsx'
clmd_pro <- '/projects/immunolife/analysis/data/metadata/immunolife/imnlf1_clmd_E20240424_M20240429.tsv'
clmd_list <- 
  list(
    # retro = read.delim(parameters$paths$clmd_retro, header = TRUE, as.is = TRUE, stringsAsFactors = FALSE),
    retro = as.data.frame(read_excel(clmd_retro, sheet = 1)),
    pro = read.delim(clmd_pro, header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
  )
clmd_list$pro$IDMGP <- clmd_list$pro$IDMGP_INRAE

if(!dim(clmd_list$retro)[1] == 482){
  warning('incorrect retrospective clinical data dimensions: should be 482, while it is:', dim(clmd_list$retro)[1])
}

clin <- as.data.frame(clmd_list$pro)[, c('IDMGP_INRAE','ATB_before_sample','ATB_60_42_merge')]
clin <- clin %>% distinct(IDMGP_INRAE, .keep_all = TRUE) 
clin <- drop_na(clin)
rownames(clin) <- clin$IDMGP_INRAE
clin$IDMGP_INRAE <- NULL
merged_df <- merge(clin, modules, by = "row.names", all = FALSE)
merged_df$ATB_before_sample <- as.factor(merged_df$ATB_before_sample)
merged_df$ATB_60_42_merge <- as.factor(merged_df$ATB_60_42_merge)




# Calcul des statistiques descriptives par groupe
summary_stats <- merged_df %>%
  group_by(ATB_60_42_merge) %>%
  summarise(across(starts_with("M"), list(mean = mean, sd = sd), .names = "{col}_{fn}"))


# Test t pour chaque module
generate_plot <- function(module, group="ATB_60_42_merge") {
  pvalue <- paste0("t-test : ", t_test_results$p_value[which(t_test_results$Module == module)])
  grob <- grobTree(textGrob(pvalue, x = 0.3, y = 0.80, hjust = 0,
                            gp = gpar(col = "black", fontsize = 8, fontface = "italic")))

  metab <- kegg$Glycolysis..Embden.Meyerhof.pathway...glucose....pyruvate[which(kegg$M00001 == module)]
  p <- ggplot(merged_df, aes(x = ATB_60_42_merge, y = get(module), fill = ATB_60_42_merge)) +
    geom_boxplot(notch = TRUE, outlier.colour = "red", outlier.shape = 8) +
    labs(title = paste0(metab), y = module, x = group) +
    geom_jitter(shape = 16, position = position_jitter(0.2)) +
    theme_minimal() +
    theme(legend.position = "top") +
    scale_fill_manual(values = c("#999999", "#E69F00")) +
    annotation_custom(grob)
  return(p)
}
t_test_results <- data.frame(Module = character(), p_value = numeric())
for (module in grep("^M", colnames(merged_df), value = TRUE)) {
  t_test <- t.test(as.formula(paste(module, "~ ATB_60_42_merge")), data = merged_df)
  t_test_results <- rbind(t_test_results, data.frame(Module = module, p_value = t_test$p.value))
}
plots <- lapply(t_test_results$Module[which(t_test_results$p_value<0.05)], generate_plot)
combined_plot <- plot_grid(plotlist = plots, ncol = 4)
print(combined_plot)


##########################################################
generate_plot2 <- function(module, group="ATB_before_sample") {
  pvalue <- paste0("t-test : ", t_test_results$p_value[which(t_test_results$Module == module)])
  grob <- grobTree(textGrob(pvalue, x = 0.3, y = 0.80, hjust = 0,
                            gp = gpar(col = "black", fontsize = 8, fontface = "italic")))

  metab <- kegg$Glycolysis..Embden.Meyerhof.pathway...glucose....pyruvate[which(kegg$M00001 == module)]
  p <- ggplot(merged_df, aes(x = ATB_before_sample, y = get(module), fill = ATB_before_sample)) +
    geom_boxplot(notch = TRUE, outlier.colour = "red", outlier.shape = 8) +
    labs(title = paste0(metab), y = module, x = group) +
    geom_jitter(shape = 16, position = position_jitter(0.2)) +
    theme_minimal() +
    theme(legend.position = "top") +
    scale_fill_manual(values = c("#999999", "#E69F00")) +
    annotation_custom(grob)
  return(p)
}
t_test_results <- data.frame(Module = character(), p_value = numeric())
for (module in grep("^M", colnames(merged_df), value = TRUE)) {
  t_test <- t.test(as.formula(paste(module, "~ ATB_before_sample")), data = merged_df)
  t_test_results <- rbind(t_test_results, data.frame(Module = module, p_value = t_test$p.value))
}

plots <- lapply(t_test_results$Module[which(t_test_results$p_value<0.05)], generate_plot2)
combined_plot <- plot_grid(plotlist = plots, ncol = 4)
print(combined_plot)





# HEATMAP
library(pheatmap)
library(grid)
# Préparation des données pour la heatmap
module_data <- merged_df %>% select(starts_with("M"))
row.names(module_data) <- merged_df$Row.names

# Heatmap
annotation_col <- data.frame(Group = as.factor(merged_df$ATB_60_42_merge))
rownames(annotation_col) <- rownames(module_data)
annotation_colors <- list(Group = c("0" = "blue","1" = "red"))
pheatmap(module_data, 
         annotation_col = annotation_col,#data.frame(Group = merged_df$ATB_60_42_merge),
         annotation_colors = annotation_colors, 
         scale = "row")


# Vérifier les colonnes constantes
constant_columns <- sapply(module_data, function(x) length(unique(x)) == 1)
null_columns <- sapply(module_data, function(x) all(is.na(x)))
module_data_clean <- module_data[, !(constant_columns | null_columns)]
pca_results <- prcomp(module_data_clean, scale. = TRUE)
scores <- as.data.frame(pca_results$x)
scores$ATB_before_sample <- merged_df$ATB_before_sample  # Ajoutez la variable de groupe pour la coloration
p1 <- ggplot(scores, aes(x = PC1, y = PC2, color = ATB_before_sample)) +
  geom_point() +
  labs(title = "PCA des modules fonctionnels", x = "PC1", y = "PC2") +
  theme_minimal()

scores$ATB_60_42_merge <- merged_df$ATB_60_42_merge  # Ajoutez la variable de groupe pour la coloration
p2 <- ggplot(scores, aes(x = PC1, y = PC2, color = ATB_60_42_merge)) +
  geom_point(label = rownames(scores)) +
  labs(title = "PCA des modules fonctionnels", x = "PC1", y = "PC2") +
  theme_minimal()

grid.arrange(p1, p2, ncol = 2)









