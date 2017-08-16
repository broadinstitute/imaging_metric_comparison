library(magrittr)
library(dplyr)
library(tidyverse)
library(stringr)

###############
# Find Correlation single cell Pearson
###############

### Find correlation single cells 10 vs 200 -> 0.998706
hit_10 <- readRDS("../../input/BBBC022_2013/single_cells_findCorrelation/hit_selected/Pearson/10000/Hit_pearson_seed_1346.rds")
variables_10 <-
  colnames(hit_10) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")

hit_200 <- readRDS("../../input/BBBC022_2013/selected_single_cell_zoom/10000/hit_selected/Pearson/BBBC022_fs_SVD_350_seed_9149.rds")
variables_200 <-
  colnames(hit_200) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")

length(intersect(variables_10, variables_200)) # 186 features in common

# find intersection of the two sets of compounds
cmpd10 <- hit_10$Metadata_broad_sample # 892 unique compounds
cmpd200 <- hit_200$Metadata_broad_sample # 894

cmpds <- intersect(cmpd10, cmpd200) # 882 compounds

length(cmpds)

cor10 <- cor(hit_10 %>% filter(Metadata_broad_sample %in% cmpds) %>% select(one_of(variables_10)) %>% as.matrix %>% t)
cor200 <- cor(hit_200 %>% filter(Metadata_broad_sample %in% cmpds) %>% select(one_of(variables_200)) %>% as.matrix %>% t)

cor10 <- cor10[lower.tri(cor10, diag = FALSE)]
cor200 <- cor200[lower.tri(cor200, diag = FALSE)]


library(lsa)
cosine(cor10, cor200)

# Find Correlation single cell ##############
######### 10,000 vs 200,000 cells 
# Pearson
# 186 features in common out of 373/386
# ~10 compounds different: 877 common compounds
# cosine distance: 0.998706
# Jaccard
# 186 features in common out of 373/386
# 872 common compounds
# 0.9983422
######### 10,000 Jaccard versus Pearson
# 382 features in common out of 386
# 764 common compounds out of ~900
# cosine of 1
######## 200,000 Jaccard versus Pearson
# 369 features in common out of 373
# 752 common compounds
# cosine of 1

# find Correlation versus SVD-entropy (350 features) on single cells
######## 10,000 Jaccard
# 248 features in common out of 354/386
# 803 commpounds in common
# 0.9682487
####### 10,000 Pearson
# 248 features in common out of 354/386
# 854 commpounds in common
# 0.973675