####################
### BBBC022 Data ###
####################

#
# check the difference between the set of features for 10,000 and 200,000
#

library(magrittr)
library(dplyr)
library(tidyverse)
library(stringr)

filename <- "../../../../input/BBBC022/Pf_Gustafsdottir.rds"

Pf <- readRDS(filename)
Pf <- Pf$data %>% 
  dplyr::rename(Metadata_Plate = Plate,
                Metadata_Well = Well,
                Metadata_broad_sample = Image_Metadata_BROAD_ID,
                Metadata_cpd_name = Image_Metadata_SOURCE_COMPOUND_NAME)

variables <-
  colnames(Pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")

metadata <-
  colnames(Pf) %>% str_subset("^Metadata_")

f1 <- "../../input/BBBC022_2013/single_cells_findCorrelation/feat_names/selected_features_10000cells_findCorr.Rda"
f2 <- "../../input/BBBC022_2013/single_cells_findCorrelation/feat_names/selected_features_200000cells_findCorr.Rda"

load(file.path(f2))
f2 <- variables.fs.sc.fc

# select the features
Pf2 <- Pf %>% select(one_of(metadata, f2))

load(file.path(f1))
f1 <- variables.fs.sc.fc

variables <- str_replace(variables, "Hoechst", "DNA")
variables <- str_replace(variables, "Ph_golgi", "AGP")
variables <- str_replace(variables, "Syto", "RNA")
cl <- colnames(Pf)
cl <- str_replace(cl, "Hoechst", "DNA")
cl <- str_replace(cl, "Ph_golgi", "AGP")
cl <- str_replace(cl, "Syto", "RNA")
colnames(Pf) <- cl

Pf1 <- Pf %>% select(one_of(metadata, f1))

cor1 <- cor(Pf1 %>% select(one_of(f1)) %>% as.matrix %>% t)
cor2 <- cor(Pf2 %>% select(one_of(f2)) %>% as.matrix %>% t)

cor1 <- cor1[lower.tri(cor1, diag = FALSE)]
cor2 <- cor2[lower.tri(cor2, diag = FALSE)]


library(lsa)
cosine(cor1, cor2)


