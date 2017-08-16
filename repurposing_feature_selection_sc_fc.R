library(magrittr)
library(dplyr)
library(tidyverse)
library(stringr)

# feature selection using single cell data and findCorrelation on the Repurposing dataset

# import data
pf.sc <- readRDS(file.path("../../input/repurposing/dmso_single_cell/Repurposing_dmso_sc.rds")) %>% scale() %>% as.data.frame() #20,491x1777

pf <- readRDS(file.path("../../input/repurposing/repurposing_normalized.rds"))
# remove features with NA
pf <- pf[ , colSums(is.na(pf)) == 0]
variables <-
  colnames(pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")
metadata <-
  colnames(pf) %>% str_subset("^Metadata_")

profiles <- pf.sc[ , colSums(is.na(pf.sc)) == 0] # remove features (columns) where there are NA

feat.keep <- intersect(colnames(profiles), colnames(pf %>% select(one_of(variables)))) # features that are not either in profiles or in pf are removed

profiles <- profiles[, names(profiles) %in% feat.keep] # remove features

variables <- colnames(profiles) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")
start.time <- Sys.time()

# remove features that are too correlated (threshold of .9)
profiles %<>%
  cytominer::select(
    sample = profiles,
    variables = variables,
    operation = "correlation_threshold"
  )
variables <-
  names(profiles) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken 

variables.fs.sc.fc <- variables
save(variables.fs.sc.fc, file="../../input/repurposing/single_cells/findCorrelation/feat_names/selected_features_findCorr.Rda")
