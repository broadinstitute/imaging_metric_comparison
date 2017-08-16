library(magrittr)
library(dplyr)
library(tidyverse)
library(stringr)

# feature selection using SVD-entropy (SR) on single cells on the Repurposing dataset

# import single cell DMSO data
pf.sc <- 
  readRDS(file.path("../../input/repurposing/dmso_single_cell/Repurposing_dmso_sc.rds")) %>% 
  scale() %>% 
  as.data.frame() #20,491x1777 -> no: 20,419x815

# import profile data
pf <- readRDS(file.path("../../input/repurposing/repurposing_normalized.rds"))

pf %<>% filter(Metadata_mmoles_per_liter == 10)

# there are two compounds that are removed because no MOA are associated and also DMSO
pf %<>% filter(!is.na(Metadata_moa))

# remove "BRD-K60230970-001-10-0@19.9991253536431" because it is a toxic compounds
pf %<>% filter(Metadata_broad_sample != "BRD-K60230970-001-10-0@19.9991253536431")

# remove features with NA
pf <- pf[ , colSums(is.na(pf)) == 0] # from 1784 to 1615

variables <-
  colnames(pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")
metadata <-
  colnames(pf) %>% str_subset("^Metadata_")


profiles <- pf.sc[ , colSums(is.na(pf.sc)) == 0] # remove features (columns) where there are NA

feat.keep <- intersect(colnames(profiles), colnames(pf %>% select(one_of(variables)))) # features that are not either in profiles or in pf are removed

profiles <- profiles[, names(profiles) %in% feat.keep] # remove features

# transpose the dataset to have featxobs (mxn)
X <- profiles %>% as.matrix() %>% t(.)
A <- tcrossprod(X, X) 

print("cross product ok")

# load the c++ function
Rcpp::sourceCpp('ranking_SVD_entropy.cpp')

print("cpp function ok")

# different type of feature selection using SVD-entropy (SR, FS1, FS2), FS2 we must select the number of features that we want.
feat.idx <- CE_entropy_SR(A) #CE_entropy_FS2_new(A, 800)

# save entropy
save(feat.idx, file = "entropy_feat_SR_repurposing_sc.RData")
print("entropy saved")

# SR is giving the entropy -> should link entropy to features name and sort them!
# should select only entropy that are positive

##############################
# in order to save all sets of 100 features.
n.feats <- seq(from = 100, to = 1600, by = 50)

for(n in n.feats){
  # names of the features to keep
  names.CEi <- rownames(X)[feat.idx[1:n]]
  
  names.CEi %<>% as.data.frame() 
  names.CEi <- rename(names.CEi, feat.name = .)
  
  # save features names (to change the path and the name)
  save(names.CEi, file=paste("../../input/repurposing/feat_names/SVD_cells_", as.character(n), "feat.Rda", sep=""))
}