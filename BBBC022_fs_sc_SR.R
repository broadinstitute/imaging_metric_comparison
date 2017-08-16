library(magrittr)
library(dplyr)
library(tidyverse)
library(stringr)

# import data
pf.sc <- 
  readRDS(file.path("../../input/BBBC022_2013/single_cells/BBBC022_sampled_10000.rds")) %>% 
  #  select(one_of(variables.fs.sc.fc)) %>% 
  scale() %>% 
  as.data.frame() #20,491x1777 -> no: 20,419x815

pf <- readRDS(file.path("../../input/BBBC022_2013/old/Pf_Gustafsdottir.rds"))

# difference of features name between old and new version:
# Hoechst -> DNA, Ph_golgi -> AGP, Syto -> RNA
pf$feat_cols <- str_replace(pf$feat_cols, "Hoechst", "DNA")
pf$feat_cols <- str_replace(pf$feat_cols, "Ph_golgi", "AGP")
pf$feat_cols <- str_replace(pf$feat_cols, "Syto", "RNA")
cl <- colnames(pf$data)
cl <- str_replace(cl, "Hoechst", "DNA")
cl <- str_replace(cl, "Ph_golgi", "AGP")
cl <- str_replace(cl, "Syto", "RNA")
colnames(pf$data) <- cl

profiles <- pf.sc[ , colSums(is.na(pf.sc)) == 0] # remove features (columns) where there are NA

feat.delete <- setdiff(colnames(profiles), pf$feat_cols) # features that are in profiles but not in pf

profiles <- profiles[, !names(profiles) %in% feat.delete] # remove features
start.time <- Sys.time()
# transpose the dataset to have featxobs (mxn)
X <- profiles %>% as.matrix() %>% t(.)
A <- tcrossprod(X, X) 

print("cross product ok")

# load the c++ function
Rcpp::sourceCpp('ranking_SVD_entropy.cpp')

print("cpp function ok")

# SR return the entropy for each feature.
feat.idx <- CE_entropy_SR(A)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken 

save(feat.idx, file = "../../input/BBBC022_2013/single_cells/SR/entropy_feat_SR_repurposing_sc.RData")
print("entropy saved")

# SR is giving the entropy -> should link entropy to features name and sort them!
# then with this features vector, can select the best top n and see whether it improves or not

##############################
x <- data.frame(feat.name = colnames(profiles), entropy = feat.idx)
x %<>% arrange(entropy)

## can add condition to select only entropy that is positive -> then no need to do the loop.

n.feats <- seq(from = 100, to = 800, by = 100)

for(n in n.feats){
  # names of the features to keep
  feat.names.CEi <- x$feat[1:n]
  
  save(feat.names.CEi, file=paste("../../input/BBBC022_2013/single_cells/SR/feat_selected/SVD_SR_", as.character(n), "feat.Rda", sep=""))
}
