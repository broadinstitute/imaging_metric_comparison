# all usefull libraries
library(magrittr)
library(dplyr)

#' Function that select randomly n features
#' @param pf dataframe
#' @param nb.feat number of selected features
#' @param seed for reproducibility
#' @return dataframe with selected features

random_feature_selection <- function(pf, nb.feat, seed=42){
  
  print("Running random feature selection...")
  
  metadata <-
    colnames(pf) %>% str_subset("^Metadata_")
  
  variables <-
    colnames(pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")
  
  set.seed(seed = seed)
  random.features <- sample(1:length(variables), nb.feat, replace=F)
  
  tmp <- pf %>% select(one_of(variables)) %>% select(random.features)
  
  pf <- bind_cols(Pf %>% select(one_of(metadata)), tmp)
  
  print("New dataset with random features selected!")
  
  return(pf)
}