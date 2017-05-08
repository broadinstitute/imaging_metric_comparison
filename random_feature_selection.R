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
  
  set.seed(seed = seed)
  random.features <- sample(1:length(pf$feat_cols), nb.feat, replace=F)
  
  tmp <- pf$data[,pf$feat_cols] %>% select(random.features)
  
  pf$data <- bind_cols(pf$data[,pf$factor_cols], tmp)
  
  pf$feat_cols <- pf$feat_cols[random.features]
  
  print("New dataset with random features selected!")
  
  return(pf)
}