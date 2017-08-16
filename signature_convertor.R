# all usefull libraries
library(magrittr)
library(dplyr)

#' Function to do signature of the compounds (by averaging the profiles over the replicate)
#' @param pf the data file
#' @return profiles

profile_generator <- function(pf){ 

  # Metadata
  metadata <-
    names(pf) %>% str_subset("^Metadata_")
  
  # Variables
  variables <-
    names(pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")

  # average along replicate (keep only id and variables)
  signature <- 
    pf %>% 
    select(one_of(variables, 'Metadata_broad_sample')) %>%
    group_by(Metadata_broad_sample) %>% 
    summarise_each(funs(mean(., na.rm=TRUE)))
  
  # add metadata information
  signature %<>%
    left_join(., pf %>% select(one_of(metadata)), by = 'Metadata_broad_sample') %>%
    group_by(Metadata_broad_sample) %>%
    slice(1)
  
  return(signature)
}