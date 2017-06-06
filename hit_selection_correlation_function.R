# all usefull libraries
library(magrittr)
library(dplyr)
library(ggplot2)
library(foreach)
library(doMC)
library(stringr)
library(tidyverse)


#' Function that perform Hit Selection.
#' Hit selection is performed based on Pearson correlation.
#'
#' @param pf the data file
#' @param n.replicate number of replicate per compounds
#' @param filename name of the method
#' @param cor.method method of correlation: "pearson", "kendall", "spearman"
#' @param feat.selected boolean if data as been feature selecteed
#' @param N number of data to make the non replicate distance distribution
#' @param seed seed number for reproductibility
#' @param nCPU number of CPU cores for parallelization
#' @return hit ratio

hit_selection_correlation <- function(pf, 
                                      n.replicate = 4,
                                      filename = "Hit_pearson", 
                                      cor.method = "pearson", 
                                      feat.selected = FALSE, 
                                      N = 5000, 
                                      seed = 42, 
                                      nCPU = 7,
                                      dir.save = "BBBC022_2013/selected_single_cells"){ 
  message(paste('Running Pearson Hit selection...', filename))
  
  # computational time
  start.time <- Sys.time()
  
  # number of CPU cores for parallelization
  registerDoMC(nCPU)
  
  # seed for the reproducibility
  set.seed(seed)
  
  variables <-
    names(pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")
  
  ## Separation of data
  
  # find the different compounds
  IDs <- unique(pf$Metadata_broad_sample) #distinct(pf, Metadata_broad_sample )
  
  ## Correlation of the data
  comp.cor.median <- foreach(i = 1:length(IDs), .combine=cbind) %dopar% {
    
    #filtering to choose only for one compound
    comp <- filter(pf, Metadata_broad_sample %in% IDs[i])
    comp <-
      comp[, variables] %>%
      as.matrix()
    
    # correlation of the features
    comp.cor <- cor(t(comp), method = cor.method)
    
    # save median of the correlation
    comp.cor.median <- median(comp.cor[lower.tri(comp.cor)],na.rm=TRUE) # 1600 values
    
  }
  
  #hist(comp.cor.median,
  #     main="Histogram for Median Replicate Correlation",
  #     xlab="Median Replicate Correlation", xlim = range(-1, 1))
  
  ## Thresholding of poor replicate correlation
  # random sequence for reproducibility
  a <- sample(1:10000, N, replace=F)
  
  # loop over N times to get a distribution
  random.replicate.cor.median <- foreach(i = 1:N, .combine=cbind) %dopar% {
    # set seed according to random sequence
    set.seed(a[i])
    
    # group by IDs
    # sample fixed number per group -> choose 7 replicates randomly from different group
    random.replicate <- 
      pf %>% 
      group_by(Metadata_broad_sample) %>% 
      sample_n(1, replace = FALSE) %>% # choose randomly 1 sample for each group of IDs
      ungroup(random.replicate)
    # choose randomly n.replicate samples each coming from a different group
    random.replicate <- sample_n(random.replicate, n.replicate, replace = FALSE)
    
    # select only the features
    comp <- random.replicate[,variables] %>% 
      as.matrix()
    
    # correlation of the features using the method defined in the parameter (pearson, spearman, kendall)
    random.replicate.cor <- cor(t(comp), method = cor.method)
    
    # median of the non replicate correlation
    random.replicate.cor.median <- median(random.replicate.cor[lower.tri(random.replicate.cor)],na.rm=TRUE) # vector length N
    
  }
  
  # histogram plot
  #hist(random.replicate.cor.median,
  #   main="Histogram for Non Replicate Median Correlation",
  #   xlab="Non Replicate Median Correlation", xlim = range(-1, 1))
  
  
  # threshold to determine if can reject H0
  thres <- quantile(random.replicate.cor.median, .95)
  
  ## Hit Selection
  # find indices of replicate median correlation > threshold
  inds <- which(comp.cor.median > thres)
  
  # find values of the median that are hit selected
  hit.select <- comp.cor.median[inds]
  
  # find component that are hit selected
  hit.select.IDs <- IDs[inds]
  
  # ratio of strong median replicate correlation
  hit.ratio <- length(hit.select)/length(comp.cor.median)
  
  ## Saving data
  
  # select high median correlation replicate 
  pf %<>% 
    filter(Metadata_broad_sample %in% hit.select.IDs)
  
  # save new dataset
  if(feat.selected){
    filename.save <- paste("../../input/",
                           dir.save,
                           "/hit_selected/Pearson/", 
                           strsplit(filename, ".rds"),
                           "_FS2_seed",
                           seed,
                           ".rds", 
                           sep = "")
  } else {
    filename.save <- paste("../../input/",
                           dir.save,
                           "/hit_selected/Pearson/",
                           strsplit(filename, ".rds"),
                           "_seed_",
                           seed,
                           ".rds", 
                           sep = "")
  }
  
  #### uncomment if want to save file
  pf %>%
    saveRDS(filename.save)
  
  end.time <- Sys.time() 
  time.taken <- end.time - start.time
  time.taken
  
  #message(paste('time to run: ', time.taken))
  
  return(hit.ratio)
}