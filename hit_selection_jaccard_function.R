# all usefull libraries
library(magrittr)
library(dplyr)
library(foreach)
library(doMC)
library(stringr)
library(tidyverse)


#' Function that perform Hit Selection.
#' Hit selection is performed based on the distance defined following the Jaccard distance, which calculate the dissimilarires between sample sets.
#' The distance is defined as the mean of the distance between 4 sets.
#' A) n.feat top features of component x
#' B) n.feat top features of component y
#' C) n.feat bottom features of component x
#' D) n.feat bottom features of component y
#' $dist(x,y) = \frac{dist_J(A,B) + dist_J(C,D)}{2}$
#' $dist_J(i,j) = \frac{|A \cup B| - |A \cap B|}{|A \cup B|}$ 
#'
#' @param df dataframe in '.rds'
#' @param filename name of the data file
#' @param n.feat number of features in the sets
#' @param feat.selected boolean if data as been feature selecteed
#' @param n.replicate number of replicate for the basline (random comparison)
#' @param N number of data to make the non replicate distance distribution
#' @param seed seed number for reproductibility
#' @param nCPU number of CPU cores for parallelization
#' @return hit ratio

hit_selection_jaccard <- function(pf, 
                                  filename = "Hit_jaccard", 
                                  n.feat = 50, 
                                  feat.selected = FALSE, 
                                  n.replicate = 4,
                                  N = 5000, 
                                  seed = 42,
                                  nCPU = 7,
                                  dir.save = "BBBC022_2013/selected_single_cells",
                                  dir.save.plus = "/hit_selected/random/Jaccard/"){ 

  message(paste('Running Jaccard Hit selection...', filename))
  start.time <- Sys.time()
  
  # number of CPU cores for parallelization
  registerDoMC(nCPU)
  
  # seed for the reproducibility
  set.seed(seed)
  
  # Variables
  variables <-
    names(pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")
  
  ## Separation of data
  # find the different compounds
  IDs <- unique(pf$Metadata_broad_sample)
  
  ## Distance of the data
  
  # load the c++ function
  Rcpp::sourceCpp('jaccard_distance_function.cpp',rebuild = FALSE)
  
  # loop over all IDs and save the median of the distance
  comp.dist.median <- foreach(i = 1:length(IDs), .combine=cbind) %dopar% {
    
    #filtering to choose only for one compound
    comp <-
      filter(pf, Metadata_broad_sample %in% IDs[i])
    comp %<>% 
      select(one_of(variables)) %>%
      as.matrix()
    
    # distance of the features
    comp.dist <- vecJaccardDistance(comp, n.feat)
    
    # save median of the distance
    comp.dist.median <- 
      median(comp.dist, na.rm=TRUE)
    
  }
  
  # hist(comp.dist.median,
  #     main="Histogram for Median Replicate Distance",
  #     xlab="Median Replicate Distance", 
  #     xlim = range(0, 1))
  
  ## Thresholding of poor replicate distance
  
  # random sequence for reproducibility
  a <- sample(1:10000, N, replace=F)
  
  # loop over N times to get a distribution
  random.replicate.dist.median <- foreach(i = 1:N, .combine=cbind) %dopar% {
    # set seed according to random sequence
    set.seed(a[i])
    
    # group by IDs
    # sample fixed number per group -> choose 4 replicates randomly from different group
    random.replicate <- 
      pf %>% 
      group_by(Metadata_broad_sample) %>% 
      sample_n(1, replace = FALSE) %>% 
      ungroup(random.replicate)
    random.replicate <- sample_n(random.replicate, n.replicate, replace = FALSE)
    
    comp <- random.replicate[,variables] %>% 
      as.matrix()
    
    # distance of the features
    comp.dist <- vecJaccardDistance(comp, n.feat)
    
    # median of the non replicate distance
    random.replicate.dist.median <- median(comp.dist,na.rm=TRUE)
    
  }
  
  # histogram plot
  #hist(random.replicate.dist.median,
  #   main="Histogram for Non Replicate Median Distance",
  #   xlab="Non Replicate Median Distance", 
  #   xlim = range(0, 1))
  
  # threshold to determine if can reject H0
  thres <- quantile(random.replicate.dist.median, .05)
  
  ## Hit Selection
  
  # find indices of replicate median distance < threshold
  inds <- which(comp.dist.median < thres)
  
  # find values of the median that are hit selected
  hit.select <- comp.dist.median[inds]
  
  # find component that are hit selected
  hit.select.IDs <- IDs[inds]
  
  # ratio of strong median replicate correlation
  hit.ratio <- length(hit.select)/length(comp.dist.median)
  
  ## Saving data
  
  # select high median correlation replicate 
  pf %<>% 
    filter(Metadata_broad_sample %in% hit.select.IDs)
  
  # save new dataset
  if(feat.selected){
    filename.save <- paste("../../input/",
                           dir.save,
                           dir.save.plus, 
                           strsplit(filename, ".rds"),
                           "_",
                           toString(n.feat), 
                           "n_FS2_seed_",
                           seed,
                           ".rds", 
                           sep = "")
  } else {
    filename.save <- paste("../../input/",
                           dir.save,
                           dir.save.plus,
                           strsplit(filename, ".rds"),
                           "_",
                           toString(n.feat), 
                           "n_seed_",
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
  
  return(hit.ratio)
}