# all usefull libraries
library(magrittr)
library(dplyr)
library(ggplot2)
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
#' @param filename name of the data file in '.rds'
#' @param n.feat number of features in the sets
#' @param feat.selected boolean if data as been feature selecteed
#' @param N number of data to make the non replicate distance distribution
#' @param seed seed number for reproductibility
#' @param nCPU number of CPU cores for parallelization
#' @return hit ratio
  
hit_selection_jaccard <- function(filename, n.feat = 50, feat.selected = FALSE, N = 5000, seed = 42, nCPU = 7){ 
  message(paste('runing Jaccard Hit Selection for file: ',filename, 'with number of features per set = ', n.feat))
  start.time <- Sys.time()
  
  # number of CPU cores for parallelization
  registerDoMC(nCPU)
  
  # seed for the reproducibility
  set.seed(seed)
  
  ## Import data
  pf <- readRDS(file.path("..", "..", "input", "BBBC022_2013", "old", filename)) # 7680x803
  
  # Remove the negative control from the data
  pf$data <- filter(pf$data, !Image_Metadata_BROAD_ID %in% "") # 6400x803
  
  ## Separation of data
  # find the different compounds
  IDs <- distinct(pf$data, Image_Metadata_BROAD_ID)
  
  ## Distance of the data
  
  # function to calculate the distance between two sets
  distJaccard <- function(x, y){
    # 4 sets: n.feat first and n.feat lasts features
    x.sort <- sort(x)
    y.sort <- sort(y)
    
    A <- names(x.sort[1:n.feat])
    B <- names(y.sort[1:n.feat])
    
    # sorting each time
    C <- names(x.sort[seq.int(to = length(x.sort), length.out = n.feat)])
    D <- names(y.sort[seq.int(to = length(y.sort), length.out = n.feat)])
    
    d.top <- 
      (length(union(A, B)) - length(intersect(A, B)))/length(union(A, B))
    d.bottom <-
      (length(union(C, D)) - length(intersect(C, D)))/length(union(C, D))
    
    return((d.top + d.bottom)/2)
  }
  
  # loop over all IDs and save the median of the distance
  comp.dist.median <- foreach(i = 1:length(IDs$Image_Metadata_BROAD_ID), .combine=cbind) %dopar% {
    
    #filtering to choose only for one compound
    comp <-
      filter(pf$data, Image_Metadata_BROAD_ID %in% IDs$Image_Metadata_BROAD_ID[i])
    comp <-
      comp[, pf$feat_cols] %>%
      as.matrix()
    
    # distance of the features
    comp.dist <-
      lapply(1:dim(comp)[1], function(x) (lapply(1:dim(comp)[1], function(y) return(distJaccard(comp[x,], comp[y,]))) %>%
                                            unlist)) %>%
      do.call(rbind, .)
    
    # save median of the correlation
    comp.dist.median <- 
      median(comp.dist[lower.tri(comp.dist)],na.rm=TRUE)
    
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
      pf$data %>% 
      group_by(Image_Metadata_BROAD_ID) %>% 
      sample_n(1, replace = FALSE) %>% 
      ungroup(random.replicate)
    random.replicate <- sample_n(random.replicate, 4, replace = FALSE)
    
    comp <- random.replicate[,pf$feat_cols] %>% 
      as.matrix()
    
    # distance of the features
    comp.dist <- 
      lapply(1:dim(comp)[1], function(x) (lapply(1:dim(comp)[1], function(y) return(distJaccard(comp[x,], comp[y,]))) %>%
                                            unlist)) %>% 
      do.call(rbind, .)
    
    # median of the non replicate distance
    random.replicate.dist.median <- median(comp.dist[lower.tri(comp.dist)],na.rm=TRUE)
    
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
  hit.select.IDs <- IDs$Image_Metadata_BROAD_ID[inds]
  
  # ratio of strong median replicate correlation
  hit.ratio <- length(hit.select)/length(comp.dist.median)
  
  message(paste('Hit ratio: ', hit.ratio))
  
  ## Saving data
    
  # select high median correlation replicate 
  pf$data %<>% 
    filter(Image_Metadata_BROAD_ID %in% hit.select.IDs)
  
  # save new dataset
  if(feat.selected){
    filename.save <- paste("../../input/BBBC022_2013/old/Hit_jaccard_", 
                           toString(n.feat), 
                           "n_fs_", 
                           toString(round(hit.ratio*10000)), 
                           ".rds", 
                           sep = "")
  } else {
    filename.save <- paste("../../input/BBBC022_2013/old/Hit_jaccard_", 
                           toString(n.feat), 
                           "n_", 
                           toString(round(hit.ratio*10000)), 
                           ".rds", 
                           sep = "")
  }
  
  #### uncomment if want to save file
  #pf %>%
  #  saveRDS(filename.save)
  
  end.time <- Sys.time() # 1.4 mins (without feature selection)
  time.taken <- end.time - start.time
  time.taken

  message(paste('time to run function agglomerative clustering: ', time.taken))
  
  return(hit.ratio)
}