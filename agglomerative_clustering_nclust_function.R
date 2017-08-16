# all usefull libraries
library(magrittr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(foreach)
library(doMC)
library(stringr)
library(tidyverse)

# This function is doing the following for BBBC022 dataset:
# 1) average along the replicate to have one single signature per compound
# 2) Agglomerative clustering using average linkage function
# 3) separate in numClust clusters
# 4) use Fisher's exact test to calculate the Odds ratio: same/different MOAs and same/different clusters classes.

# pf: profile of BBBC022 after hit selection
# numClust: number of clusters
agglomarative_clustering <- function(pf, numClust){
  message(paste('runing agglomerative clustering with k = ', numClust))
  start.time <- Sys.time()
  
  # find the different compounds
  IDs <- distinct(pf$data, Image_Metadata_BROAD_ID)
  #dim(IDs)
  
  # mapping IDs-Compounds: give access to the compound knowing the IDs
  map.names <-
    pf$data %>%
    dplyr::select(one_of(c("Image_Metadata_BROAD_ID","Image_Metadata_SOURCE_COMPOUND_NAME"))) %>%
    unique %>%
    mutate_if(is.factor, as.character) %>%
    as.data.frame()
  rownames(map.names) <- map.names$Image_Metadata_BROAD_ID
  
  
  # average along replicate (keep only id and variables)
  sum.comp <- 
    pf$data %>% 
    group_by(Image_Metadata_BROAD_ID) %>% 
    summarise_each(funs(mean(., na.rm=TRUE)), -c(Image_Metadata_SOURCE_COMPOUND_NAME, Well, Plate)) %>%
    mutate_if(is.factor, as.character)
  
  # add column of compound and convert tu uppercase (because some compounds same name but in lower/upper case)
  sum.comp %<>%
    dplyr::left_join(., map.names, by = "Image_Metadata_BROAD_ID") %>% 
    mutate_if(is.character, funs(toupper))
  #dim(sum.comp)
  
  # keep only one componment (first one)
  sum.comp %<>%
    group_by(Image_Metadata_SOURCE_COMPOUND_NAME) %>%
    slice(1) %>%
    ungroup
  #dim(sum.comp)
  
  # correlation compound-compound
  cor.comp <-
    sum.comp[, pf$feat_cols] %>% 
    as.matrix() %>%
    t() %>%
    cor()
  
  # metric for clustering: 1 - correlation
  met.clust <- as.dist(1 - cor.comp)
  
  # Agglomerative clustering using a linkage function
  clust.res <- hclust(met.clust, method = "average")
  
  #plot(clust.res, cex = 0.01, hang = -1)
  #rect.hclust(clust.res, k = numClust, border = "red")
  
  
  # find the cluster for each compound
  sum.comp.clust <- 
    sum.comp %>%
    mutate(cluster = cutree(clust.res, k = numClust))
  
  
  # import MOAs data
  moa <- 
    read.csv("../../input/BBBC022_2013/MOAs.csv", na.strings = c("", "NA")) %>%
    mutate_if(is.factor, as.character) %>%
    plyr::rename(c("Name" = "Image_Metadata_SOURCE_COMPOUND_NAME")) 
  # compounds name to upper case
  moa$Image_Metadata_SOURCE_COMPOUND_NAME <-
    lapply(moa[, "Image_Metadata_SOURCE_COMPOUND_NAME"], stringr::str_to_upper) %>%
    unlist
  moa %<>% 
    group_by(Image_Metadata_SOURCE_COMPOUND_NAME) %>%
    slice(1) %>%
    ungroup
  
  # join moa data to cluster data 
  final.sum.comp.clust <-
    sum.comp.clust %>%
    select(Image_Metadata_SOURCE_COMPOUND_NAME, cluster) %>%
    dplyr::left_join(., moa, by = "Image_Metadata_SOURCE_COMPOUND_NAME") 
  #dim(final.sum.comp.clust)
  # remove row where MOA is NA
  final.sum.comp.clust.moa <-
    final.sum.comp.clust[complete.cases(final.sum.comp.clust[,"MOA"]),]
  #dim(final.sum.comp.clust.moa)
  
  # select only the MOA and the cluster column
  df.clust.moa <-
    final.sum.comp.clust.moa %>% 
    select(cluster, MOA)
  #dim(df.clust.moa)
  
  for (i in 1:dim(df.clust.moa)[1]){
    # if there are more than 1 moa associated
    if (str_detect(df.clust.moa$MOA[i], ",")){
      t1 <- str_trim(str_split(df.clust.moa$MOA[i], ",")[[1]])
      df.clust.moa$MOA[i] <- t1[1]
      new.row <- df.clust.moa[i,]
      new.row$MOA <- t1[2]
      df.clust.moa <- rbind(df.clust.moa, new.row)
    }
    
  }
  
  dim.before <- dim(df.clust.moa)[1]
  
  # remove compound that are clustered on their own!
  df.clust.moa %<>%
    group_by(cluster) %>% 
    mutate(n = n()) %>% # count number of element in each cluster
    ungroup(cluster) %>%
    filter(n > 1) # remove element that have only themself in their cluster
  
  dim.after <- dim(df.clust.moa)[1]
  num.singleton <- dim.before - dim.after
  threshold <- .05*dim(IDs)[1]
  
  if(num.singleton > threshold){
    message(paste('number of singletons: ', num.singleton, ' threshold: ', threshold))
    stop('number of singleton too big! ')
  }
  
  
  # remove row where MOA is NA
  moa.ext <-
    moa[complete.cases(moa[,"MOA"]),]
  #dim(moa.ext)
  
  # duplicate row where multiple moas for same compound
  for (i in 1:dim(moa.ext)[1]){
    # if there are more than 1 moa associated
    if (str_detect(moa.ext$MOA[i], ",")){
      t1 <- str_trim(str_split(moa.ext$MOA[i], ",")[[1]])
      moa.ext$MOA[i] <- t1[1]
      new.row <- moa.ext[i,]
      new.row$MOA <- t1[2]
      moa.ext <- rbind(moa.ext, new.row)
    }
    
  }
  #dim(moa.ext)
  
  # create matrix of comparison
  moa.mat <- outer(df.clust.moa$MOA, df.clust.moa$MOA, function(x, y) x==y)
  clust.mat <- outer(df.clust.moa$cluster, df.clust.moa$cluster, function(x, y) x==y)
  # transform into a vector with the indices and combine both moa and cluster information together, keep only the upper triangle of the matrix 
  tmp1 <- melt(moa.mat) %>%
    plyr::rename(c("value" = "same.moa"))
  tmp2 <- melt(clust.mat)
  tmp1$same.clust <- tmp2$value
  tmp1 %<>% filter(Var1 < Var2)
  
  # contingency table
  contingency.table <- 
    tmp1 %>%
    group_by(same.moa, same.clust) %>%
    summarise(cnt = n()) %>%
    xtabs(cnt ~ same.moa+same.clust, data = .)
  contingency.mat <-
    matrix(c(contingency.table[2,2], contingency.table[2,1], contingency.table[1,2], contingency.table[1,1]), 
           nrow = 2, 
           ncol = 2, 
           byrow = TRUE)
  print(contingency.table)
  
  # fisher test
  f.test.res <- fisher.test(contingency.mat, alternative = "greater")
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste('time to run function agglomerative clustering: ',time.taken))

  return(f.test.res)
}