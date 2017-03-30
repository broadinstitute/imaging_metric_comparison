
library(magrittr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(foreach)
library(doMC)
library(stringr)
library(tidyverse)

#' Alternative method to agglomerative clustering. Look at the top top.corr percentage compound pair correlation.
#' Perform a fisher exact's test for top top.corr % and the rest versus same MOAs and different MOAs.
#'
#' @param filename name of the data file in '.rds'
#' @param top.corr top percentage of the correlation
#' @return result of the Fisher's exact test


alternative_method <- function(filename, top.corr = 5){
  
  message(paste('runing alternative method for file: ',filename, ' with top ', top.corr, '% correlation.'))
  start.time <- Sys.time()
  
  # laoding data
  pf <- readRDS(file.path("..", "..", "input", "BBBC022_2013", "old", filename))

  profiles <- pf$data
  
  variables <- pf$feat_cols
  
  ## Summary of the compound
  # find the different compounds
  IDs <- distinct(profiles, Image_Metadata_BROAD_ID)
  
  # mapping IDs-Compounds: give access to the compound knowing the IDs
  map.names <-
    profiles %>%
    dplyr::select(one_of(c("Image_Metadata_BROAD_ID","Image_Metadata_SOURCE_COMPOUND_NAME"))) %>% # select ID and Compound
    unique %>% # keep unique rows, one column can have same value
    mutate_if(is.factor, as.character) %>% # if is a factor transform into a character
    as.data.frame() # be sure it is a dataframe
  rownames(map.names) <- map.names$Image_Metadata_BROAD_ID # transform the names of the rows into the IDs
  
  # average along replicate (keep only id and variables)
  sum.comp <- 
    profiles %>% 
    group_by(Image_Metadata_BROAD_ID) %>% 
    summarise_each(funs(mean(., na.rm=TRUE)), -c(Image_Metadata_SOURCE_COMPOUND_NAME, Well, Plate)) %>% # mean along replicate
    mutate_if(is.factor, as.character)
  
  # add column of compound and convert to uppercase (because some compounds same name but in lower/upper case)
  sum.comp %<>%
    dplyr::left_join(., map.names, by = "Image_Metadata_BROAD_ID") %>% 
    mutate_if(is.character, funs(toupper)) # same compound can be writen in upper and lower case and looks different
  
  # keep only one componment (first one) (this is a safe way to perform)
  sum.comp %<>%
    group_by(Image_Metadata_SOURCE_COMPOUND_NAME) %>%
    slice(1) %>%
    ungroup %>%
    as.data.frame()
  
  row.names(sum.comp) <- sum.comp$Image_Metadata_BROAD_ID
  
  # correlation compound-compound
  cor.comp <-
    sum.comp[, variables] %>% 
    as.matrix() %>%
    t() %>%
    cor()
  
  # extract upper part of the compound compound correlation matrix and sort
  ind <- which( upper.tri(cor.comp) , arr.ind = TRUE ) # indices of the upper triangle of the correlation matrix
  
  cor.compound.pair <- 
    data.frame(col = dimnames(cor.comp)[[2]][ind[,2]],
               row = dimnames(cor.comp)[[1]][ind[,1]],
               val = cor.comp[ind]) %>%
    dplyr::arrange(desc(val))
  
  # threshold of the top 5%
  quant <- quantile(cor.compound.pair$val, 1.0 - top.corr/100)
  
  ## MOA data
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
  
  # duplicate rows when multiple MOAs for one sample
  for (i in 1:dim(moa)[1]){
    # if there are more than 1 moa associated
    if (!is.na(moa$MOA[i]) & str_detect(moa$MOA[i], ",")){
      t1 <- str_trim(str_split(moa$MOA[i], ",")[[1]])
      moa$MOA[i] <- t1[1]
      new.row <- moa[i,]
      new.row$MOA <- t1[2]
      moa <- rbind(moa, new.row)
    }
    
  }
  # joining moa and compounds information
  tmp <- moa %>%
    dplyr::left_join(., sum.comp, by = "Image_Metadata_SOURCE_COMPOUND_NAME")
  dim(tmp)
  # binary indicator matrix of ID vs MOA 
  comp_vs_moa <- tmp %>%
    select(MOA, Image_Metadata_BROAD_ID) 
  comp_vs_moa <-  as.matrix(as.data.frame.matrix(table(comp_vs_moa)))
  
  # number of moa in common for each ID pairs
  comp_vs_comp <-
    t(comp_vs_moa)  %*% comp_vs_moa %>% # calculate matrix compound ID - compound ID relation
    melt(.) %>% # transform matrix into column
    filter(value > 0) %>% # filter value bigger than 0
    filter(Var1 != Var2)
  
  cor.compound.pair %<>% 
    mutate_if(is.factor, as.character) %>%
    dplyr::left_join(., 
                     mutate_if(comp_vs_comp, is.factor, as.character), 
                     by = c("col" = "Var1", "row" = "Var2"))
  
  # replace NA with 0
  cor.compound.pair[is.na(cor.compound.pair)] <- 0
  
  # top 5% correlated compound, number of same moa for all pair
  top5.same.moa <- 
    cor.compound.pair %>% 
    filter(val >= quant) %>% 
    select(value) %>% 
    sum()
  
  # top 5% correlated compound, number of pairs not sharing same MOA
  top5.diff.moa <-
    cor.compound.pair %>% 
    filter(val >= quant) %>% 
    filter(value == 0) %>% 
    nrow()
  
  # 95% rest, number of same moa for all pair
  rest95.same.moa <- 
    cor.compound.pair %>%
    filter(val < quant) %>%
    select(value) %>%
    sum()
  
  # 95% rest, number of pairs not sharing same MOA
  rest95.diff.moa <- 
    cor.compound.pair %>%
    filter(val < quant) %>%
    filter(value == 0) %>%
    nrow()
  
  contingency.mat <-
    matrix(c(top5.same.moa, 
             top5.diff.moa, 
             rest95.same.moa, 
             rest95.diff.moa), 
           nrow = 2, 
           ncol = 2, 
           byrow = TRUE)
  
  rownames(contingency.mat) <- c("top 5%", "rest")
  colnames(contingency.mat) <- c("same moa", "diff moa")
  
  contingency.mat
  print(contingency.mat)
  
  # fisher's exact test
  f.test.res <- fisher.test(contingency.mat, alternative = "greater")
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste('time to run function: ',time.taken))
  
  return(f.test.res)
}