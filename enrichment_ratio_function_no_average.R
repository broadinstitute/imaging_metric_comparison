library(magrittr)
library(dplyr)
library(foreach)
library(doMC)
library(reshape2)
library(stringr)

#' Function that perform enrichment ratio without averaging the profiles
#'
#' @param pf the data file
#' @param top.x top percentage of matching compound
#' @param seed seed number for reproductibility
#' @param nCPU number of CPU cores for parallelization
#' @param N number of data to make the non replicate distance distribution
#' @param filename name of the dataframe
#' @param method name of the method used for hit selection
#' @return enrichment ratio

enrichment_ratio <- function(pf,
                             top.x = 0.02,
                             seed = 42,
                             nCPU = 7,
                             N = 1000,
                             filename,
                             method = "Pearson"){ 
  
  # for reproductibility
  set.seed(seed)
  
  # number of CPU cores for parallelization
  registerDoMC(nCPU)
  
  variables <- names(pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")
  # Metadata
  metadata.pf <- names(pf) %>% str_subset("^Metadata_")
  
  top.percentage.matching.moa <- function(cor.cmpd, n.moa.cmpd.pair){
    cor.cmpd.pair <- 
      melt(cor.cmpd) %>%
      rename(cmpd1 = Var1, cmpd2 = Var2, corr = value) %>% # rename columns
      filter(cmpd1 != cmpd2) %>% # remove column were same compound
      mutate(cmpd1_n = lapply(str_split(cmpd1, "_"), function(x) x[1]) %>% unlist) %>%
      mutate(cmpd2_n = lapply(str_split(cmpd2, "_"), function(x) x[1]) %>% unlist) %>%
      select(one_of(c("corr", "cmpd1_n", "cmpd2_n"))) %>%
      group_by(cmpd1_n, cmpd2_n) %>%
      summarise_each(funs(median(., na.rm=TRUE))) %>% # do average of all correlation TODO: maybe median also?
      left_join(.,
                n.moa.cmpd.pair,
                by = c("cmpd1_n" = "Var1", "cmpd2_n" = "Var2")) %>%
      filter(cmpd1_n != cmpd2_n) # remove NA
    
    top.moa.matching <- 
      cor.cmpd.pair %>%
      group_by(cmpd1_n) %>% # group by compound
      arrange(cmpd1_n, desc(corr)) %>% # sort correlation from higher to lower
      filter(corr > quantile(corr, 1.0-top.x)) %>% # look at the top 10% correlation in each group
      summarise(p = sum(value)/n()) %>% # percentage of similar moa
      ungroup()
    
    final.number <-
      top.moa.matching %>%
      filter(p > 0) %>% # p bigger than 0 mean that at least there is one moa in common
      summarise(n = n()/nrow(top.moa.matching)) %>% # number of compound that have a least one MOA in common divided by the total number of compounds
      as.numeric()
    
    return(final.number)
  }
  ## Metadata 
  # import MOAs data
  moa <- 
    read.csv("../../input/MOAs.csv", na.strings = c("", "NA")) %>%
    mutate_if(is.factor, as.character) # has to do that to duplicate rows where multiple MOA
  
  moa %<>% 
    filter(!is.na(MOA_group)) %>% # remove rows where moa = NA
    rename(Image_Metadata_SOURCE_COMPOUND_NAME = Name) # rename compound name
  
  # remove duplicate compounds
  moa$Image_Metadata_SOURCE_COMPOUND_NAME <- 
    toupper(moa$Image_Metadata_SOURCE_COMPOUND_NAME) %>% # put compound names in upper character (to make it comparable)
    str_replace_all(fixed(" "), "") # remove white space
  
  moa %<>%
    group_by(Image_Metadata_SOURCE_COMPOUND_NAME) %>%
    slice(1) %>%
    ungroup
  
  # duplicate rows where multiple MOA associated to one compounds
  for (i in 1:nrow(moa)){
    # if there are more than 1 moa associated
    if (str_detect(moa$MOA_group[i], ",")){
      t1 <- str_trim(str_split(moa$MOA_group[i], ",")[[1]])
      moa$MOA_group[i] <- t1[1]
      new.row <- moa[i,]
      new.row$MOA_group <- t1[2]
      moa <- rbind(moa, new.row)
    }
  }
  
  pf$Metadata_cpd_name %<>% 
    as.character() %>% 
    str_replace_all(fixed(" "), "") # remove white space
  
  # mapping IDs-Compounds: give access to the compound knowing the IDs
  id.cmpds <-
    pf %>%
    dplyr::select(one_of(c("Metadata_broad_sample","Metadata_cpd_name"))) %>% # select ID and Compound
    mutate_if(is.character, funs(toupper)) %>% # same compound can be writen in upper and lower case and looks different
    group_by(Metadata_cpd_name) %>%
    slice(1) %>%
    ungroup
  
  # metadata: compound, MOA, target, ID
  metadata <- 
    moa %>%
    left_join(., id.cmpds, by = c("Image_Metadata_SOURCE_COMPOUND_NAME"="Metadata_cpd_name"))
  
  # select only rows that have a BROAD_ID
  metadata %<>% filter(!is.na(Metadata_broad_sample))
  
  # select MOA that are appearing more than once (meaning at least two compounds are related to it)
  n.MOA <- table(metadata$MOA_group) %>% as.data.frame() %>% filter(Freq != 1)
  metadata %<>%  filter(MOA_group %in% n.MOA$Var1)
  
  pf.cmpds <- 
    pf %>% 
    filter(Metadata_broad_sample %in% metadata$Metadata_broad_sample) %>% # select ID that have a unique compound
    select(one_of(variables, 'Metadata_broad_sample')) %>%
    group_by(Metadata_broad_sample) %>% 
    mutate(rep = row_number()) %>%
    mutate(Metadata_broad_sample_rep = paste(Metadata_broad_sample, rep, sep="_")) %>% # to have a unique name broad sample ID
    #summarise_each(funs(mean(., na.rm=TRUE))) %>%
    as.data.frame()
  

  # keep track of ID of the compounds
  row.names(pf.cmpds) <- pf.cmpds$Metadata_broad_sample_rep
  
  pf.cmpd.meta <- 
    metadata %>%
    dplyr::left_join(., pf.cmpds, by = "Metadata_broad_sample") %>%
    as.data.frame
  
  # binary indicator matrix of ID vs MOA 
  n.moa.ID <- 
    pf.cmpd.meta %>%
    select(MOA_group, Metadata_broad_sample) %>%
    group_by(Metadata_broad_sample) %>%
    slice(1) %>% #### remove duplicate
    table %>%
    as.data.frame.matrix %>%
    as.matrix
  
  # number of moa in common for each ID pairs
  n.moa.cmpd.pair <-
    t(n.moa.ID)  %*% n.moa.ID %>% # calculate matrix compound ID - compound ID relation
    melt(.) %>% # transform matrix into column
    filter(Var1 != Var2)
  
  
  ## Correlation compound-compound 
  
  # correlation compound-compound
  cor.cmpd <-
    pf.cmpds[, variables] %>% 
    as.matrix() %>%
    t() %>%
    cor()
  
  ## Percentage
  
  final.number <- top.percentage.matching.moa(cor.cmpd, n.moa.cmpd.pair)
  
  
  ## Baseline
  
  v <- rownames(cor.cmpd) # extract names of the rows and columns to shuffle
  
  set.seed(seed)
  seeds <- sample(1:10000, N, replace=F)
  
  
  random.percent <- foreach(i = 1:N, .combine=cbind) %dopar% {
    # for reproducibility
    set.seed(seeds[i])
    # randomly shuffle names of the compounds
    t <- sample(v)
    
    # shuflle in the same random way the names of the rows and the columns
    cor.comp.random <- cor.cmpd
    rownames(cor.comp.random) <- t
    colnames(cor.comp.random) <- t
    
    # knn for random
    random.percent <- top.percentage.matching.moa(cor.comp.random, n.moa.cmpd.pair)
  }
  
  ## Results
  enrichment.ratio <- data.frame(mean = mean(random.percent), 
                                 quant = quantile(random.percent, .95), 
                                 percent = final.number, 
                                 filename = filename,
                                 method = method)
  return(enrichment.ratio)
}