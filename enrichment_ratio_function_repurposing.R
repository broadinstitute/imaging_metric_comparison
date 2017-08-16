library(magrittr)
library(dplyr)
library(foreach)
library(doMC)
library(reshape2)
library(stringr)

#' Function that perform the Enrichment ratio
#' This method is based on the MOA information, here separation with a vertical bar |
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
      left_join(.,
                n.moa.cmpd.pair,
                by = c("cmpd1" = "Var1", "cmpd2" = "Var2"))
    
    top.moa.matching <- 
      cor.cmpd.pair %>%
      group_by(cmpd1) %>% # group by compound
      arrange(cmpd1, desc(corr)) %>% # sort correlation from higher to lower
      filter(corr > quantile(corr, 1.0-top.x)) %>% # look at the top x% correlation in each group
      summarise(p = sum(value)/n()) %>% # percentage of similar moa
      ungroup()
    
    final.number <-
      top.moa.matching %>%
      filter(p > 0) %>% # p bigger than 0 mean that at least there is one moa in common
      summarise(n = n()/nrow(top.moa.matching)) %>% # number of compound that have a least one MOA in common divided by the total number of compounds
      as.numeric()
    
    return(final.number)
  }

  # MOAs data
  all_metadata <- 
    pf %>%
    select(one_of(metadata.pf))
  
  all_metadata %<>%
    group_by(Metadata_broad_sample) %>%
    slice(1) %>%
    ungroup
  
  for (i in 1:nrow(all_metadata)){
    # if there are more than 1 moa associated
    if (str_detect(all_metadata$Metadata_moa[i], "\\|")){
      t1 <- str_trim(str_split(all_metadata$Metadata_moa[i], "\\|")[[1]])
      all_metadata$Metadata_moa[i] <- t1[1]
      for(j in 2:length(t1)){
        new.row <- all_metadata[i,]
        new.row$Metadata_moa <- t1[j]
        all_metadata <- rbind(all_metadata, new.row)
      }
    }
  }
  
  # select MOA that are appearing more than once (meaning at least two compounds are related to it)
  n.MOA <- table(all_metadata$Metadata_moa) %>% as.data.frame() %>% filter(Freq != 1)
  all_metadata %<>%  filter(Metadata_moa %in% n.MOA$Var1)
  
  # average along replicate (keep only id and variables)
  pf.cmpds <- 
    pf %>% 
    filter(Metadata_broad_sample %in% all_metadata$Metadata_broad_sample) %>% # select ID that have a unique compound
    select(one_of(variables, 'Metadata_broad_sample')) %>%
    group_by(Metadata_broad_sample) %>% 
    summarise_each(funs(mean(., na.rm=TRUE))) %>%
    as.data.frame()
  
  # keep track of ID of the compounds
  row.names(pf.cmpds) <- pf.cmpds$Metadata_broad_sample
  
  # Attention: since some compounds have more than one MOAs, few rows have same compounds name.
  pf.cmpd.meta <- 
    all_metadata %>%
    dplyr::left_join(., pf.cmpds, by = "Metadata_broad_sample") %>%
    as.data.frame
  
  # binary indicator matrix of ID vs MOA 
  n.moa.ID <- 
    pf.cmpd.meta %>%
    select(Metadata_moa, Metadata_broad_sample) %>%
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