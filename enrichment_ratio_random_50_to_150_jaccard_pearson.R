library(magrittr)
library(dplyr)
library(foreach)
library(doMC)
library(reshape2)
library(stringr)
library(ggplot2)


## Parameters
# top percentage of matching compound.
top.x <- 0.02 # 0.1 = 10%, 0.05 = 5%, 0.02 = 2%

# for reproductibility
seed <- 42
set.seed(seed)

# number of CPU cores for parallelization
registerDoMC(7)

# dataframe of result
enrichment.ratio <- data.frame(mean = numeric(0), quant = numeric(0), percent = numeric(0), filename = character(0), method = character(0))

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


## Import data
data.filenames.p <- list.files(path = "../../input/BBBC022_2013/selected_single_cell_zoom/hit_selected/Pearson/")
data.filenames.p <- str_subset(data.filenames.p, "random")
data.filenames.j <- list.files(path = "../../input/BBBC022_2013/selected_single_cell_zoom/hit_selected/Jaccard/")
data.filenames.j <- str_subset(data.filenames.j, "random")

start.time <- Sys.time()

for(filename in data.filenames.p){
  print(filename)
  pf <- readRDS(file.path("..", "..", "input", "BBBC022_2013", "selected_single_cell_zoom", "hit_selected", "Pearson", filename))
  profiles <- pf$data
  variables <- pf$feat_cols
  
  ## Metadata 
  # import MOAs data
  moa <- 
    read.csv("../../input/BBBC022_2013/MOAs.csv", na.strings = c("", "NA")) %>%
    mutate_if(is.factor, as.character) # has to do that to duplicate rows where multiple MOA
  
  moa %<>% 
    filter(!is.na(MOA_group)) %>% # remove rows where moa = NA
    rename(Image_Metadata_SOURCE_COMPOUND_NAME = Name) # rename compound name
  
  # remove duplicate compounds
  moa$Image_Metadata_SOURCE_COMPOUND_NAME <- toupper(moa$Image_Metadata_SOURCE_COMPOUND_NAME) # put compound names in upper character (to make it comparable)
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
  
  # mapping IDs-Compounds: give access to the compound knowing the IDs
  id.cmpds <-
    pf$data %>%
    dplyr::select(one_of(c("Image_Metadata_BROAD_ID","Image_Metadata_SOURCE_COMPOUND_NAME"))) %>% # select ID and Compound
    mutate_if(is.character, funs(toupper)) %>% # same compound can be writen in upper and lower case and looks different# be sure it is a dataframe
    group_by(Image_Metadata_SOURCE_COMPOUND_NAME) %>%
    slice(1) %>%
    ungroup
  
  # metadata: compound, MOA, target, ID
  metadata <- 
    moa %>%
    left_join(., id.cmpds, by = "Image_Metadata_SOURCE_COMPOUND_NAME")
  
  # select only rows that have a BROAD_ID
  metadata %<>% filter(!is.na(Image_Metadata_BROAD_ID))
  
  # select MOA that are appearing more than once (meaning at least two compounds are related to it)
  n.MOA <- table(metadata$MOA_group) %>% as.data.frame() %>% filter(Freq != 1)
  metadata %<>%  filter(MOA_group %in% n.MOA$Var1)
  
  # average along replicate (keep only id and variables)
  pf.cmpds <-
    pf$data %>%
    filter(Image_Metadata_BROAD_ID %in% metadata$Image_Metadata_BROAD_ID) %>% # select ID that have a unique compound
    group_by(Image_Metadata_BROAD_ID) %>% 
    summarise_each(funs(mean(., na.rm=TRUE)), -c(Image_Metadata_SOURCE_COMPOUND_NAME, Well, Plate)) %>% # mean along replicate
    ungroup %>%
    as.data.frame()
  
  # keep track of ID of the compounds
  row.names(pf.cmpds) <- pf.cmpds$Image_Metadata_BROAD_ID
  
  # Attention: since some compounds have more than one MOAs, few rows have same compounds name.
  pf.cmpd.meta <- 
    metadata %>%
    dplyr::left_join(., pf.cmpds, by = "Image_Metadata_BROAD_ID") %>%
    as.data.frame
  
  # binary indicator matrix of ID vs MOA 
  n.moa.ID <- 
    pf.cmpd.meta %>%
    select(MOA_group, Image_Metadata_BROAD_ID) %>%
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
  N <- 1000 # number of time reproduce the same analysis
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
  enrichment.ratio <- bind_rows(enrichment.ratio, 
                                data.frame(mean = mean(random.percent), 
                                           quant = quantile(random.percent, .95), 
                                           percent = final.number, 
                                           filename = filename,
                                           method = "Pearson_random"))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken # 50 secs

save(enrichment.ratio,file="enrichment_ratio_p_random.Rda")

start.time <- Sys.time()

for(filename in data.filenames.j){
  print(filename)
  pf <- readRDS(file.path("..", "..", "input", "BBBC022_2013", "selected_single_cell_zoom", "hit_selected", "Jaccard", filename))
  profiles <- pf$data
  variables <- pf$feat_cols
  
  ## Metadata 
  # import MOAs data
  moa <- 
    read.csv("../../input/BBBC022_2013/MOAs.csv", na.strings = c("", "NA")) %>%
    mutate_if(is.factor, as.character) # has to do that to duplicate rows where multiple MOA
  
  moa %<>% 
    filter(!is.na(MOA_group)) %>% # remove rows where moa = NA
    rename(Image_Metadata_SOURCE_COMPOUND_NAME = Name) # rename compound name
  
  # remove duplicate compounds
  moa$Image_Metadata_SOURCE_COMPOUND_NAME <- toupper(moa$Image_Metadata_SOURCE_COMPOUND_NAME) # put compound names in upper character (to make it comparable)
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
  
  # mapping IDs-Compounds: give access to the compound knowing the IDs
  id.cmpds <-
    pf$data %>%
    dplyr::select(one_of(c("Image_Metadata_BROAD_ID","Image_Metadata_SOURCE_COMPOUND_NAME"))) %>% # select ID and Compound
    mutate_if(is.character, funs(toupper)) %>% # same compound can be writen in upper and lower case and looks different# be sure it is a dataframe
    group_by(Image_Metadata_SOURCE_COMPOUND_NAME) %>%
    slice(1) %>%
    ungroup
  
  # metadata: compound, MOA, target, ID
  metadata <- 
    moa %>%
    left_join(., id.cmpds, by = "Image_Metadata_SOURCE_COMPOUND_NAME")
  
  # select only rows that have a BROAD_ID
  metadata %<>% filter(!is.na(Image_Metadata_BROAD_ID))
  
  # select MOA that are appearing more than once (meaning at least two compounds are related to it)
  n.MOA <- table(metadata$MOA_group) %>% as.data.frame() %>% filter(Freq != 1)
  metadata %<>%  filter(MOA_group %in% n.MOA$Var1)
  
  # average along replicate (keep only id and variables)
  pf.cmpds <-
    pf$data %>%
    filter(Image_Metadata_BROAD_ID %in% metadata$Image_Metadata_BROAD_ID) %>% # select ID that have a unique compound
    group_by(Image_Metadata_BROAD_ID) %>% 
    summarise_each(funs(mean(., na.rm=TRUE)), -c(Image_Metadata_SOURCE_COMPOUND_NAME, Well, Plate)) %>% # mean along replicate
    ungroup %>%
    as.data.frame()
  
  # keep track of ID of the compounds
  row.names(pf.cmpds) <- pf.cmpds$Image_Metadata_BROAD_ID
  
  # Attention: since some compounds have more than one MOAs, few rows have same compounds name.
  pf.cmpd.meta <- 
    metadata %>%
    dplyr::left_join(., pf.cmpds, by = "Image_Metadata_BROAD_ID") %>%
    as.data.frame
  dim(pf.cmpd.meta)
  
  # binary indicator matrix of ID vs MOA 
  n.moa.ID <- 
    pf.cmpd.meta %>%
    select(MOA_group, Image_Metadata_BROAD_ID) %>%
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
  N <- 1000 # number of time reproduce the same analysis
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
  enrichment.ratio <- bind_rows(enrichment.ratio, 
                                data.frame(mean = mean(random.percent), 
                                           quant = quantile(random.percent, .95), 
                                           percent = final.number, 
                                           filename = filename,
                                           method = "Jaccard_random"))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken # 50 secs

save(enrichment.ratio,file="enrichment_ratio_random.Rda")



###### plot
enr.ratio <- enrichment.ratio %>% mutate(ratio = percent/mean)
enr.ratio %<>% mutate(n.feat = lapply(str_split(enr.ratio[,"filename"], "_"), function(x) x[5]) %>% unlist)
enr.ratio[,"n.feat"] %<>% as.numeric()

tmp <- 
  enr.ratio %>%
  group_by(n.feat, method) %>%
  summarise(ratio.mean = mean(ratio), 
            ratio.sd = sd(ratio)/sqrt(n()))

ggplot(tmp, aes(x=n.feat, y=ratio.mean, colour=method)) + 
  geom_errorbar(aes(ymin=ratio.mean-ratio.sd, ymax=ratio.mean+ratio.sd), width=.1) +
  geom_line() +
  geom_point() +
  labs(y = "ratio", x = "nb features") 

save(enr.ratio, file = "final_enr_ratio.Rda")
