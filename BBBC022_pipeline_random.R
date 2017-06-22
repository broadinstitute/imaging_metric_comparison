####################
### BBBC022 Data ###
####################

library(magrittr)
library(dplyr)
library(tidyverse)
library(stringr)

filenames <- list.files(path = "../../input/BBBC022_2013/Profile/feat_selected/random/")

filenames <- c("200","250","300","350","400","500","600","700","750")

all <- 1:length(filenames)
for(n in all){
  f <- list.files(path = paste("../../input/BBBC022_2013/Profile/feat_selected/random/", filenames[n], "/", sep = ""))
  numfeat <- filenames[n]
  
  hit.ratio.p <- c()
  hit.ratio.j <- c()
  
  set.seed(42)
  N <- 10
  seeds <- sample(1:10000, N, replace=F)
  j <- 0
  for(i in f){
    j <- j + 1
    Pf <- readRDS(paste("../../input/BBBC022_2013/Profile/feat_selected/random/", 
                        filenames[n], "/",
                        i, 
                        sep = ""))
    Pf <- Pf$data %>%
      rename(Metadata_Plate = Plate,
             Metadata_Well = Well,
             Metadata_broad_sample = Image_Metadata_BROAD_ID,
             Metadata_cpd_name = Image_Metadata_SOURCE_COMPOUND_NAME)
    
    variables <-
      colnames(Pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")
    
    metadata <-
      colnames(Pf) %>% str_subset("^Metadata_")
  
    #################
    # Hit Selection #
    #################
    
    # uplodaing functions
    source("hit_selection_correlation_function.R")
    source("hit_selection_jaccard_function.R")
    
    start.time <- Sys.time()
    
    # reproductibility
    
    dir.save <- "BBBC022_2013/Profile"
    
    hit <- hit_selection_correlation(Pf, 
                                     n.replicate = 4, 
                                     filename = filenames[n],
                                     cor.method = "pearson", 
                                     feat.selected = F, 
                                     seed = seeds[j], 
                                     N = 5000, 
                                     dir.save = dir.save)
    hit.ratio.p <- cbind(hit.ratio.p, hit)
    
    num.feat <- round(0.05*length(variables))
    
    hit <- hit_selection_jaccard(Pf, 
                                 n.replicate = 4,
                                 filename = filenames[n],
                                 n.feat = num.feat, 
                                 feat.selected = F, 
                                 seed = seeds[j], 
                                 N = 5000,
                                 dir.save = dir.save)
    hit.ratio.j <- cbind(hit.ratio.j, hit)
  }
  print("mean hit ratio: ")
  print(mean(hit.ratio.p))
  print("standard deviation hit ratio: ")
  print(sd(hit.ratio.p)/n)
  
  print("mean hit ratio: ")
  print(mean(hit.ratio.j))
  print("standard deviation hit ratio: ")
  print(sd(hit.ratio.j)/N)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  
  ####################
  # Enrichment ratio #
  ####################
  
  source("enrichment_ratio_function.R")
  
  data.filenames.p <- list.files(path = "../../input/BBBC022_2013/Profile/hit_selected/random/Pearson/")
  data.filenames.j <- list.files(path = "../../input/BBBC022_2013/Profile/hit_selected/random/Jaccard/")
  
  
  ######
  data.filenames.p <- data.filenames.p[grep(as.character(numfeat),data.filenames.p)]
  data.filenames.j <- data.filenames.j[grep(as.character(numfeat),data.filenames.j)]
  
  ######
  
  # dataframe of result
  enrichment.ratio <- data.frame(mean = numeric(0), quant = numeric(0), percent = numeric(0), filename = character(0), method = character(0))
  
  seeds2 <- sample(1:10000, length(data.filenames.p), replace=F)
  
  for(i in 1:length(data.filenames.p)){
    print(i)
    seed <- seeds2[i]
    
    pf.p <- readRDS(file.path("..", "..", "input", "BBBC022_2013", "Profile", "hit_selected", "random", "Pearson", data.filenames.p[i]))
    enrichment.ratio <- bind_rows(enrichment.ratio, 
                                  enrichment_ratio(pf.p, top.x = 0.02, seed = seed,
                                                   nCPU = 7, N = 1000, filename = data.filenames.p[i], method = "Pearson"))
    
    
    pf.j <- readRDS(file.path("..", "..", "input", "BBBC022_2013", "Profile", "hit_selected", "random", "Jaccard", data.filenames.j[i]))  
    enrichment.ratio <- bind_rows(enrichment.ratio, 
                                  enrichment_ratio(pf.j, top.x = 0.02, seed = seed, 
                                                   nCPU = 7, N = 1000, filename = data.filenames.j[i], method = "Jaccard"))
  }
  
  
  ###### result
  enr.ratio <- enrichment.ratio %>% mutate(ratio = percent/mean)
  
  enr.ratio %<>%
    group_by(method) %>%
    summarise(ratio.mean = mean(ratio), 
              ratio.sd = sd(ratio)/sqrt(n()))
  enr.ratio
  
  filename.enr.ratio <- paste("../../input/BBBC022_2013/Profile/enrichment_ratio/random/enr_ratio_", 
                              filenames[n], 
                              "feat",
                              ".Rda",
                              sep = "")
  
  save(enr.ratio, file=filename.enr.ratio)
}