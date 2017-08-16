####################
### BBBC022 Data ###
####################

library(magrittr)
library(dplyr)
library(tidyverse)
library(stringr)
library(parallel)

# uplodaing functions
source("hit_selection_correlation_function.R")
source("hit_selection_jaccard_function.R")
source("enrichment_ratio_function.R")

################################
# Loading data + preprocessing #
################################

nCPU <- detectCores() -1

# number of replicates
nrep <- 4

# number of samples 
nsample <- 5000

dir.save <- "BBBC022_2013/selected_single_cell_zoom"
inputDir <- paste("../../input/", dir.save, "/", sep="")
filenames <- list.files(path = paste(inputDir, "feat_selected/", sep = ""))

all <- 1:length(filenames)
for(n in all){
  numfeat <- str_split(filenames[n], ".rds") %>% unlist
  numfeat <- str_split(numfeat, "_") %>% unlist
  numfeat <- numfeat[length(numfeat)-1]
  Pf <- readRDS(paste(inputDir, 
                      "feat_selected/",
                      filenames[n], 
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
  
  start.time <- Sys.time()
  
  # reproductibility
  set.seed(42)
  N <- 10
  seeds <- sample(1:10000, N, replace=F)
  
  hit.ratio.p <- c()
  i <- 0
  for(s in seeds){
    print(i)
    i <- i + 1
    hit <- hit_selection_correlation(Pf, 
                                     n.replicate = nrep, 
                                     filename = filenames[n],
                                     cor.method = "pearson", 
                                     feat.selected = F, 
                                     seed = s, 
                                     N = nsample, 
                                     nCPU = nCPU,
                                     dir.save = dir.save,
                                     dir.save.plus = "/hit_selected/Pearson/")
    hit.ratio.p <- cbind(hit.ratio.p, hit)
  }
  print("mean hit ratio: ")
  print(mean(hit.ratio.p))
  print("standard deviation hit ratio: ")
  print(sd(hit.ratio.p)/N)
  
  hit.ratio.j <- c()
  
  num.feat <- round(0.05*length(variables))
  
  i <- 0
  for(s in seeds){
    print(i)
    i <- i + 1
    hit <- hit_selection_jaccard(Pf, 
                                 n.replicate = nrep,
                                 filename = filenames[n],
                                 n.feat = num.feat, 
                                 feat.selected = F, 
                                 seed = s, 
                                 N = nsample,
                                 nCPU = nCPU,
                                 dir.save = dir.save,
                                 dir.save.plus = "/hit_selected/Jaccard/")
    hit.ratio.j <- cbind(hit.ratio.j, hit)
  }
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
  
  data.filenames.p <- list.files(path = paste(inputDir, "hit_selected/Pearson/", sep=""))
  data.filenames.j <- list.files(path = paste(inputDir, "hit_selected/Jaccard/", sep=""))
  
  
  ######
  #tmp <- paste("BBC022_fs_SVD_", as.character(numfeat), sep = "")
  #data.filenames.p <- data.filenames.p[grep(tmp,data.filenames.p)]
  #data.filenames.j <- data.filenames.j[grep(tmp,data.filenames.j)]
  ######
  
  # dataframe of result
  enrichment.ratio <- data.frame(mean = numeric(0), quant = numeric(0), percent = numeric(0), filename = character(0), method = character(0))
  
  for(i in 1:length(data.filenames.p)){
    print(i)
    seed <- str_split(data.filenames.p[i], "_") %>% unlist
    seed <- str_split(seed, ".rds") %>% unlist
    seed <- seed[length(seed)-1]
    
    pf.p <- readRDS(paste(inputDir, "hit_selected/Pearson/", data.filenames.p[i], sep=""))
    enrichment.ratio <- bind_rows(enrichment.ratio, 
                                  enrichment_ratio(pf.p, top.x = 0.02, seed = seed,
                                                   nCPU = nCPU, N = 1000, filename = data.filenames.p[i], method = "Pearson_SVD_single_cell"))
    
    
    pf.j <- readRDS(paste(inputDir, "hit_selected/Jaccard/", data.filenames.j[i]))
    enrichment.ratio <- bind_rows(enrichment.ratio, 
                                  enrichment_ratio(pf.j, top.x = 0.02, seed = seed, 
                                                   nCPU = nCPU, N = 1000, filename = data.filenames.j[i], method = "Jaccard_SVD_single_cell"))
  }
  
  
  ###### result
  enr.ratio <- enrichment.ratio %>% mutate(ratio = percent/mean)
  
  enr.ratio %<>%
    group_by(method) %>%
    summarise(ratio.mean = mean(ratio), 
              ratio.sd = sd(ratio)/sqrt(n()))
  enr.ratio
  
  filename.enr.ratio <- paste(inputDir,
                              "enrichment_ratio/enr_ratio_", 
                              data.filenames.p[i], 
                              ".Rda",
                              sep = "")
  
  save(enr.ratio, file=filename.enr.ratio)
}