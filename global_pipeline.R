library(magrittr)
library(dplyr)
library(tidyverse)
library(stringr)
library(parallel)

# uplodaing functions
source("hit_selection_correlation_function.R")
source("hit_selection_jaccard_function.R")
source("enrichment_ratio_function.R") ## ATTENTION!!! 2 different type
#source("enrichment_ratio_function_repurposing.R") #different MOA annotation (separation)

#dir.save <- "repurposing/single_cells/SVD"
dir.save <- "BBBC022_2013/single_cells/SR"
n.feats <- seq(from = 100, to = 700, by = 100)
# number of replicates
n.rep <- 4
nCPU <- detectCores()
#n.feats <- seq(from = 1400, to = 1500, by = 50)
#n.rep <- 5

# method that is usingSVD-entropy feature selection on single cell with SR
methodName <- "single_cells_SVD_SR"

filenames <- list.files(paste("../../input/", dir.save, "/feat_names/", sep=""))
#total.Pf <- readRDS("../../input/repurposing/profile/normalized_fc_sc.rds")
#metadata <-
#  colnames(total.Pf) %>% str_subset("^Metadata_")
# import data
#pf <- readRDS(file.path("../../input/repurposing/repurposing_normalized.rds"))
pf <- readRDS(file.path("../../input/BBBC022_2013/old/Pf_Gustafsdottir.rds"))

# difference of features name between old and new version: (this is for BBBC022)
# Hoechst -> DNA, Ph_golgi -> AGP, Syto -> RNA
pf$feat_cols <- str_replace(pf$feat_cols, "Hoechst", "DNA")
pf$feat_cols <- str_replace(pf$feat_cols, "Ph_golgi", "AGP")
pf$feat_cols <- str_replace(pf$feat_cols, "Syto", "RNA")
cl <- colnames(pf$data)
cl <- str_replace(cl, "Hoechst", "DNA")
cl <- str_replace(cl, "Ph_golgi", "AGP")
cl <- str_replace(cl, "Syto", "RNA")
colnames(pf$data) <- cl

total.Pf <- pf$data %>%
  rename(Metadata_Plate = Plate,
         Metadata_Well = Well,
         Metadata_broad_sample = Image_Metadata_BROAD_ID,
         Metadata_cpd_name = Image_Metadata_SOURCE_COMPOUND_NAME)

metadata <-
  colnames(total.Pf) %>% str_subset("^Metadata_") # loosing 1 metadata.

variables <-
  colnames(total.Pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")

for(n in n.feats){
  print(n)
  tmp <-  paste("SR_", as.character(n), "feat", sep = "")
  data.filename <- filenames[grep(tmp, filenames)]
  load(paste("../../input/", dir.save, "/feat_names/", data.filename, sep=""))
  features <- as.vector(feat.names.CEi)
  
  Pf <- total.Pf %>% select(one_of(metadata, features))
  
  variables <-
    colnames(Pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")
  
  set.seed(42)
  N <- 5
  seeds <- sample(1:10000, N, replace=F)
  
  hit.ratio.p <- c()
  i <- 0
  for(s in seeds){
    print(i)
    i <- i + 1
    hit <- hit_selection_correlation(Pf, 
                                     n.replicate = n.rep, 
                                     filename = data.filename,
                                     cor.method = "pearson", 
                                     seed = s, 
                                     nCPU = nCPU, 
                                     N = 4000, 
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
                                 n.replicate = n.rep,
                                 filename = data.filename,
                                 n.feat = num.feat,
                                 seed = s, 
                                 nCPU = dnCPU,
                                 N = 4000, 
                                 dir.save = dir.save,
                                 dir.save.plus = "/hit_selected/Jaccard/")
    hit.ratio.j <- cbind(hit.ratio.j, hit)
  }
  print("mean hit ratio: ")
  print(mean(hit.ratio.j))
  print("standard deviation hit ratio: ")
  print(sd(hit.ratio.j)/N)
  
  data.filenames.p <- list.files(path = paste("../../input/", dir.save, "/hit_selected/Pearson/", sep = ""))
  data.filenames.j <- list.files(path = paste("../../input/", dir.save, "/hit_selected/Jaccard/", sep = ""))
  
  
  ######
  tmp <- paste("SR_", as.character(n), "feat", sep = "")
  data.filenames.p <- data.filenames.p[grep(tmp,data.filenames.p)]
  data.filenames.j <- data.filenames.j[grep(tmp,data.filenames.j)]
  ######
  
  # dataframe of result
  enrichment.ratio <- data.frame(mean = numeric(0), quant = numeric(0), percent = numeric(0), filename = character(0), method = character(0))
  
  for(i in 1:length(data.filenames.p)){
    print(i)
    seed <- 42
    
    pf.p <- readRDS(paste("../../input/", dir.save, "/hit_selected/Pearson/", data.filenames.p[i], sep = ""))
    e.r.p <- enrichment_ratio(pf.p, top.x = 0.02, seed = seed,
                              nCPU = detectCores(),
                              N = 1000,
                              filename = data.filenames.p[i], method = "Pearson")
    print("enrichment ratio pearson")
    print(e.r.p$percent/e.r.p$mean)
    enrichment.ratio <- bind_rows(enrichment.ratio, e.r.p)
    
    
    pf.j <- readRDS(paste("../../input/", dir.save, "/hit_selected/Jaccard/", data.filenames.j[i], sep = ""))  
    e.r.j <- enrichment_ratio(pf.j, top.x = 0.02, seed = seed, 
                              nCPU = detectCores(),
                              N = 1000,
                              filename = data.filenames.j[i], method = "Jaccard")
    print("enrichment ratio jaccard")
    print(e.r.j$percent/e.r.j$mean)
    enrichment.ratio <- bind_rows(enrichment.ratio, e.r.j)
  }
  
  enr.ratio <- enrichment.ratio %>% mutate(ratio = percent/mean)
  
  enr.ratio %<>%
    group_by(method) %>%
    summarise(ratio.mean = mean(ratio), 
              ratio.sd = sd(ratio)/sqrt(n()))
  enr.ratio
  
  filename.enr.ratio <- paste("../../input/", dir.save, "/enrichment_ratio/", 
                              data.filenames.p[i], 
                              #".Rda",
                              sep = "")
  
  enr.ratio[enr.ratio=="Jaccard"] <- paste("Jaccard_", methodName, sep="")
  enr.ratio[enr.ratio=="Pearson"] <- paste("Pearson_", methodName, sep="")
  
  enr.ratio %<>% mutate(n.feat = n)
  
  save(enr.ratio, file=filename.enr.ratio) 
}
