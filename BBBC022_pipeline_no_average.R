####################
### BBBC022 Data ###
####################

library(magrittr)
library(dplyr)
library(tidyverse)
library(stringr)

################################
# Loading data + preprocessing #
################################
filename <- "../../../../input/BBBC022/Pf_Gustafsdottir.rds"

Pf <- readRDS(filename)
Pf <- Pf$data %>% 
  dplyr::rename(Metadata_Plate = Plate,
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
set.seed(42)
N <- 10
seeds <- sample(1:10000, N, replace=F)
dir.save <- "BBBC022_2013/Profile_without_fs"

hit.ratio.p <- c()
i <- 0
for(s in seeds){
  print(i)
  i <- i + 1
  hit <- hit_selection_correlation(Pf, 
                                   n.replicate = 4, 
                                   #filename = filenames[n],
                                   cor.method = "pearson", 
                                   feat.selected = F, 
                                   seed = s, 
                                   N = 5000, 
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
                               n.replicate = 4,
                               #filename = filenames[n],
                               n.feat = num.feat, 
                               feat.selected = T, 
                               seed = s, 
                               N = 5000,
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

source("enrichment_ratio_function_no_average.R")

data.filenames.p <- list.files(path = "../../input/BBBC022_2013/Profile_without_fs/hit_selected/Pearson/")
data.filenames.j <- list.files(path = "../../input/BBBC022_2013/Profile_without_fs/hit_selected/Jaccard/")

# dataframe of result
enrichment.ratio <- data.frame(mean = numeric(0), quant = numeric(0), percent = numeric(0), filename = character(0), method = character(0))

for(i in 1:length(data.filenames.p)){
  print(i)
  seed <- str_split(data.filenames.p[i], "_") %>% unlist
  seed <- str_split(seed, ".rds") %>% unlist
  seed <- seed[length(seed)-1]
  
  pf.p <- readRDS(file.path("..", "..", "input", "BBBC022_2013", "Profile_without_fs", "hit_selected", "Pearson", data.filenames.p[i]))
  enrichment.ratio <- bind_rows(enrichment.ratio, 
                                enrichment_ratio(pf.p, top.x = 0.02, seed = seed,
                                                 nCPU = 7, N = 1000, filename = data.filenames.p[i], method = "Pearson"))
  
  
  pf.j <- readRDS(file.path("..", "..", "input", "BBBC022_2013", "Profile_without_fs", "hit_selected", "Jaccard", data.filenames.j[i]))  
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

filename.enr.ratio <- paste("../../input/BBBC022_2013/Profile_without_fs/enrichment_ratio/enr_ratio_no_fs_median_", 
                            #data.filenames.p[i], 
                            ".Rda",
                            sep = "")

enr.ratio[enr.ratio=="Jaccard"] <- "Jaccard_profile_no_average_median"
enr.ratio[enr.ratio=="Pearson"] <- "Pearson_profile_no_average_median"

enr.ratio %<>% mutate(n.feat = dim(Pf)[2])

save(enr.ratio, file=filename.enr.ratio)


#Whole dataset:

#  method    ratio.mean  ratio.sd
#1 Jaccard   1.824503    0.015004642
#2 Pearson   1.874721    0.004085244
