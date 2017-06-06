#################
### CDRP Data ###
#################

library(magrittr)
library(dplyr)
library(tidyverse)
library(stringr)

################################
# Loading data + preprocessing #
################################
filename <- "../../../../input/cdrp/Pf_bio_new_all.rds"

Pf <- readRDS(filename)

variables <-
  names(Pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")

# load metadatafile for compound names
cmpd.file <- 
  read.csv("../../../../input/cdrp/cdrp.cpd.meta.csv") %>%
  filter( !grepl("BRD",CPD_NAME)) %>% # filter only compound that are bioactive (that have a compound name) (1699)
  filter(!is.na(CPD_ID)) %>% # filter NA rows (1696)
  select(one_of(c("BROAD_CPD_ID", "CPD_NAME"))) %>%
  rename(Metadata_cpd_name = CPD_NAME) # rename CPD_NAME

# join Pf and cmpd.file
Pf %<>%
  left_join(., cmpd.file, by = c("Metadata_pert_id"="BROAD_CPD_ID")) %>%
  filter(!is.na(Metadata_cpd_name)) # remove NA (12793x1645)

metadata <-
  names(Pf) %>% str_subset("^Metadata_")

######################
# Features Selection #
######################




#################
# Hit Selection #
#################

# uplodaing functions
source("hit_selection_correlation_function.R")
source("hit_selection_jaccard_function.R")

start.time <- Sys.time()

# reproductibility
set.seed(42)
N <- 20
seeds <- sample(1:10000, N, replace=F)

hit.ratio.p <- c()
i <- 0
for(s in seeds){
  print(i)
  i <- i + 1
  hit <- hit_selection_correlation(Pf, n.replicate = 8, 
                                   cor.method = "pearson", 
                                   feat.selected = F, 
                                   seed = s, 
                                   N = 5000, 
                                   dir.save = "CDRP/Profile")
  hit.ratio.p <- cbind(hit.ratio.p, hit)
}
print("mean hit ratio: ")
print(mean(hit.ratio.p))
print("standard deviation hit ratio: ")
print(sd(hit.ratio.p)/N)

hit.ratio.j <- c()
i <- 0
for(s in seeds){
  print(i)
  i <- i + 1
  hit <- hit_selection_jaccard(Pf, 
                               n.replicate = 8,
                               n.feat = 100, 
                               feat.selected = F, 
                               seed = s, 
                               N = 5000,
                               dir.save = "CDRP/Profile")
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

source("enrichment_ratio_function.R")

data.filenames.p <- list.files(path = "../../input/CDRP/Profile/hit_selected/Pearson/")
data.filenames.j <- list.files(path = "../../input/CDRP/Profile/hit_selected/Jaccard/")

# dataframe of result
enrichment.ratio <- data.frame(mean = numeric(0), quant = numeric(0), percent = numeric(0), filename = character(0), method = character(0))

for(i in 1:length(data.filenames.p)){
  seed <- str_split(data.filenames.p[i], "_") %>% unlist
  seed <- str_split(seed, ".rds") %>% unlist
  seed <- seed[4]
  
  pf.p <- readRDS(file.path("..", "..", "input", "CDRP", "Profile", "hit_selected", "Pearson", data.filenames.p[i]))
  enrichment.ratio <- bind_rows(enrichment.ratio, 
                                enrichment_ratio(pf.p, top.x = 0.02, seed = seed,
                                                 nCPU = 7, N = 1000, filename = data.filenames.p[i], method = "Pearson"))

  
  pf.j <- readRDS(file.path("..", "..", "input", "CDRP", "Profile", "hit_selected", "Jaccard", data.filenames.j[i]))  
  enrichment.ratio <- bind_rows(enrichment.ratio, 
                                enrichment_ratio(pf.j, top.x = 0.02, seed = seed, 
                                                 nCPU = 7, N = 1000, filename = data.filenames.j[i], method = "Jaccard"))
}


###### plot
enr.ratio <- enrichment.ratio %>% mutate(ratio = percent/mean)

enr.ratio %<>%
  group_by(method) %>%
  summarise(ratio.mean = mean(ratio), 
            ratio.sd = sd(ratio)/sqrt(n()))
