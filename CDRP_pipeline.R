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


######################
# Features Selection #
######################
# remove features that are too correlated (threshold of .9)
Pf %<>%
  cytominer::select(
    sample = Pf,
    variables = variables,
    operation = "correlation_threshold"
  )

#########################
# End of feat selection #
#########################
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
                                   dir.save = "CDRP/FindCorrelation")
  hit.ratio.p <- cbind(hit.ratio.p, hit)
}
print("mean hit ratio: ")
print(mean(hit.ratio.p))
print("standard deviation hit ratio: ")
print(sd(hit.ratio.p)/N)

hit.ratio.j <- c()

num.feat <- round(0.16*length(variables))

i <- 0
for(s in seeds){
  print(i)
  i <- i + 1
  hit <- hit_selection_jaccard(Pf, 
                               n.replicate = 8,
                               n.feat = num.feat, 
                               feat.selected = F, 
                               seed = s, 
                               N = 5000,
                               dir.save = "CDRP/FindCorrelation")
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

data.filenames.p <- list.files(path = "../../input/CDRP/FindCorrelation/hit_selected/Pearson/")
data.filenames.j <- list.files(path = "../../input/CDRP/FindCorrelation/hit_selected/Jaccard/")


###### TO REMOVE
#data.filenames.j <- data.filenames.j[grep(as.character(num.feat),data.filenames.j)]
######

# dataframe of result
enrichment.ratio <- data.frame(mean = numeric(0), quant = numeric(0), percent = numeric(0), filename = character(0), method = character(0))

for(i in 1:length(data.filenames.p)){
  seed <- str_split(data.filenames.p[i], "_") %>% unlist
  seed <- str_split(seed, ".rds") %>% unlist
  seed <- seed[4]
  
  pf.p <- readRDS(file.path("..", "..", "input", "CDRP", "FindCorrelation", "hit_selected", "Pearson", data.filenames.p[i]))
  enrichment.ratio <- bind_rows(enrichment.ratio, 
                                enrichment_ratio(pf.p, top.x = 0.02, seed = seed,
                                                 nCPU = 7, N = 1000, filename = data.filenames.p[i], method = "Pearson"))

  
  pf.j <- readRDS(file.path("..", "..", "input", "CDRP", "FindCorrelation", "hit_selected", "Jaccard", data.filenames.j[i]))  
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


'#
CDRP:
For 1,627 features (without any feature selection)
	method	ratio.mean	ratio.sd
1	Jaccard	2.059577		0.003770476
2	Pearson	2.313524		0.008606781

-> need independent features (here lot of them are redundant)

For findCorrelation: (687 features)
	method	ratio.mean	ratio.sd
1	Pearson	2.437757		0.008883476
2	Jaccard	2.173824		0.008427467 (5%)
3	Jaccard	2.240242		0.005313548 (10%)
4	Jaccard	2.42721	  	0.005514428 (12%)
5	Jaccard	2.462866		0.005273884 (15%)
6	Jaccard 2.21458 		0.05082299 (16%)
7	Jaccard	2.538451 	  0.01156294 (17%)
8	Jaccard	2.207021		0.00996389 (18%)
9	Jaccard	2.37801	  	0.008864711 (20%)
'