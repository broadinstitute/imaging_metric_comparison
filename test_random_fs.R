# select random features
source("random_feature_selection.R")

# random seed for different 
N <- 1#0
set.seed(seed = 42)
seeds <- sample(1:10000, N, replace=F)

# number of features
n.feats <- seq(from = 100, to = 1500, by = 100)

#data.filename <- "Pf_Gustafsdottir.rds"
#data.filename <- "../../input/repurposing/repurposing_normalized.rds"
#pf <- readRDS(file.path("..", "..", "input", "BBBC022_2013", "old", data.filename)) 
#Pf <- readRDS(data.filename)

# select dose at 10um
#Pf %<>% filter(Metadata_mmoles_per_liter == 10)

# there are two compounds that are removed because no MOA are associated and also DMSO
#Pf %<>% filter(!is.na(Metadata_moa))

# remove "BRD-K60230970-001-10-0@19.9991253536431" because it is a toxic compounds
#Pf %<>% filter(Metadata_broad_sample != "BRD-K60230970-001-10-0@19.9991253536431")

# remove features with NA
#Pf <- Pf[ , colSums(is.na(Pf)) == 0] # from 1784 to 1615

#metadata <-
#  colnames(Pf) %>% str_subset("^Metadata_") # loosing 1 metadata.

#variables <-
#  colnames(Pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")

pf <- readRDS(file.path("../../input/repurposing/repurposing_normalized.rds"))

# select dose at 10um
#pf %<>% filter(Metadata_mmoles_per_liter == 10)

# there are two compounds that are removed because no MOA are associated and also DMSO
pf %<>% filter(!is.na(Metadata_moa))

# remove "BRD-K60230970-001-10-0@19.9991253536431" because it is a toxic compounds
pf %<>% filter(Metadata_broad_sample != "BRD-K60230970-001-10-0")

pf %<>% 
  mutate(Metadata_broad_sample_id = Metadata_broad_sample) 
pf %<>%
  mutate(Metadata_broad_sample = paste(Metadata_broad_sample, Metadata_mmoles_per_liter, sep="@"))

# remove features with NA
pf <- pf[ , colSums(is.na(pf)) == 0] # from 1784 to 1615

# add rank to the dose, because each compounds does not have the same doses
rank <- pf[!duplicated(pf$Metadata_broad_sample),]

rank %<>%
  arrange(Metadata_mmoles_per_liter) %>%
  group_by(Metadata_broad_sample_id) %>%
  mutate(Metadata_rank=row_number()) %>%
  ungroup %>%
  select(one_of("Metadata_rank", "Metadata_broad_sample"))

pf %<>%
  left_join(., rank, by = "Metadata_broad_sample")

# select dose
dose = 6
Pf <- pf %>% filter(Metadata_rank == dose)

metadata <-
  colnames(Pf) %>% str_subset("^Metadata_") # loosing 1 metadata.

variables <-
  colnames(Pf) %>% str_subset("^Cells_|^Cytoplasm_|^Nuclei_")

for(n in n.feats){
  print("number of features")
  print(n)
  for(s in seeds){
    pf.random <- random_feature_selection(Pf, nb.feat=n, seed=s)
    
    filename.save <- 
      paste("../../input/repurposing/profile/random/feat_names/random_",
            n,
            "_dose",
            dose,
            "_seed_",
            s,
            ".rds", 
            sep = "")
    
    pf.random %>%
      saveRDS(filename.save)
  }
}
