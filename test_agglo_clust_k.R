
# uplodaing function
source("agglomerative_clustering_old_function.R")

#library needed
library(magrittr)
library(dplyr)
library(ggplot2)


k <- seq(from = 10, to = 180, by = 1)
files <- c("Hit_pearson_5925.rds", "Hit_pearson_fs_5975.rds", "Hit_jaccard_50n_6225.rds", "Hit_jaccard_30n_fs_6206.rds")

fisher.test.res <- c()
for(filename in files){
  for(n.clust in k){
    fisher.test.res <- cbind(fisher.test.res, agglomarative_clustering(filename, n.clust))
  }  
  
}

value <- fisher.test.res %>% as.data.frame()
value <- value["estimate",] %>% unlist()

test <- data.frame(k = rep(k, 4), 
                   val = value, 
                   method = rep(c("Pearson", "Pearson with feat. sel", "Jaccard", "Jaccard with feat. sel"), each = 5))

ggplot(data = test, aes(x=k, y=val)) + 
  geom_line(aes(colour=method)) + 
  labs(y = "odds ratio")
