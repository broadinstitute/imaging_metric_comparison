# script to plot for multiple files and multiple k the odds ratio of the Fisher's exact test


# uplodaing function
source("agglomerative_clustering_nclust_function.R")

#library needed
library(magrittr)
library(dplyr)
library(ggplot2)


k <- seq(from = 10, to = 160, by = 1)
files <- c("Hit_pearson_5925.rds","Hit_pearson_FS1_6369.rds","Hit_pearson_FS2_svd_6438.rds",
           "Hit_jaccard_50n_6225.rds","Hit_jaccard_30n_FS1_svd_6400.rds","Hit_jaccard_30n_FS2_svd_5950.rds") #c("Hit_pearson_5925.rds", "Hit_pearson_fs_5975.rds", "Hit_jaccard_50n_6225.rds", "Hit_jaccard_30n_fs_6206.rds")

fisher.test.res <- c()
for(filename in files){
  for(n.clust in k){
    fisher.test.res <- cbind(fisher.test.res, agglomarative_clustering(filename, n.clust))
  }  
  
}

value <- fisher.test.res %>% as.data.frame()
value <- value["estimate",] %>% unlist()

test <- data.frame(#k = rep(k, 4), 
                   k = rep(k, 6), 
                   val = value, 
                   method = rep(c("Pearson", "Pearson with FS1", "Pearson with FS2", 
                                  "Jaccard","Jaccard with FS1", "Jaccard with FS2"), each = length(k)))
#                   method = rep(c("Pearson", "Pearson with feat. sel", "Jaccard", "Jaccard with feat. sel"), each = length(k)))
#test <- data.frame(k = rep(k, 2), 
#                   val = value,
#                   method = rep(c("Pearson with feat. sel SVD", "Jaccard with feat. sel SVD"), each = length(k)))

odds_vs_k_plot <- ggplot(data = test, aes(x=k, y=val)) + 
  geom_line(aes(colour=method)) + 
  labs(y = "odds ratio") +
  ylim(0, 5)

ggsave(filename = 'odds_ratio_vs_k_svd_250_compa.png', plot = odds_vs_k_plot, width = 10, height = 6)
