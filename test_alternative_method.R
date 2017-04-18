# script to plot for multiple files and multiple k the odds ratio of the Fisher's exact test


# uplodaing function
source("moa_compound_cor_function.R")

#library needed
library(magrittr)
library(dplyr)
library(ggplot2)


top.corr <- seq(from = 1, to = 15, by = 1)
files <- c("Hit_pearson_5925.rds","Hit_pearson_FS1_6369.rds","Hit_pearson_FS2_svd_6438.rds",
           "Hit_jaccard_50n_6225.rds","Hit_jaccard_30n_FS1_svd_6400.rds","Hit_jaccard_30n_FS2_svd_5950.rds") #c("Hit_pearson_5925.rds", "Hit_pearson_fs_5975.rds", "Hit_jaccard_50n_6225.rds", "Hit_jaccard_30n_fs_6206.rds")

fisher.test.res <- c()
for(filename in files){
  for(tc in top.corr){
    fisher.test.res <- cbind(fisher.test.res, alternative_method(filename, top.corr = tc))
  }  
  
}

value <- fisher.test.res %>% as.data.frame()
value <- value["estimate",] %>% unlist()

test <- data.frame(#top.corr = rep(top.corr, 2),
                   top.corr = rep(top.corr, 6),
                   val = value, 
                   method = rep(c("Pearson","Pearson with feat. sel SVD FS1", "Pearson with feat. sel SVD FS2",
                                  "Jaccard","Jaccard with feat. sel SVD FS1", "Jaccard with feat. sel SVD FS2"), each = length(top.corr)))
                   #method = rep(c("Pearson with feat. sel SVD", "Jaccard with feat. sel SVD"), each = length(top.corr)))
                   #method = rep(c("Pearson", "Pearson with feat. sel", "Jaccard", "Jaccard with feat. sel"), each = length(top.corr)))

odds_vs_top_connection_plot <- ggplot(data = test, aes(x=top.corr, y=val)) + 
  geom_line(aes(colour=method)) + 
  labs(y = "odds ratio", x = "% top connection") +
  ylim(0, 8)

ggsave(filename = 'odds_ratio_vs_top_connection_svd_250_compa.png', plot = odds_vs_top_connection_plot, width = 10, height = 6)
