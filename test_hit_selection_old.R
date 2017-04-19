# script to plot for multiple files and multiple method the hit ratio
# In order to see stabilization of the hit ratio for different seeds (random sampling)

# uplodaing functions
source("hit_selection_jaccard_distance_old_function.R")
source("hit_selection_correlation_old_function.R")

#library needed
library(magrittr)
library(dplyr)
library(ggplot2)


files <- c("Pf_Gustafsdottir.rds", "Pf_Gustafsdottir_fs_svd_FS1_250.rds", "Pf_Gustafsdottir_fs_svd_FS2_250.rds") #c("Pf_Gustafsdottir.rds", "Pf_Gustafsdottir_fs.rds")
feat.sel <- c(FALSE, TRUE, TRUE) #c(FALSE, TRUE)
n.feat <- c(50, 30, 30) #c(50, 30)
seeds <- seq(from = 1, to = 10, by = 1)

hit.ratio <- c()

for(i in 1:length(files)){
  tmp <- c()
  for(s in seeds){
    hit <- hit_selection_correlation(files[i], cor.method = "pearson", feat.selected = feat.sel[i], seed = s, N = 5000)
    hit.ratio <- cbind(hit.ratio, hit)
    tmp <- c(tmp, hit)
  }
  print("for the files: ")
  print(files[i])
  print("mean hit ratio: ")
  print(mean(tmp))
  print("standard deviation hit ratio: ")
  print(sd(tmp))
}

for(i in 1:length(files)){
  tmp <- c()
  for(s in seeds){
    hit <- hit_selection_jaccard(files[i], n.feat = n.feat[i], feat.selected = feat.sel[i], seed = s, N = 5000)
    hit.ratio <- cbind(hit.ratio, hit)
    tmp <- c(tmp, hit)
  }
  print("for the files: ")
  print(files[i])
  print("mean hit ratio: ")
  print(mean(tmp))
  print("standard deviation hit ratio: ")
  print(sd(tmp))
}

test <- data.frame(seed = rep(seeds, 6), 
                   hit.ratio = hit.ratio[1,], 
                   methods = rep(c("Pearson", "Pearson with FS1", "Pearson with FS2",
                                   "Jaccard", "Jaccard with FS1", "Jaccard with FS2"), each = length(seeds)))
#                   methods = rep(c("Pearson", "Pearson with feat. sel", "Jaccard", "Jaccard with feat. sel"), each = length(seeds)))

#test <- data.frame(seed = rep(seeds, 2), 
#                  hit.ratio = hit.ratio[1,],
#                   methods = rep(c("Pearson with feat sel SVD-entropy", "Jaccard with feat sel SVD-entropy"), each = length(seeds)))

hit_ratio_vs_seed_plot <- ggplot(data = test, aes(x=methods, y=hit.ratio)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  labs(y = "hit ratio") +
  ylim(0, 1)

ggsave(filename = 'hit_ratio_vs_seed_svd_250_compa.png', plot = hit_ratio_vs_seed_plot, width = 9, height = 6)
