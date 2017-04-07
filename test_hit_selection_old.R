# script to plot for multiple files and multiple method the hit ratio
# In order to see stabilization of the hit ratio for different seeds (random sampling)

# uplodaing functions
source("hit_selection_jaccard_distance_old_function.R")
source("hit_selection_correlation_old_function.R")

#library needed
library(magrittr)
library(dplyr)
library(ggplot2)


files <- "Pf_Gustafsdottir_fs_svd.rds" #c("Pf_Gustafsdottir.rds", "Pf_Gustafsdottir_fs.rds")
feat.sel <- TRUE #c(FALSE, TRUE)
n.feat <- 30 #c(50, 30)
seeds <- seq(from = 1, to = 100, by = 1)

hit.ratio <- c()

for(i in 1:length(files)){
  for(s in seeds){
    hit.ratio <- cbind(hit.ratio, hit_selection_correlation(files[i], cor.method = "pearson", feat.selected = feat.sel[i], seed = s, N = 5000))
  }
}

for(i in 1:length(files)){
  for(s in seeds){
    hit.ratio <- cbind(hit.ratio, hit_selection_jaccard(files[i], n.feat = n.feat[i], feat.selected = feat.sel[i], seed = s, N = 5000))
  }
}

#test <- data.frame(seed = rep(seeds, 4), 
#                   hit.ratio = hit.ratio[1,], 
#                   methods = rep(c("Pearson", "Pearson with feat. sel", "Jaccard", "Jaccard with feat. sel"), each = length(seeds)))

test <- data.frame(seed = rep(seeds, 2), 
                   hit.ratio = hit.ratio[1,],
                   methods = rep(c("Pearson with feat sel SVD-entropy", "Jaccard with feat sel SVD-entropy"), each = length(seeds)))

hit_ratio_vs_seed_plot <- ggplot(data = test, aes(x=methods, y=hit.ratio)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  labs(y = "hit ratio") +
  ylim(0, 1)

ggsave(filename = 'hit_ratio_vs_seed_svd.png', plot = hit_ratio_vs_seed_plot, width = 6, height = 6)