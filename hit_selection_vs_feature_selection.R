# Script that calculate hit selection for different number of features selected with SVD-entropy FS2.
# and plot of the result

# uplodaing functions
source("hit_selection_jaccard_distance_old_function.R")
source("hit_selection_correlation_old_function.R")

#library needed
library(dplyr)
library(ggplot2)

filenames <- c("Pf_Gustafsdottir_sc_FS2_100.rds", "Pf_Gustafsdottir_sc_FS2_200.rds", "Pf_Gustafsdottir_sc_FS2_300.rds",
               "Pf_Gustafsdottir_sc_FS2_400.rds", "Pf_Gustafsdottir_sc_FS2_500.rds", "Pf_Gustafsdottir_sc_FS2_600.rds",
               "Pf_Gustafsdottir_sc_FS2_700.rds", "Pf_Gustafsdottir.rds")


feat.sel <- c(T, T, T, T, T, T, T, F)

n.feat <- c(30, 30, 30, 30, 30, 30, 30, 50)

seeds <- seq(from = 1, to = 10, by = 1)

hit.ratio <- c()
mean.hit.ratio <- c()
sd.hit.ratio <- c()
for(i in 1:length(filenames)){
  tmp <- c()
  for(s in seeds){
    hit <- hit_selection_correlation(filenames[i], cor.method = "pearson", feat.selected = feat.sel[i], seed = s, N = 5000, repository = "single_cells")
    hit.ratio <- cbind(hit.ratio, hit)
    tmp <- c(tmp, hit)
  }
  mean.hit.ratio <- cbind(mean.hit.ratio, mean(tmp))
  sd.hit.ratio <- cbind(sd.hit.ratio, sd(tmp))
  
  print("for the files: ")
  print(filenames[i])
  print("mean hit ratio: ")
  print(mean(tmp))
  print("standard deviation hit ratio: ")
  print(sd(tmp))
}

for(i in 1:length(filenames)){
  tmp <- c()
  for(s in seeds){
    hit <- hit_selection_jaccard(filenames[i], n.feat = n.feat[i], feat.selected = feat.sel[i], seed = s, N = 5000, repository = "single_cells")
    hit.ratio <- cbind(hit.ratio, hit)
    tmp <- c(tmp, hit)
  }
  mean.hit.ratio <- cbind(mean.hit.ratio, mean(tmp))
  sd.hit.ratio <- cbind(sd.hit.ratio, sd(tmp))
  
  print("for the files: ")
  print(filenames[i])
  print("mean hit ratio: ")
  print(mean(tmp))
  print("standard deviation hit ratio: ")
  print(sd(tmp))
}

test <- data.frame(feat.sel = rep(c(100, 200, 300, 400, 500, 600, 700, 799), 2), 
                   hit.ratio = mean.hit.ratio[1,], 
                   sd = sd.hit.ratio[1,],
                   methods = rep(c("Pearson", "Jaccard"), each = length(filenames)))

p <- ggplot(data = test, aes(x=feat.sel, y=hit.ratio, group=methods, color=methods)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=hit.ratio-sd, ymax=hit.ratio+sd), width=.2,
                position=position_dodge(0.05)) +
  labs(y = "hit ratio", x = "nb of features selected")

print(p)

ggsave(filename = 'hit_ratio_vs_feat_sel.png', plot = p, width = 9, height = 6)

save(test, file="data_feat_sel.Rda")
