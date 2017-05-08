source("hit_selection_jaccard_distance_old_function.R")
source("hit_selection_correlation_old_function.R")
source("random_feature_selection.R")

start.time <- Sys.time()

pf <- readRDS(file.path("..", "..", "input", "BBBC022_2013", "old", "Pf_Gustafsdottir.rds"))
N <- 20
set.seed(42)
a <- sample(1:10000, N, replace=F)
nb.feat <- c(100, 200, 300, 400, 500, 600, 700)

seeds <- seq(from = 1, to = 5, by = 1)

mean.hit.ratio.person <- c()
sd.hit.ratio.pearson <- c()
mean.hit.ratio.jaccard <- c()
sd.hit.ratio.jaccard <- c()

# for each number of features
for(j in nb.feat){
  print("Feature selection round: ")
  print(j)
  tmp.pearson <- c()
  tmp.jaccard <- c()
  for(i in 1:N) {
    print("number of random pf: ")
    print(i)
    # select randomly nb.feat features
    pf_selected <- random_feature_selection(pf, j, seed = a[i])
    for(s in seeds){
      hit.pearson <- hit_selection_correlation(pf_selected, cor.method = "pearson", feat.selected = TRUE, seed = s, N = 5000)
      hit.jaccard <-hit_selection_jaccard(pf_selected, n.feat = 30, feat.selected = TRUE, seed = s, N = 5000)
      tmp.pearson <- c(tmp.pearson, hit.pearson)
      tmp.jaccard <- c(tmp.jaccard, hit.jaccard)
    }
  }
  
  # in case of crash save mean(tmp) and sd(tmp)
  filename.save <- paste(toString(j), 
                         "mean_and_sd_hit_ratio_pearson",
                         ".Rda", 
                         sep = "")
  mean.tmp.pearson <- mean(tmp.pearson)
  sd.tmp.pearson <- sd(tmp.pearson)
  save(mean.tmp.pearson, sd.tmp.pearson, file=filename.save)
  filename.save <- paste(toString(j), 
                         "mean_and_sd_hit_ratio_jaccard",
                         ".Rda", 
                         sep = "")
  mean.tmp.jaccard <- mean(tmp.jaccard)
  sd.tmp.jaccard <- sd(tmp.jaccard)
  save(mean.tmp.jaccard, sd.tmp.jaccard, file=filename.save)
  
  mean.hit.ratio.person <- cbind(mean.hit.ratio.person, mean(tmp.pearson))
  sd.hit.ratio.pearson <- cbind(sd.hit.ratio.pearson, sd(tmp.pearson))
  mean.hit.ratio.jaccard <- cbind(mean.hit.ratio.jaccard, mean(tmp.jaccard))
  sd.hit.ratio.jaccard <- cbind(sd.hit.ratio.jaccard, sd(tmp.jaccard))
    
  print("number of feature: ")
  print(j)
  print("mean hit ratio: ")
  print(mean(tmp.pearson))
  print(mean(tmp.jaccard))
  print("standard deviation hit ratio: ")
  print(sd(tmp.pearson))
  print(sd(tmp.jaccard))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

save(mean.hit.ratio.person, sd.hit.ratio.pearson, mean.hit.ratio.person, sd.hit.ratio.pearson, file="hit_random_fs_pearson_jaccard.Rda")
