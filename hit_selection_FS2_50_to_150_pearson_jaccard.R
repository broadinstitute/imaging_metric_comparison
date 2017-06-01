# hit selection function structured version

# uplodaing functions
source("hit_selection_jaccard_distance_old_function.R")
source("hit_selection_correlation_old_function.R")

start.time <- Sys.time()

data.filenames <- list.files(path = "../../input/BBBC022_2013/selected_single_cell_zoom/feat_selected/")

# reproductibility
set.seed(42)
N <- 20
seeds <- sample(1:10000, N, replace=F)

hit.ratio <- c()

for(i in 1:length(data.filenames)){
  
  pf <- readRDS(paste("../../input/BBBC022_2013/selected_single_cell_zoom/feat_selected/", data.filenames[i], sep = ""))
  
  tmp <- c()
  for(s in seeds){
    hit <- hit_selection_correlation(pf, filename = data.filenames[i], cor.method = "pearson", feat.selected = T, seed = s, N = 5000)
    hit.ratio <- cbind(hit.ratio, hit)
    tmp <- c(tmp, hit)
  }
  print("for the files: ")
  print(data.filenames[i])
  print("mean hit ratio: ")
  print(mean(tmp))
  print("standard deviation hit ratio: ")
  print(sd(tmp))
}

hit.ratio <- c()

for(i in 1:length(data.filenames)){
  
  pf <- readRDS(paste("../../input/BBBC022_2013/selected_single_cell_zoom/feat_selected/", data.filenames[i], sep = ""))
  
  tmp <- c()
  for(s in seeds){
    hit <- hit_selection_jaccard(pf, filename = data.filenames[i], n.feat = 30, feat.selected = T, seed = s, N = 5000)
    hit.ratio <- cbind(hit.ratio, hit)
    tmp <- c(tmp, hit)
  }
  print("for the files: ")
  print(data.filenames[i])
  print("mean hit ratio: ")
  print(mean(tmp))
  print("standard deviation hit ratio: ")
  print(sd(tmp))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken