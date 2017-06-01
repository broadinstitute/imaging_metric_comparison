# hit selection function structured version

# uplodaing functions
source("hit_selection_jaccard_distance_old_function.R")
source("hit_selection_correlation_old_function.R")

start.time <- Sys.time()

data.filenames <- list.files(path = "../../input/BBBC022_2013/selected_single_cell_zoom/feat_selected/")

data.filenames <- str_subset(data.filenames, "random")

hit.ratio <- c()

for(i in 1:length(data.filenames)){
  
  pf <- readRDS(paste("../../input/BBBC022_2013/selected_single_cell_zoom/feat_selected/", data.filenames[i], sep = ""))
  
  hit <- hit_selection_correlation(pf, filename = data.filenames[i], cor.method = "pearson", feat.selected = T, seed = 42, N = 5000)
  hit.ratio <- cbind(hit.ratio, hit)
  print("for the files: ")
  print(data.filenames[i])
  print("hit ratio: ")
  print(hit)
}

hit.ratio <- c()

for(i in 1:length(data.filenames)){
  
  pf <- readRDS(paste("../../input/BBBC022_2013/selected_single_cell_zoom/feat_selected/", data.filenames[i], sep = ""))
  
  hit <- hit_selection_jaccard(pf, filename = data.filenames[i], n.feat = 30, feat.selected = T, seed = 42, N = 5000)
  hit.ratio <- cbind(hit.ratio, hit)


  print("for the files: ")
  print(data.filenames[i])
  print("mean hit ratio: ")
  print(hit)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken