library(ggplot2)
library(dplyr)
library(magrittr)

final.enr.ratio <- data.frame()

## data single cell fs
files <- list.files(path = "../../input/BBBC022_2013/selected_single_cell_zoom/10000/enrichment_ratio/")
files.split <- lapply(strsplit(files, "_"), function(x) x[6]) %>% as.numeric() %>% unlist

for(i in 1:length(files)){
  load(file = paste("../../input/BBBC022_2013/selected_single_cell_zoom/10000/enrichment_ratio/",
                    files[i],
                    sep="")) 
  tmp <- data.frame(n.feat = c(files.split[i], files.split[i]))
  enr.ratio %<>% bind_cols(., tmp)
  final.enr.ratio %<>% bind_rows(., enr.ratio)
}

load(file ="../../input/BBBC022_2013/single_cells_findCorrelation/enrichment_ratio/enr_ratio_feat_sel_findCorr_10000.Rda")
final.enr.ratio %<>% bind_rows(., enr.ratio)

final.enr.ratio$method <- paste(final.enr.ratio$method, "_10000", sep="")

files <- list.files(path = "../../input/BBBC022_2013/selected_single_cell_zoom/200000/enrichment_ratio/")
files.split <- strsplit(files, ".Rda") %>% unlist
files.split <- lapply(strsplit(files.split, "_"), function(x) x[4]) %>% as.numeric() %>% unlist

for(i in 1:length(files)){
  load(file = paste("../../input/BBBC022_2013/selected_single_cell_zoom/200000/enrichment_ratio/",
                    files[i],
                    sep=""))
  enr.ratio[,4] <- sapply(enr.ratio[,4], as.numeric)
  final.enr.ratio %<>% bind_rows(., enr.ratio)
}

load(file ="../../input/BBBC022_2013/single_cells_findCorrelation/enrichment_ratio/enr_ratio_feat_sel_findCorr_200000.Rda")
final.enr.ratio %<>% bind_rows(., enr.ratio)

#load(file = "../../input/BBBC022_2013/Profile_without_fs/enrichment_ratio/enr_ratio_no_fs.Rda")
#enr.ratio[,4] <- 799
#final.enr.ratio %<>% bind_rows(., enr.ratio)


ggplot(final.enr.ratio, aes(x=n.feat, y=ratio.mean, colour=method)) + 
  geom_errorbar(aes(ymin=ratio.mean-ratio.sd, ymax=ratio.mean+ratio.sd), width=.1) +
  geom_line() +
  geom_point() +
  labs(y = "enrichment ratio", x = "number of features") + 
  xlim(100, 800) +
  ylim(1, 3)

ggsave("10k_vs_200kcells.png", width = 12)
