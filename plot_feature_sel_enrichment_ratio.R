

#cell_cnt <- read_csv("workspace/scratch/code/Cell_counts.csv")

# mean across the plates
#cell_cnt %<>% select(-Image_Metadata_Plate) %>% group_by(Image_Metadata_Well) %>% summarise(Image_Count_Cells = median(Image_Count_Cells))

#split Image_Metadata_Well
#cell_cnt %<>% mutate(r = substring(Image_Metadata_Well, 1, 1)) %>% mutate(c = substring(Image_Metadata_Well, 2))

#ggplot(data = cell_cnt, aes(x = c, y = r)) +
#  geom_tile(aes(fill = Image_Count_Cells))


# for each plate


##############
library(ggplot2)
library(dplyr)
library(magrittr)

files <- list.files(path = "../../input/BBBC022_2013/Profile/enrichment_ratio/SVD_entropy/")
files.split <- lapply(strsplit(files, "_"), function(x) x[6]) %>% as.numeric() %>% unlist

final.enr.ratio <- data.frame()

for(i in 1:length(files)){
  load(file = paste("../../input/BBBC022_2013/Profile/enrichment_ratio/SVD_entropy/",
                           files[i],
                           sep="")) 
  tmp <- data.frame(n.feat = c(files.split[i], files.split[i]))
  enr.ratio %<>% bind_cols(., tmp)
  final.enr.ratio %<>% bind_rows(., enr.ratio)
}
final.enr.ratio[final.enr.ratio=="Jaccard"] <- "Jaccard_SVD_profile"
final.enr.ratio[final.enr.ratio=="Pearson"] <- "Pearson_SVD_profile"

## data random fs
file.random <- list.files(path = "../../input/BBBC022_2013/Profile/enrichment_ratio/random/")
files.split.random <- 
  lapply(strsplit(file.random, "feat.Rda"), function(x) x[1]) %>% unlist
files.split.random <- lapply(strsplit(files.split.random, "_"), function(x) x[3]) %>% as.numeric() %>% unlist

for(i in 1:length(file.random)){
  load(file = paste("../../input/BBBC022_2013/Profile/enrichment_ratio/random/",
                    file.random[i],
                    sep="")) 
  tmp <- data.frame(n.feat = c(files.split.random[i], files.split.random[i]))
  enr.ratio %<>% bind_cols(., tmp)
  final.enr.ratio %<>% bind_rows(., enr.ratio)
}
final.enr.ratio[final.enr.ratio=="Jaccard"] <- "Jaccard_random"
final.enr.ratio[final.enr.ratio=="Pearson"] <- "Pearson_random"

## data single cell fs
files <- list.files(path = "../../input/BBBC022_2013/selected_single_cell_zoom/enrichment_ratio/")
files.split <- lapply(strsplit(files, "_"), function(x) x[6]) %>% as.numeric() %>% unlist

for(i in 1:length(files)){
  load(file = paste("../../input/BBBC022_2013/selected_single_cell_zoom/enrichment_ratio/",
                    files[i],
                    sep="")) 
  tmp <- data.frame(n.feat = c(files.split[i], files.split[i]))
  enr.ratio %<>% bind_cols(., tmp)
  final.enr.ratio %<>% bind_rows(., enr.ratio)
}


load(file = "../../input/BBBC022_2013/Profile_selected/enrichment_ratio/enr_ratio_feat_sel_findCorr.Rda") 
final.enr.ratio %<>% bind_rows(., enr.ratio)

load(file ="../../input/BBBC022_2013/single_cells_findCorrelation/enrichment_ratio/enr_ratio_feat_sel_findCorr_10000.Rda")
final.enr.ratio %<>% bind_rows(., enr.ratio)

load(file ="../../input/BBBC022_2013/single_cells_findCorrelation/enrichment_ratio/enr_ratio_feat_sel_findCorr_200000.Rda")
final.enr.ratio %<>% bind_rows(., enr.ratio)

load(file = "../../input/BBBC022_2013/Profile_without_fs/enrichment_ratio/enr_ratio_no_fs.Rda")
enr.ratio[,4] <- 799
final.enr.ratio %<>% bind_rows(., enr.ratio)


ggplot(final.enr.ratio, aes(x=n.feat, y=ratio.mean, colour=method)) + 
  geom_errorbar(aes(ymin=ratio.mean-ratio.sd, ymax=ratio.mean+ratio.sd), width=.1) +
  geom_line() +
  geom_point() +
  labs(y = "enrichment ratio", x = "number of features")

ggsave("feat_sel_enrich_ratio_all.png", width = 12)





## data from 50 to 150 for random feature selection and for SVD on single cell
#load(file = "../Metric_Comparison/workspace/software/results/master/2017-06-01_20b07eb6/final_enr_ratio copy.Rda")

#enr.ratio %<>% select(one_of("method", "ratio", "n.feat"))

#tmp <- 
#  enr.ratio %>%
#  group_by(n.feat, method) %>%
#  summarise(ratio.mean = mean(ratio), 
#            ratio.sd = sd(ratio)/sqrt(n()))

#tmp <- tmp[c("method", "ratio.mean", "ratio.sd", "n.feat")]
#tmp[tmp == "Jaccard"] <- "Jaccard_SVD_single_cell"
#tmp[tmp == "Pearson"] <- "Pearson_SVD_single_cell"

#final.enr.ratio %<>% bind_rows(., tmp)

#ggplot(final.enr.ratio, aes(x=n.feat, y=ratio.mean, colour=method)) + 
#  geom_errorbar(aes(ymin=ratio.mean-ratio.sd, ymax=ratio.mean+ratio.sd), width=.1) +
#  geom_line() +
#  geom_point() +
#  labs(y = "ratio", x = "nb features")

