# Optimizing the Image-Based Cell Profiling Pipeline

In this [project](https://github.com/broadinstitute/imaging_metric_comparison/blob/master/Duc_Marie_thesis.pdf), we explore images of cancer cells reacting to some treatments (compounds) in a high-throughput assay.
The aim is to improve the morphological profiling pipeline steps which is converting cell images into a signature for 
each treatment.

We compare different methods to select appropriate features, to select compounds showing a strong
phenotype (names as hit selection), and to evaluate whether the profiles having the same Mechanism of Action (MOA) 
are indeed more similar.

We found out that an unsupervised feature selection based on SVD-entropy combined with a similarity metric based on 
Jaccard distance for the selection of the hits was improving the specificity of the signatures of each compound 
compared to the actual pipeline by 30\%. To evaluate this specificity, we defined a metric called the enrichment 
ratio. It uses the MOA information in order to evaluate the global accuracy of the signatures. 
We additionally defined a visualization method called ''waterfalls plot'', 
that helps to determine the optimal dose that should be used for each drug.

The feature selection based on SVD-entropy can be found in the [following paper](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/22/14/10.1093/bioinformatics/btl214/2/btl214.pdf?Expires=1503015222&Signature=ZrUCThRkrF8LTH-HaGasjkxGwTEdIEoRAYTWZkitXCDy9u9q8d8DDz8DLn95uGot6Zw0dzTGMtWz63Spx7Iw2RlRvtJDBZ6qC9duheYyOP585Bd8Eiqz0zfAwUJmi-NHuEIodDPBh4sTmrzNjzBI~Nbbt-fPTbtPqwN6GiyARBJfWJRvVHRf0kUo3R0JJiN6O8ISihw51UE6Y2C7c5Et4WzlCN5iMKLqQ9GVpm-Vf2VOF5RO5~8P4s6QLk01ApOlPHBfXlLW5DZb10-wGtb1es7FYhbevsdVENLvQXMBCSG4RicSCkvBMy3~IgeTREmvctwqBZw2ETk~uP3ZdXHB0A__&Key-Pair-Id=APKAIUCZBIA4LVPAVW3Q): 
"Novel Unsupervised Feature Filtering of Biological Data" by R. Varshavsky and al.


## Getting Started

These files give examples in order to perform the cell profiling pipeline on two different dataset (the BBBC022 dataset
and the Repurposig dataset).

The BBBC022 images can be found [here](https://data.broadinstitute.org/bbbc/BBBC022/).

The Repurposing dataset is private.

The code is each time for the two different datasets, because the MOAs annotations files is different
(separation based on `,` or `|`). Moreover the metadata are also different (the BBBC022 data is older).

Abbreviation: 
* sc = single cell
* fs = feature selection
* SVD = SVD-entropy feature selection
* fc = findCorrelation feature selection

### Prerequisites

You need to use R. 
The R libraries should be installed using: ```install.packages('name of the package')```.

The following packages are needed:
* dplyr / magrittr / ggplot2 / reshape2 / foreach / doMC / stringr / tidyverse / caret / Rcpp / RcppArmadillo

You also need to install [cytominer](https://github.com/cytomining/cytominer).

### Most important files

* **BBBC022_waterfallsPlot.Rmd**: example to create waterfalls plot on BBBC022 dataset

* **enrichment_ratio_function.R**: function to calculate the enrichment ratio on BBBC022 dataset

* **enrichment_ratio_function_repurposing.R**: function to calculate the enrichment ratio on the Repurposing dataset

* **feature_selection_SVD_entropy_cpp_old.Rmd**: example of feature selection using SVD-entropy

* **feature_selection_findCorrelation.Rmd**: example of feature selection using findCorrelation on BBBC022

* **global_pipeline.R**: example of pipeline going from hit selection (with two different metrics: 
Pearson's correlation and the Jaccard distance) to the enrichment ratio.

* **hit_selection_correlation_function.R**: hit selection function using correlation metric

* **hit_selection_jaccard_function.R**: hit selection function using the Jaccard distance metric

* **jaccard_distance_function.cpp**: c++ code for the Jaccard distance

* **random_feature_selection.R**: code to select randomly features in a dataset

* **ranking_SVD_entropy.cpp**: c++ code to do the SVD-entropy feature selection

* **repurposing_feature_selection_sc_fc.R**: example of feature selection with findCorrelation on Repurposing dataset

* **repurposing_waterfallPlot.Rmd**: example to create waterfalls plot on the Repurposing dataset

## Author

* Marie Duc

## Acknowledgments

* Dr. Mohammad Rohban
* Dr. Anne Carpenter
* All the Imaging Platform at the Broad Institute
