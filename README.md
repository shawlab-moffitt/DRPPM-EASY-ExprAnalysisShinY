# RShinyAnalysisGenerator

# Introduction

# Requirements

* `R` - https://cran.r-project.org/src/base/R-4/
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

# Installation

1. 

# R Dependencies

* shiny_1.6.0
* shinythemes_1.2.0
* shinyjqui_0.4.0
* pheatmap_1.0.12
* RColorBrewer_1.1-2
* clusterProfiler_4.0.5
* dplyr_1.0.7
* DT_0.18
* enrichplot_1.12.2
* shinycssloaders_1.0.0
* ggplot2_3.3.5
* ggpubr_0.4.0
* msigdbr_7.4.1
* reshape2_1.4.4
* tibble_3.1.3
* plotly_4.9.4.1
* readr_2.0.1
* limma_3.48.3
* enrichR_3.0
* ggrepel_0.9.1


# Required Files

## Required By User

* **Expression Matrix:** Data in a tab delimited format with a header where the first column consists of the gene symbols and the following columns should start with the sample name followed by the expression data. It is requested that you filter out the lowly expressed genes prior to running this analysis, as when the files are too large it can cause memory issues for some computers.
* **Meta Data:** This should be a tab delimited file with two columns, header is not necessary but you may dictate whether there is a header or not when adding your file to the `app.R` code. The first column should consist of column names and the second column should be the phenotype grouping for each sample.
* **GMT file or Gene Set Data (optional):** There is a choice for the user to user their own .gmt file which should follow the format designated by the Broad Institute as seen [here](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29). Additionally, there is an option for the user to input a gene set file which would be a tab delimited, two column file. The first column being the gene set name repeating while the second column has a single gene symbol per row matched with the gene set name that it is within. If the user does not opt to include their own .gmt file or gene set data, the user may designate whether the samples are from a human or a mouse model and an MSigDB gmt file will be loaded consisting of all the MSigDB collections. Please note, this is a larger data set and  running the GSEA function can take a few minutes. Lastly, there is a commented out segment where the user can derive their own MSigDB gene set from the `msigdbr()` function if they may request a species other than the human or mouse model. To retrieve a list of available species users can run the `msigdbr_species()` function.
