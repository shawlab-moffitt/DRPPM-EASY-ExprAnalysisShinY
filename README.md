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
* tidyr_1.1.3
* GSVA_1.40.1


# Required Files

## Required By User

* **Expression Matrix (.tsv/.txt):** This is required for generating the enriched signatures table in the GettingStartedScript, `GSEA_Sig_Table_Gen.R`, as well for the Shiny app. The data should be in a tab delimited format with a header where the first column consists of the gene symbols and the following columns should start with the sample name followed by the expression data. It is requested that you filter out the lowly expressed genes prior to running this analysis, as when the files are too large it can cause memory issues for some computers.
* **Meta Data (.tsv/.txt):** This is required for generating the enriched signatures table in the GettingStartedScript, `GSEA_Sig_Table_Gen.R`, as well for the Shiny app. This should be a tab delimited file with two columns, header is not necessary but you may dictate whether there is a header or not when adding your file to the `app.R` code. The first column should consist of column names and the second column should be the phenotype grouping for each sample.
* **GMT file or Gene Set Data (optional):** This is an optional input when generating the enriched signatures table in the GettingStartedScript, `GSEA_Sig_Table_Gen.R`. There is a choice for the user to user their own .gmt file which should follow the format designated by the Broad Institute as seen [here](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29). Furthermore, this .gmt file could be used as input when interacting with the Shiny app within the GSEA tab. Additionally, there is an option for the user to input a gene set file which would be a tab delimited, two column file. The first column being the gene set name repeating while the second column has a single gene symbol per row matched with the gene set name that it is within. 
  *  If the user does not opt to include their own .gmt file or gene set data, the user may designate whether the samples are from a human or a mouse model and an MSigDB gmt file will be loaded consisting of all the MSigDB collections. Please note, this is a larger data set and  running the GSEA function can take a few minutes. Lastly, there is a commented out segment where the user can derive their own MSigDB gene set from the `msigdbr()` function if they may request a species other than the human or mouse model. To retrieve a list of available species users can run the `msigdbr_species()` function.

## Required and Provided

* **MSigDB Files:** These files are provided in the /data/ folder which begin with "msigdb_" and are formated for the GSEA analysis within the GettingStartedScript, `GSEA_Sig_Table_Gen.R` and the Shiny app, as well as are used for the user interface tables to select gene sets within the GSEA tab. The msigdb_gsNsym.tsv files are zipped due to size constraints, when unzipped they are around 140mb-160mb each. The MSigDB gene sets are available for *Homo sapiens* and *Mus musculus* which can either both be downloaded or only the species you need. You choose if the samples are human or not when setting up the scripts.
* **Tab2 Gene Set GMT file:** This is a temporary file which is a subset of the MSigDB gene set collection. In the future it will be replaced with relevant gene sets that users may use when interacting with their data. Please keep in mind that this is from the Homo Sapien collection and will not generate results if you are observing mouse model data. The current gene set in this place is also included in the main MSigDB file, so idealy this tab can be ignored until it contains a more unique gene set.

# Getting Started Script

One of the required file inputs for the Shiny App is the GSEA enriched signatures table for your samples. The production of this table can take several minutes depending on the number of gene sets that are being ranked, but to save time we wrote a [script](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/GettingStartedScript/GSEA_Sig_Table_Gen.R) that requires just a few file inputs and does the rest for you! We use all of the collections from the Molecular Signatures Database (MSigDB) for the initial ranking but the getting started script and the Shiny App both allow for user input of gene set/gmt files.

## User Input for Getting Started Script

More details on these file inputs in the [Required Files](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/README.md#required-files) section.

* `expr_file` - Expression matrix file
* `meta_file` - Meta File
  * Please note `TRUE` or `FALSE` in `header` if there is a header in the meta file
* `OutPath` - Text input of user defined outfile path - This is where your enriched signature table(s) will be written to. The specific file naming is taken care of in the script.
* Gene Set File
  * `MSigDB_file` - MSigDB Gene Set
    * If you would like to use an [MSigDB gene set](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/tree/main/data) you must download and unzip the appropriate msigdb_gsNsym_\*.tsv file for either *Homo sapiens* (HS) or *Mus musculus* (MM).
    * There is code included but commented out in regards to the msigdbr package which allows you to retreive gene sets based on certain species as well as from specific collections.
    * Please note `TRUE` or `FALSE` if the sample data is human or not, so the correct gene set is located and used.
  * User Provided Gene Set
    * `GeneSet_file.u.gmt` - user provided .gmt file
    * `GeneSet_file.u.gs` - user provided geneset file
      * Please note `TRUE` or `FALSE` in `header.gs` if there is a header in the gene set file.

## Generating Enriched Signature Table(s)

Once these files are input in the script it can be run in its entirety. The script consists of reading the files and correcting for any formatting issues, generating groups based off of the meta file, and running through the signal-to-noise calculation. Based off the groups found in the meta file, an enriched signatures table will be made for each combination of the groupings and writen to a file containing the group names being compared.


# Prepping the R Shiny App

Similar to generating the Enriched Signatures Table, the Shiny App requires just a few user inputs and it can be up and running.


