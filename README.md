# DRPPM Expression Analysis Shiny App

# Introduction

With the rapid generation of large datasets as a result of advancing next-generation sequecing (NGS) technologies, the need for quick and reproducible data analysis tools has become paramount. The ability to analyze both genomic and proteomic data sets can allow for the detection of compelling genes and features that may be of interest for further analysis. The General Expression Analysis App hosted in R Shiny does not require an extensive computation background to run these analyses. With user input of expression and meta data, along with gene set we have sourced and provided, the user may produce useful anaylsis and visualizations on a web-based interface within minutes. This method of reproducible data analysis generates a number of visualization to view Gene Set Enrichment (GSEA) and RNAseq analysis, as well as the ability to download various tables and gene sets that are produced by the App for further use.

# Installation

* Download ZIP file from https://github.com/shawlab-moffitt/GeneralExpressionAnalysisApp
* Unzip and load into directory as a project in R Studio
* Open the ‘App.R’ script and write in user input files and options as directed at the top of the script
  * ‘App.R’ script begins with example files loaded in from the ExampleData folder
* Press ‘Run App’ button in R Studio to run in application or browser window and enjoy!

# Requirements

* `R` - https://cran.r-project.org/src/base/R-4/
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

# R Dependencies

|  |  |  |  |  |
| --- | --- | --- | --- | --- |
| shiny_1.6.0 | shinythemes_1.2.0 | shinyjqui_0.4.0 | shinycssloaders_1.0.0 | tools_4.1.0 |
| dplyr_1.0.7 | tidyr_1.1.3 | readr_2.0.1 | tibble_3.1.3 | DT_0.18 |
| ggplot2_3.3.5 | plotly_4.9.4.1 | enrichplot_1.12.2 | pheatmap_1.0.12 | ggrepel_0.9.1 |
| enrichR_3.0 | limma_3.48.3 | clusterProfiler_4.0.5 | limma_3.48.3 | GSVA_1.40.1 |
| BiocManager_1.30.16 | reshape2_1.4.4 | ggpubr_0.4.0 |  |  |

# Required Files

## Required By User

* **Expression Matrix (.tsv/.txt):** This is required for generating the enriched signatures table in the GettingStartedScript, `GSEA_Sig_Table_Gen.R`, as well for the Shiny app. The data should be in a tab delimited format with a header where the first column consists of the gene symbols and the following columns should start with the sample name followed by the expression data. It is requested that you filter out the lowly expressed genes prior to running this analysis, as when the files are too large it can cause memory issues for some computers.
* **Meta Data (.tsv/.txt):** This is required for generating the enriched signatures table in the GettingStartedScript, `GSEA_Sig_Table_Gen.R`, as well for the Shiny app. This should be a tab delimited file with two columns, header is not necessary but you may dictate whether there is a header or not when adding your file to the `app.R` code. The first column should consist of column names and the second column should be the phenotype grouping for each sample.
* **GMT file or Gene Set Data (optional):** This is an optional input when generating the enriched signatures table in the GettingStartedScript, `GSEA_Sig_Table_Gen.R`. There is a choice for the user to user their own .gmt file which should follow the format designated by the Broad Institute as seen [here](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29). Furthermore, this .gmt file could be used as input when interacting with the Shiny app within the GSEA tab. Additionally, there is an option for the user to input a gene set file which would be a tab delimited, two column file. The first column being the gene set name repeating while the second column has a single gene symbol per row matched with the gene set name that it is within. 
  *  If the user does not opt to include their own .gmt file or gene set data, the user may designate whether the samples are from a human or a mouse model and an MSigDB gmt file will be loaded consisting of all the MSigDB collections. Please note, this is a larger data set and  running the GSEA function can take a few minutes. Lastly, there is a commented out segment where the user can derive their own MSigDB gene set from the `msigdbr()` function if they may request a species other than the human or mouse model. To retrieve a list of available species users can run the `msigdbr_species()` function.

## Required and Provided

* **MSigDB Files:** These files are provided in the /GeneSets/ folder which begin with "msigdb_" and are formated for the GSEA analysis within the GettingStartedScript, `GSEA_Sig_Table_Gen.R` and the Shiny app, as well as are used for the user interface tables to select gene sets within the GSEA tab. The msigdb_gsNsym.tsv files are zipped due to size constraints, when unzipped they are around 140mb-160mb each. The MSigDB gene sets are available for *Homo sapiens* and *Mus musculus* which can either both be downloaded or only the species you need. You choose if the samples are human or not when setting up the scripts.
* **Tab2 Gene Set GMT file:** This is a temporary file which is a subset of the MSigDB gene set collection. In the future it will be replaced with relevant gene sets that users may use when interacting with their data. Please keep in mind that this is from the Homo Sapien collection and will not generate results if you are observing mouse model data. The current gene set in this place is also included in the main MSigDB file, so idealy this tab can be ignored until it contains a more unique gene set.

# Getting Started Scripts

One of the required file inputs for the Shiny App is the GSEA enriched signatures table for your samples. The production of this table can take several minutes depending on the number of gene sets that are being ranked, but to save time we wrote a [script](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/GettingStartedScripts/GSEA_Sig_Table_Gen.R) that requires just a few file inputs and does the rest for you! We use all of the collections from the Molecular Signatures Database (MSigDB) for the initial ranking but the getting started script and the Shiny App both allow for user input of gene set/gmt files.

### User Input for GSEA_Sig_Table_Gen.R Script

More details on these file inputs in the [Required Files](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/README.md#required-files) section.

* `expr_file` - Expression matrix file
* `meta_file` - Meta File
  * Please note `TRUE` or `FALSE` in `header` if there is a header in the meta file
* `OutPath` - Text input of user defined outfile path - This is where your enriched signature table(s) will be written to. The specific file naming is taken care of in the script.
* Gene Set File
  * `MSigDB_file` - MSigDB Gene Set
    * If you would like to use an [MSigDB gene set](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/tree/main/GeneSets) you must download and unzip the appropriate msigdb_gsNsym_\*.tsv file for either *Homo sapiens* (HS) or *Mus musculus* (MM).
    * There is code included but commented out in regards to the msigdbr package which allows you to retreive gene sets based on certain species as well as from specific collections.
    * Please note `TRUE` or `FALSE` if the sample data is human or not, so the correct gene set is located and used.
  * User Provided Gene Set
    * `GeneSet_file.u.gmt` - user provided .gmt file
    * `GeneSet_file.u.gs` - user provided geneset file
      * Please note `TRUE` or `FALSE` in `header.gs` if there is a header in the gene set file.

### Generating Enriched Signature Table(s)

Once these files are input in the script it can be run in its entirety. The script consists of reading the files and correcting for any formatting issues, generating groups based off of the meta file, and running through the signal-to-noise calculation. Based off the groups found in the meta file, an enriched signatures table will be made for each combination of the groupings and writen to a file containing the group names being compared.

### Generating a Gene Set RData List for ssGSEA

If you choose to use your own Gene Sets either in .gmt or tab delimited format as described above, in order to perform ssGSEA analysis the Gene Set must be converted into an RData list object when loaded into the app. This list can take several minutes to generate depending on the size of the gene set file, so there is a separate script to perform this with [GeneSetRDataListGen.R](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/GettingStartedScripts/GeneSetRDataListGen.R). The only user input is the Gene Set file path and name, whether or not it has a header, and the desired outfile path and name. Once they are input, the code can be run as a whole and it will produce and save an RData list which can be input to the R Shiny app. These RData lists have already been generated for provided Gene Sets and can be found in the [data folder](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/tree/main/GeneSets) of this repsitory.

# Prepping the R Shiny App

Similar to generating the Enriched Signatures Table, the Shiny App requires just a few user inputs and it can be up and running. Once the [app](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/App/app.R) script is loaded all that is needed is to update the user input and then it can be run!

### User Input for R Shiny App

Below is the begining chunk of code for the Shiny App and where the user will designate which files to use for the analysis. Here you designate the name for the app, your expression and meta data, a single or list of enriched signature tables to display, as well as your gene set seletion for the second tab of the GSEA section. Along with ansering a could `TRUE or FALSE` statements in this chunk, it should be the only secion that requires editing by the user, though to note, there is a section below in the script that calls on the MSigDB gene sets based on the 'Human' `TRUE or FALSE` statement. These files so not have paths so the app should either be ran in the same working directory as the data input files or you may add the paths to those file names.

```{r}

####----User Data Input----####

#Input desired project name for webpage - will be followed by 'RNAseq Analysis'
ProjectName <- "USP7 Human Demo"

##--User Input File Names--##

#expression data
expr_file <- "htseq_gene_level_fpkm_T_geneName_max_1cutoff_v2.txt"

#meta data
meta_file <- "USP7_meta.tsv"
#Is there a header?
header <- TRUE

#Enriched Signatures data table
ES_tables <- c("USP7_Enrich_Sig.tsv")

#If human: set TRUE
#If mouse: set FALSE
human <- TRUE

##--User Gene Set Input--##

#write in the name of your gene set list for shiny UI
userGSlist_name <- 'CellMarker Gene Sets'

#path to your gene set file .gmt or .txt/.tsv
userGS_file <- 'CellMarker_gsNsym_HS.tsv'
#Does gene set file have header?
header.gs <- TRUE

#path to your R data list object for ssGSEA
userRData_file <- 'CellMarker_GS_HS.RData'

```

Once the file names and choices are correctly written you should be able to hit the "Run App" button in your R Studio desktop application load up the web page if you decide to save the code to a server.

# Shiny App Features

## Gene Set Enrichment Analysis

In app, the user may select two comparison groups within the side panel which are based off of the meta file that was provided, as well a gene set by selecting from one of the gene set tables in the side panel. If the user uploads their own gene set data a table will appear based on the file input. 

**Please note:** Depending on the designated P-value selected as well as the genes in the expression matrix, there will be some gene sets that may show errors. This is due to the genes within the gene set not being significantly enriched or available to produce any data or visualizations.

![alt text](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/App_Demo_Pictures/GSEA_Shiny_1.png?raw=true)

1. Select comparison groups base on types of samples from the meta data file
2. Select stat comparison method for ssGSEA boxlots
    * With 2 sample types it allows for Wilcos and t-test, with 3 or more types it allows for Kruskal-Wallis and ANOVA
3. You may use your own gene set either in .gmt, .tsv. or .txt format as described in the [Required Files](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/README.md#required-files) section
4. You may search for a collection or gene set of interest
5. The Normalized Enrichment Score (NES) and P-value appear at the top of the enrichment plot, as well as a message stating whether the selected gene set is up or down regulated in a sample type
6. Below the plot will show a ranked list of leading edge genes based on a signal-2-noise calculation, which may be downloaded for the users desired use

![alt text](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/App_Demo_Pictures/GSEA_Shiny_2.png?raw=true)

1. Based on the list of enriched signatures tables loaded into the app, their basenames will show up as a list that you may select from to view.
    * Due to this, it is helpful to be descriptive in naming the files. The table generator script does take this into account when nameing the output files

![alt text](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/App_Demo_Pictures/GSEA_Shiny_3.png?raw=true)

1. This tab allows you to upload your own gene set file and when uploaded it will display the list of gene sets within the file below. These can be selected to be viewed with other visualizations as well, such as the enrichment or volcano plot
2. This screenshot is showcasing the enriched signature table generation in app if you choose to perform it there rather than prior to loading the app
    * This is useful if you have a large number of sample types or if you have a specific comparison you are interested in with the gene set of your choice uploaded
    * Please keep in mind that depending on the size of the gene set data or number of samples this may take several minutes

![alt text](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/App_Demo_Pictures/GSEA_Shiny_8.png?raw=true)

1. When viewing the ssGSEA box plot for a specific gene set, the ssGSEA enrichment score table may be downloaded
    * This can be used further in the ssGSEA analysis app that is being developed

## RNAseq Analysis

![alt text](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/App_Demo_Pictures/GSEA_Shiny_4.png?raw=true)

These will be the RNAseq analysis parameters shown on the main side bar of this tab.

1. Here you may select comparison sample types for the volcano and MA plots, as well as the pathway analysis and DEG table
2. You may select the number of top genes to see displayed on the main heatmap
3. Here you have the choice of variance measure between MAD, CV, and VAR
4. The samples are able to to be clustered with multiple methods in a drop down box
5. The Most Variable Genes will appear in the side bar below the sample cluster results
    * This table may be downloaded as a tab delimited rank file or in .gmt format to be used as a gene set.
6. Below the Most Variable Gene table are parameters to download a .gmt file based on the the Differentially Expressed Genes table

![alt text](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/App_Demo_Pictures/GSEA_Shiny_5.png?raw=true)

1. Users may select specific genes from a list generated from the expression matrix.
    * The starting values selected are genes involved in cytokines.
    * The starting values may be deleted with backspace and new gene options will appear in a dropdown box when the box is selected or a gene of choice is typed in the box
2. The user may also paste a list of genes in a text input box and those genes will be visualized in the heatmap
    * Please note, the genes selected in the first and second box of this side panel are stacked and will show the genes selected in both boxes
3. The heatmap will begin with all samples displayed but the user may delete ones that they do not wish to see in the Sample Selection box

![alt text](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/App_Demo_Pictures/GSEA_Shiny_6.png?raw=true)

1. To adjust the logFC the user must enter the absolute value of the desired fold change and it will adjust the vertical dashed line on the volcano plot and the horizontal dashed lines on the MA plot
    * The data point colors will also be adjusted based on where the points lie
2. The user may adjust the significance threshold by entering the desired P-value, where the value will be calculated as -log10 and adjust the horizontal line on the volcano plot accordingly
3. The number of top hits begins at 10 genes in each fold change direction which can be seen labeled on the volcano and MA plots
    * Please note that increasing the number of the labeled genes may decrease the visibility of the gene names depending on over lapping
4. The user may select additional genes to be labeled on the volcano and MA plots by selecting from the available list or typing in a gene of interest
5. At the top of the volcano and MA plots there will be a description based on the sample types selected.
    * The user may select different sample types to compare in the 'RNAseq Parameters' tab of the side panel
6. Below the plot will appear a text box based on the data point that the users mouse may be hovering over
    * It will not appear if the mouse is not hovering a data point

![alt text](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/App_Demo_Pictures/GSEA_Shiny_7.png?raw=true)

1. Users may select two genes from the expression matrix to view in a scatter plot
2. These are from the base expression plot, so there is an option to log2+1 transform the data
3. A table will show below with expression data for those genes which can be downloaded as a .tsv file

![alt text](https://github.com/shawlab-moffitt/RShinyAnalysisGenerator/blob/main/App_Demo_Pictures/GSEA_Shiny_9.png?raw=true)

1. The user may select pathways to visualize with Enrichr for their samples
    * Below the upregulated visualization and table shows the downregulated pathway
3. The pathway analysis table can be downloaded as a .tsv file or it can be formated into a .gmt file for the selected pathway based on the genes listed in the table

# Quesions and Comments

Please email Alyssa Obermayer at alyssa.obermayer@moffitt.org if you have any further comments or questions in regards to the r Shiny Application.
