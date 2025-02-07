# DRPPM Expression Analysis Shiny App

# Introduction

With the rapid generation of high throughput datasets, the need for quick and reproducible data analysis tools has become paramount. Thus we developed DRPPM-EASY, an **E**xpression **A**nalysis App hosted in R **S**hin**Y**. The application minimizes the need for extensive computational experience, providing an easy plug-and-play process to allow the user to explore their expression data on a web-based platform. The application provides many of the typical exploratory analyses and visualization, such as unsupervised hierarchical clustering, Gene Set Enrichment (GSEA), and differential gene expression. Results can be downloaded as a comprehensive table or gmt gene set for further downstream analysis.

The user can access the easy app through https://shawlab-moffitt.shinyapps.io/drppm_easy_url_app/

The ability to analyze omics data (such as genomic and proteomic data) is of significant interest. Thus we also generated DRPPM-EASY-Integration to facilitate the integration between multiple expression data set. [DRPPM-EASY-Integration](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration) allows users to compare expression data between two matrices. The user can use the platform to identify key features being shared or differentated between the two data sets. We have implemented a pipeline to perform reciprocal GSEA and gene set comparison, to examine the "connectivity" between two data matrices.

Finally, we developed an interface for exploration of complex hetergeneous data such as CCLE or CPTAC. [DRPPM-EASY-LargeProject](https://github.com/shawlab-moffitt/DRPPM-EASY-LargeProject) utilizes a sample selection tab which allows users to select expression and meta data and the [DRPPM-EASY-LargeProject-Integration](https://github.com/shawlab-moffitt/DRPPM-EASY-LargeProject-Integration) allows the user to compare there expression data with samples they choose to select from a large project dataset. As an example, we have implemented an instance of querying the Cancer Cell Line Encyclopedia (CCLE) data based on cancer type or lineage and sample type for analysis within the DRPPM-EASY app http://shawlab.science/shiny/DRPPM_EASY_LargeProject_CCLE/. Separate instance of the Large Project implementation is http://shawlab.science/shiny/DRPPM_EASY_LargeProject_LSCC_CPTAC/. Backup instance of these server are made available in 

The flow chart below gives a layout of the EASY app's family infrastructure which we will describe in further detail.

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/App_Demo_Pictures/EASY_App_Overview.svg?raw=true)

# Setup

* Users may download a zip file of the GitHub repository or cloan it to their local workspace using the command line.
* Please check to make sure all of the pacakge dependencies are installed.
* As long as file structure and naming isn't effected, the app can be ran without any additional coding and files cant be input from the interface.

## Advanced Setup

* Under the `User Data Input` section of the app in the first section of the app.R script, users may include startup files and parameters and the app will upload the data when starting up.
* If this app will be hosted on a server accessible to others, we have included a feature where you may specify a password in the script.
  * This password is not acsessible from the interface, only the app.R script.

# Requirements

* `R` - https://cran.r-project.org/src/base/R-4/

# R Dependencies

|  |  |  |  |  |
| --- | --- | --- | --- | --- |
| shiny_1.6.0 | shinythemes_1.2.0 | shinyjqui_0.4.0 | shinycssloaders_1.0.0 | tools_4.1.0 |
| dplyr_1.0.7 | tidyr_1.1.3 | readr_2.0.1 | tibble_3.1.3 | DT_0.18 |
| ggplot2_3.3.5 | plotly_4.9.4.1 | enrichplot_1.12.2 | pheatmap_1.0.12 | ggrepel_0.9.1 |
| enrichR_3.0 | limma_3.48.3 | clusterProfiler_4.0.5 | limma_3.48.3 | GSVA_1.40.1 |
| BiocManager_1.30.16 | reshape2_1.4.4 | ggpubr_0.4.0 |  |  |

# Required Files

## Required User Input Files

* **Expression Matrix:**
  * A tab or comma delimited table with gene symbols or other features located in the first column with subsequent columns consiting of the sample name as the header and expression data down the column.
  *  The current App expects lowly expressed genes filtered out and normalized data either to FPKM or TMM.
     * Larger files might inflict memory issues.
  * [Example](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/Example_Data/TCGA_CHOL_Expression_PatientID.txt)
* **Meta Data:**
  * A tab or comma delimited with two or more columns. First column consisting of sample names and following columns to be any additional meta or clinical data defining the sample.
  * [Example](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/Example_Data/TCGA_CHOL_Clinical_PatientID.txt)

## Required Files for App Setup - Provided

* **GeneSets:** 
  * This folder contains currated gene sets from public recourses such as MSigDB, LINCS L1000, Cell Marker, as well as gene sets annotated in our lab defining immune signatures and ER stress signatures.

# Shiny App Features

## Differential Gene expression Analysis

Differential gene expression analysis is performed between two groups of samples chosen by the user within the apps interface. The expression data for these groups is log transformed (log2 + 1) followed by a model matrix designed from the meta groups chosen. The matrix is then used for LIMMA linear modeling followed by Bayes statistics to rank the genes based on statistical significance. Figures and tables produce through these methods include hetmaps, volcano, MA, and box plots, as well as a gene expression scatter plot and pathway analysis with enrichR.

### Heatmaps

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/App_Demo_Pictures/EASY_Heatmap.png?raw=true)

1. The option the display a custom heatmap of selected genes and samples
2. Download sample cluster results based on cut tree and the clustering method chosen
3. Download most variable genes list as a ranked table or as a .gmt file for GSEA
   * Variance measure can be selected by user

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/App_Demo_Pictures/EASY_CustomHeatmap.png?raw=true)

1. User may select or type in genes they wish to view
   * This is initially filled with gene involved with tumor and immun cell interaction
2. Specific samples may also be selected to view
   * All samples are shows to start but within the box you may delete the name to un-select the sample

### Volcano Plot

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/App_Demo_Pictures/EASY_Volcano.png?raw=true)

1. User may select sample groups by meta type to compare in the LIMMA analysis
2. Log2FC and -log10(P.Value) may be adjusted to move the dashed lines on the plot
3. The number of top hits displayed may be adjusted and the user may select hits to annotate by selecting or writing in a gene name
   * Depending on how many genes being displayed, name may overlap
4. Hovering over a poin/gene on the plot will display a text box below showing the name of the gene and statistical values

### Box Plot

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/App_Demo_Pictures/EASY_BoxPlot.png?raw=true)

1. Users may search for and select a gene to view log2 transformed expression values as a box plot

### Gene Expression Scatter Plot

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/App_Demo_Pictures/EASY_GeneExpScatter.png?raw=true)

1. Users may select two genes they want to compare
2. The expression data show will be originally un-logged, so the user may choose to log2 transform
3. The user may download the table of expression values
   * The table will be of un-logged values

### EnrichR Pathway Analysis

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/App_Demo_Pictures/EASY_GeneExpScatter.png?raw=true)

1. Users may select the comparison groups for analysis
2. Adjusted P-Value and logFC cutoff may be adjusted based on the users preference of significance
   * Depending on the thresholds given an error might show if there are no genes that make this cutoff
3. The user may select a pathway from the EnrichR database
4. The up and down regulated pathway tables may be downloaded as they are or in GMT format

### Differentially Expressed Genes Top Table

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/App_Demo_Pictures/EASY_GeneExpScatter.png?raw=true)

1. The table displayed may be downloaded as a gmt file of one gene set based on the number of top hits input
   * At the bottom of the table there is a download button for the entire top table

## Gene Set Enrichment Analysis

Gene set enrichment analysis (GSEA) is performed through signal-to-noise calculations on the expression matrix. This scales the difference of means by standard deviation, producing a ranked list of genes. The GSEA function identifies the enrichment score of a gene set by taking the ranked genes list into account.

### Enrichment Plot

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/App_Demo_Pictures/EASY_EnrichPlot.png?raw=true)

1. The user may select the P-Value cutoff for the significance of the gene when running through the ranked list during GSEA
   * **Please note:** Depending on the designated P-value selected as well as the genes in the expression matrix, there will be some gene sets that may show errors. This is due to the genes within the gene set not being significantly enriched or available to produce any data or visualizations.
2. Users may select a gene set to view in the enrichment plot from MSigDB or another database that is loaded in
   * Initially the second tab has the Cell Marker database but that can be adjusted by the user
3. The user may upload their own gene set file as one of the formats described in the gene set bullet point of [reqired files](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY#required-files)
4. The Normalized Enrichement Score (NES) and P-value will be annotated above the plot
5. The leading edge gene list of the plot may be downloaded

### GSEA Heatmap

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/App_Demo_Pictures/EASY_GSEAheatmap.png?raw=true)

1. The user may adjust figure parameters
2. The user must choose a gene set to display
3. The user may also upload there own, similar to the enrichment plot

### Enriched Signatures Table

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/App_Demo_Pictures/EASY_EnrichTable.png?raw=true)

1. If multiple tables loaded into the app, the user may choose which to view based on file name
2. The user may generate their own enriched signatures table with a gmt or gene set file of their choice

### Single Sample GSEA Box Plot

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/App_Demo_Pictures/EASY_ssGSEAboxplot.png?raw=true)

1. A single sample GSEA (ssGSEA) method may be chosen and a stat compare method to view in the boxplot
2. The user may choose a gene set to view
3. If the user chooses to use their own gene set they will have to upload here and select the button to generate an RData list
4. The ssGSEA score file may be downloaded here
   * This table may be used in the [DRPPM-EASY-Integration App](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration) wher you can compare two of these types of files to each other, or you can compare one of the files to an expression matrix to generate a correlation ranking

### Reference Citation ###
Obermayer A, Dong L, Hu Q, Golden M, Noble JD, Rodriguez P, Robinson TJ, Teng M, Tan AC, Shaw TI*. DRPPM-EASY: A Web-Based Framework for Integrative Analysis of Multi-Omics Cancer Datasets. Biology (Basel). 2022 Feb 8;11(2):260. doi: 10.3390/biology11020260. PMID: 35205126 PMCID: PMC8869715

# Quesions and Comments

Please email Alyssa Obermayer at alyssa.obermayer@moffitt.org if you have any further comments or questions in regards to the R Shiny Application.
