# Easy-App - Expression Analysis

# Introduction

With the rapid generation of high throughput datasets, the need for quick and reproducible data analysis tools has become paramount. Thus we developed DRPPM-EASY, an **E**xpression **A**nalysis App hosted in R **S**hin**Y**. The application minimizes the need for extensive computational experience, providing an easy plug-and-play process to allow the user to explore their expression data on a web-based platform. The application provides many of the typical exploratory analyses and visualization, such as unsupervised hierarchical clustering, Gene Set Enrichment (GSEA), and differential gene expression. Results can be downloaded as a comprehensive table or gmt gene set for further downstream analysis.

The user can access the easy app through https://shawlab-moffitt.shinyapps.io/drppm_easy_url_app/

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

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/App_Demo_Pictures/EASY_App_Overview.svg?raw=true)

### Reference Citation ###
Obermayer A, Dong L, Hu Q, Golden M, Noble JD, Rodriguez P, Robinson TJ, Teng M, Tan AC, Shaw TI*. DRPPM-EASY: A Web-Based Framework for Integrative Analysis of Multi-Omics Cancer Datasets. Biology (Basel). 2022 Feb 8;11(2):260. doi: 10.3390/biology11020260. PMID: 35205126 PMCID: PMC8869715

# Quesions and Comments

Please email Alyssa Obermayer at alyssa.obermayer@moffitt.org if you have any further comments or questions in regards to the R Shiny Application.



Copyright 2022 Moffitt Cancer Center Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
