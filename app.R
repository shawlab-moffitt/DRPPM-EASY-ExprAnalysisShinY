type_id <- paste0("v2.0.20250214")

# User Data Input --------------------------------------------------------------
# Project Name
ProjectName <- ''
# Expression Matrix
expression_file <- ''
# Meta data
meta_file <- ''
#If human: set TRUE
#If mouse: set FALSE
human <- TRUE

## Advanced User Input ---------------------------------------------------------
set.seed(10242022)
# Subsetting Feature
Subsetting_Feature <- ''
# Starting Feature
Feature_Selected <- ''
# Does user want password Protection?
Password_Protected <- FALSE
PasswordSet <- 'password'


## Provided Input --------------------------------------------------------------
## User make sure paths are correct
ExampleExpr_File <- "Example_Data/TCGA_CHOL_Expression_PatientID.txt"
ExampleClin_File <- "Example_Data/TCGA_CHOL_Clinical_PatientID.txt"
GeneSetHS_File <- "GeneSets/GeneSet_List_HS_v6.RData"
GeneSetTableHS_File <- "GeneSets/GeneSet_Categories_HS.txt"
GeneSetMM_File <- "GeneSets/GeneSet_List_MM.RData"
GeneSetTableMM_File <- "GeneSets/GeneSet_Categories_MM.txt"
MM_HS_Conversion_File <- "GeneSets/hs_mm_conversion_20240701.txt"


#increase file upload size
options(shiny.maxRequestSize=5000*1024^2)
#options("recount3_organism_human_project_homes_URL_http://duffel.rail.bio/recount3" = c("data_sources/sra", "data_sources/gtex", "data_sources/tcga"))

# Load Libraries ---------------------------------------------------------------
library("shiny")
library("shinythemes")
library("shinyjqui")
library("pheatmap")
library("BiocManager")
library("RColorBrewer")
library("dplyr")
library("DT")
library("ggplot2")
library("ggpubr")
library("reshape2")
library("tibble")
library("viridis")
library("scales")
library("plotly")
library("readr")
library("enrichR")
library("ggrepel")
library("tidyr")
library("tools")
library("shinycssloaders")
library("Hmisc")
library("stringr")
library("glue")
library("data.table")
library("clusterProfiler")
library("GSVA")
library("limma")
library("enrichplot")
library("ComplexHeatmap")
library("GEOquery")
library("edgeR")
library("DESeq2")
library("svglite")

if ("recount" %in% rownames(installed.packages()) && "recount3" %in% rownames(installed.packages())) {
  library("recount")
  library("recount3")
  options("recount3_organism_human_project_homes_URL_http://duffel.rail.bio/recount3" = c("data_sources/sra", "data_sources/gtex", "data_sources/tcga"))
  RecountAvail <- TRUE
}

#'biomaRt', 'sva', 'ConsensusTME', 'quantiseqr','singscore'
immudecon_check <- "immunedeconv" %in% rownames(installed.packages())
if (immudecon_check == TRUE) {
  library(immunedeconv)
}


date_to_gene <- function(vec,mouse = FALSE) {
  if (is.null(vec)) exit("Please provide vector argument")
  if (!mouse) {
    vec <- gsub('(\\d+)-(Mar)','MARCH\\1',vec)
    vec <- gsub('(\\d+)-(Sep)','SEPT\\1',vec)
    vec <- gsub('(\\d+)-(Dec)','DEC\\1',vec)
  } else {
    vec <- gsub('(\\d+)-(Mar)','March\\1',vec)
    vec <- gsub('(\\d+)-(Sep)','Sept\\1',vec)
    vec <- gsub('(\\d+)-(Dec)','Dec\\1',vec)
  }
  return(vec)
}


## Gene Set data ---------------------------------------------------------------

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
CTKgenes <- c("IL2","IL12A","IL12B","IL17A","IFNA1","IFNB1","IFNG","IFNGR","CD11b",
              "ITGAM","CD33","ENTPD1","ICOSLG","CD275","CD278","TNFSF9","TNFRSF9",
              "CD40","CD40LG","CD70","CD27","TNFSF18","TNFRSF18","TNFSF14","TNFRSF14",
              "TNFSF4","TNFRSF4","HLA-A","CD3","CEACAM1","CD80","CD86","CTLA4","CD276",
              "VTCN1","PVR","CD226","TIGIT","CD96","LGALS3","LGALS3BP","LGALS9","LGALS9C",
              "HAVCR2","HHLA2","TMIGD2","CD274","PDCD1LG2","PDCD1","VSIR")
CTKgenes_hs <- c("IL2","IL12A","IL12B","IL17A","IFNA1","IFNB1","IFNG","IFNGR","CD11b",
                 "ITGAM","CD33","ENTPD1","ICOSLG","CD275","CD278","TNFSF9","TNFRSF9",
                 "CD40","CD40LG","CD70","CD27","TNFSF18","TNFRSF18","TNFSF14","TNFRSF14",
                 "TNFSF4","TNFRSF4","HLA-A","CD3","CEACAM1","CD80","CD86","CTLA4","CD276",
                 "VTCN1","PVR","CD226","TIGIT","CD96","LGALS3","LGALS3BP","LGALS9","LGALS9C",
                 "HAVCR2","HHLA2","TMIGD2","CD274","PDCD1LG2","PDCD1","VSIR")

CTKgenes_mm <- c("Il2","Il12a","Il12b","Il17a","Ifna13","Ifnb1","Ifng","Ifngr1","Cd11b","Itgam",
                 "Cd33","Entpd1","Icosl","Icos","Tnfsf9","Tnfrsf9","Cd40","Cd40lg","Cd70","Cd27",
                 "Tnfsf18","Tnfrsf18","Tnfsf14","Tnfrsf14","Tnfsf4","Tnfrsf4","H2-K1","CD3G",
                 "Ceacam1","Cd80","Cd86","Ctla4","Cd276","Vtcn1","Pvr","Cd226","Tigit","Cd96","Lgals3",
                 "Lgals3bp","Lgals9","Lgals9c","Havcr2","Hhla2","Cd274","Pdcd1lg2","Pdcd1","Vsir")

if (!file.exists(expression_file)) {
  FileProvided <- FALSE
} else { FileProvided <- TRUE }

if (!isTruthy(ProjectName)) {
  ProjectNameHeader <- paste("{ Expression Analysis }")
} else {
  ProjectNameHeader <- paste("{",ProjectName,"Expression Analysis }")
}

if (!isTruthy(Subsetting_Feature)) {
  Subsetting_Feature <- "Select All Samples"
}
if (!isTruthy(Feature_Selected)) {
  Feature_Selected <- 1
}

mouse <- ifelse(!as.logical(human),TRUE,FALSE)


# Functions --------------------------------------------------------------------
## Function to identify numeric or character columns
#column_type_check <- function(data, type = "numeric", max_na_prop = 0.1) {
#  # Check if the 'data' argument is a data frame
#  if (!is.data.frame(data)) {
#    stop("Input must be a data frame.")
#  }
#
#  # Initialize a vector to store column names
#  num_cols <- c()
#  chr_cols <- c()
#
#  for (col in names(data)) {
#    max_na <- (length(data[,col][!is.na(data[,col])]) - length(data[,col][is.infinite(data[,col])])) * max_na_prop
#    if (!is.logical(data[[col]])) {
#      if (!all(is.na(as.numeric(data[[col]]))) | sum(is.na(as.numeric(data[[col]]))) <= max_na) {
#        num_cols <- c(num_cols,col)
#      } else {
#        chr_cols <- c(chr_cols,col)
#      }
#    } else {
#      chr_cols <- c(chr_cols,col)
#    }
#  }
#
#  if (tolower(type) == "numeric") {
#    return(num_cols)
#  } else if (tolower(type) == "character") {
#    return(chr_cols)
#  } else {
#    stop("Invalid 'type' argument. Use 'numeric' or 'character'.")
#  }
#
#}

get_feat_cols <- function(data, keep_num = TRUE) {
  if (keep_num) {
    data <- data[sapply(data, function(x) length(unique(x))>1)]
  } else {
    data <- data[sapply(data, function(x) length(unique(x))>1)]
    data <- data[sapply(data, function(x) length(unique(x))<length(x))]
  }
  return(colnames(data))
}

# Function to identify numeric or character columns
column_type_check <- function(data, type = "numeric", max_na_prop = 0.1) {
  # Check if the 'data' argument is a data frame
  if (!is.data.frame(data)) {
    stop("Input must be a data frame.")
  }
  
  # Initialize a vector to store column names
  num_cols <- c()
  chr_cols <- c()
  
  for (col in names(data)) {
    max_na <- (length(data[,col][!is.na(data[,col])]) - length(data[,col][is.infinite(data[,col])])) * max_na_prop
    if (!is.logical(data[[col]])) {
      if (!all(is.na(as.numeric(data[[col]]))) | sum(is.na(as.numeric(data[[col]]))) <= max_na) {
        num_cols <- c(num_cols,col)
      } else {
        chr_cols <- c(chr_cols,col)
      }
    } else {
      chr_cols <- c(chr_cols,col)
    }
  }
  
  if (tolower(type) == "numeric") {
    #return(num_cols)
    return(colnames(data))
  } else if (tolower(type) == "character") {
    #return(chr_cols)
    return(colnames(data))
  } else {
    stop("Invalid 'type' argument. Use 'numeric' or 'character'.")
  }
  
}

cv <- function(x){
  (sd(x)/mean(x))*100
}

ExprFilter2 <- function(vec,criteria,proportion) {
  Samp2meet <- length(vec)*proportion
  meet <- sum(vec>criteria)
  if (meet > Samp2meet) {
    return(TRUE)
  } else { return(FALSE) }
}

boxopt <- c("none","wilcox.test","t.test","kruskal.test","anova")

detect_species <- function(genes) {
  # check the capitalization pattern
  is_human_gene <- function(gene) {
    return(grepl("^[A-Z]+$", gene))
  }
  is_mouse_gene <- function(gene) {
    return(grepl("^[A-Z][a-z]*$", gene))
  }
  # Count the number of human and mouse gene patterns
  human_count <- sum(sapply(genes, is_human_gene))
  mouse_count <- sum(sapply(genes, is_mouse_gene))
  # Determine the majority match
  if (human_count > mouse_count) {
    return("human")
  } else if (mouse_count > human_count) {
    return("mouse")
  } else {
    return("undetermined") # If counts are equal or if there are no matches
  }
}


# Password Table ---------------------------------------------------------------
# user database for logins
if (Password_Protected) {
  user_base <- tibble::tibble(
    user = "user",
    password = PasswordSet,
    permissions = "admin",
    name = "User"
  )
}

# UI Tabs ----------------------------------------------------------------------

## Login Tab -------------------------------------------------------------------
login_tab <- tabPanel(
  title = icon("lock"),
  value = "login",
  loginUI("login")
)

## About Tab -------------------------------------------------------------------
About_tab <- tabPanel("About",
                      fluidPage(
                        mainPanel(
                          h3("Introduction"),
                          p("With the rapid generation of large datasets as a result of advancing next-generation sequecing (NGS) technologies,
                            the need for quick and reproducible data analysis tools has become paramount. The ability to analyze both genomic
                            and proteomic data sets can allow for the detection of compelling genes and features that may be of interest for further analysis.
                            The General Expression Analysis App hosted in R Shiny does not require an extensive computation background to run these analyses.
                            With user input of expression and meta data, along with gene set we have sourced and provided, the user may produce useful anaylsis
                            and visualizations on a web-based interface within minutes. This method of reproducible data analysis generates a number of
                            visualization to view Gene Set Enrichment (GSEA) and differential gene expression analysis, as well as the ability to download various
                            tables and gene sets that are produced by the App for further use. Additional apps are being developed for the EASY family.
                            DRPPM-EASY-Integraction allows users to further analyze data obtained from the main DRPPM-EASY app and compare expression data between two matrices.
                            DRPPM-EASY-CCLE integrates a sample selection tab which allows users to select expression and meta data from the Cancer Cell Line Encyclopedia (CCLE)
                            based on cancer type or lineage and sample type for analysis within the DRPPM-EASY app. The flow chart below gives a layout of the
                            EASY app's family infrastructure which we will describe in further detail."),
                          h3("Unsupervised Clustering Methods"),
                          p("Unsupervised clustering is performed by calculating the top variable genes through MAD, CV, or VAR as the variance measure.
                            The choice of variance measure and number to top genes/probes can be determined by the user, as well as the number of clusters with k-means.
                            The identified most variable genes can be downloaded as a table for the user. The “complete” method was chosen as the default
                            algorithm for clustering, but other methods such as Ward, average, and centroid could be chosen based on the ‘hclust()’ function in R."),
                          h3("Differential Gene Expression Methods"),
                          p("Differential gene expression analysis is performed on the expression data between Limma: Two groups defined by the provided metadata file.
                            The samples chosen are log-transformed (log2 + 1). A model matrix is then designed for LIMMA linear modeling followed by empirical
                            Bayes statistics to moderate gene variance and modeling of the global characteristic of the data. Of note, the current pipeline expects
                            the matrix input quantification are pre-normalized, such as in CPM, TMM, or FPKM."),
                          h3("Gene Set Enrichment Analysis (GSEA) Methods"),
                          p("Gene Set Enrichment Analysis (GSEA) is performed through signal-to-noise calculations on the expression matrix.
                            This calculation requires at least two types of sample groups and at least three samples for each grouping.
                            The signal-to-noise calculation scales the difference of means by standard deviation and generates a ranked list of genes,
                            effectively capturing the differences between phenotypes. The GSEA identifies the enrichment score (ES) for a gene set,
                            which indicates the extent that the gene set is overrepresented at the beginning or end of the list of ranked genes.
                            The enrichment score is calculated by walking down the ranked gene list and increasing the running-sum statistic when a gene
                            is in the gene set and decreasing when it is not. The maximum deviation from zero that is encountered through walking the list is the ES,
                            where a positive ES indicates the gene set is enriched at the top of the ranked list and a negative ES indicates the gene set is enriched
                            at the bottom of the ranked list. A leading-edge gene list is provided displaying a subset of the genes in the gene set that contribute
                            the most to the ES."),
                          h4("Gene Set Sources"),
                          p("The Molecular Signatures Database (MSigDB) is a collection of over 32,000 annotated gene sets divided into 9 major collections
                            as well as various sub-collections [PMID: 16199517, PMID: 12808457]. The MSigDB  gene sets were downloaded with the msigdbr package
                            and processed through R using Tidyverse packages to combine the gene sets into a data frame."),
                          p("The Cell Marker gene sets were derived from a comprehensive list of cell markers from cell types in various tissues obtained
                            through manually curating over 100,000 published papers [PMID: 30289549]. This data was obtained as an extensive data frame and
                            parsed to make gene sets with the tidyverse packages."),
                          p("The Library of Integrated Network-based Cellular Signatures (LINCS) L1000 gene sets were derived from the Connectivity Map (CMAP)
                            projects data which treated 4 human cell lines with ~1300 drugs followed by profiling genome-wide mRNA expression following small
                            molecule perturbation [PMID: 24906883].")
                        )
                      )
)

## Data Input Tab --------------------------------------------------------------
DataInput_tab <- tabPanel("Data Input",
                          ### Sidebar ------------------------------------------
                          fluidPage(
                            sidebarPanel(
                              width = 3,
                              id = "DataInputPanel",
                              p(),
                              textInput("UserProjectName","Project Name:", value = "Expression Analysis"),
                              uiOutput("rendExprFileInput"),
                              fluidRow(
                                column(6, style = 'margin-top:-20px;',
                                       if (immudecon_check == TRUE) {
                                         checkboxGroupInput("RawCountQuantNorm",label = NULL,
                                                            choices = c("Input data is log-transformed",
                                                                        "Normalize Raw Counts",
                                                                        "Quantile Normalization",
                                                                        "Filter Matrix",
                                                                        "Perform Immune Deconvolution"))
                                       } else {
                                         checkboxGroupInput("RawCountQuantNorm",label = NULL,
                                                            choices = c("Input data is log-transformed",
                                                                        "Normalize Raw Counts",
                                                                        "Quantile Normalization",
                                                                        "Filter Matrix"))
                                       }
                                ),
                                column(6, style = 'margin-top:-20px;',
                                       uiOutput("rendRawCountNorm"),
                                       uiOutput("rendDESeqDesignCol"),
                                       uiOutput("rendDESeqDesignColRef")
                                )
                              ),
                              conditionalPanel(condition = "input.RawCountQuantNorm.includes('Perform Immune Deconvolution')",
                                               radioButtons("ImmDeconvAnalysis","Perform Expression Analysis On:",choices = c("Gene Expression","Immune Deconvolution"),
                                                            selected = "Gene Expression", inline = T)
                                               ),
                              #if (immudecon_check == TRUE) {
                              #  radioButtons("ImmDeconvAnalysis","Perform Expression Analysis On:",choices = c("Gene Expression","Immune Deconvolution"),
                              #               selected = "Gene Expression", inline = T)
                              #},
                              conditionalPanel("input.RawCountQuantNorm.includes('Filter Matrix')",
                                               fluidRow(
                                                 column(6, style = 'margin-bottom:-15px;',
                                                        numericInput("FilterNum","Filter value:",value = 1, min = 0)
                                                 ),
                                                 column(6, style = 'margin-bottom:-15px;',
                                                        numericInput("FilterProp","Proportion of samples (%):",value = 10, max = 100, min = 0)
                                                 )
                                               )
                              ),
                              hr(),
                              uiOutput("rendClinFileInput"),
                              div(uiOutput("rendMetaColOfInterest"), style = "margin-top:-20px"),
                              fluidRow(
                                column(12, style = 'margin-top:-10px;margin-bottom:-15px;',
                                       actionButton("UseExpData","Load Example Data"),
                                       tags$a(href="http://shawlab.science/shiny/PATH_SURVEYOR_ExampleData/PATH_SURVEYOR_App/", "Download example data", target='_blank'),
                                )
                              ),
                              hr(),
                              tabsetPanel(
                                tabPanel("Subest Data",
                                         p(),
                                         selectizeInput("SubsetCol","Subset Samples By:",choices = NULL,selected = NULL),
                                         #selectizeInput("SubsetCrit","Sample Criteria:",choices = NULL,selected = 1),
                                         #uiOutput("rendSubsetCol"),
                                         uiOutput("rendSubsetCrit")
                                ),
                                tabPanel("Split Columns",
                                         p(),
                                         fluidRow(
                                           column(8, style = 'padding-right:2px',
                                                  selectizeInput("ColToSplit","Column to split:",choices = NULL, selected = "")
                                           ),
                                           column(4, style = 'padding-left:2px',
                                                  textInput("SplitDelim","Split character:")
                                           )
                                         ),
                                         textInput("SplitNewColNames","New column names (comma delim)",
                                                   placeholder = "NewColumn1,NewColumn2,NewColumn3"),
                                         actionButton("SplitColumn","Click to split column")
                                ),
                                tabPanel("Merge Columns",
                                         p(),
                                         fluidRow(
                                           column(8, style = 'padding-right:2px',
                                                  selectizeInput("ColsToMerge","Columns to merge:",choices = NULL, selected = "", multiple = T)
                                           ),
                                           column(4, style = 'padding-left:2px',
                                                  textInput("MergeDelim","Merge character:")
                                           )
                                         ),
                                         textInput("MergeNewColName","New column name:",
                                                   placeholder = "New_Merged_Column"),
                                         actionButton("MergeColumns","Click to merge columns")
                                )
                              )
                            ),
                            ### Main Panel -------------------------------------
                            mainPanel(
                              verbatimTextOutput("FileCheckAlerts"),
                              fluidRow(
                                column(2, style = 'margin-top:15px;padding-right:2px',
                                       uiOutput("renddownload_expr")
                                ),
                                column(2, style = 'margin-top:15px;padding-left:2px;padding-right:2px',
                                       uiOutput("renddownload_meta")
                                ),
                                column(2, style = 'margin-top:15px;padding-left:2px;padding-right:2px',
                                       uiOutput("renddownload_notes")
                                )
                              ),
                              uiOutput("rendExprFilePrevHeader"),
                              uiOutput("rendExprHead"),
                              div(DT::dataTableOutput("ExprFile_Preview"), style = "font-size:10px"),
                              uiOutput("rendClinFilePrevHeader"),
                              uiOutput("rendClinHead"),
                              div(DT::dataTableOutput("ClinFile_Preview"), style = "font-size:10px")
                            )
                          ),
                          tagList( #nolint
                            tags$head(
                              tags$style(
                                HTML(
                                  "
          .info_box {
            width: auto;
            height: auto;
            color: #000000;
            background-color: #f5f5f5;
            padding: 3px 8px;
            font-size: 12px;
            z-index : 9999;
          }
          ",
          glue::glue(
            "
            #{'AppVersion'} {{
              position: {'fixed'};
              top: 0;
              right: 0;
            }}
            ")
                                )
                              )
                            ),
          div(id = "AppVersion", class = "info_box", type_id)
                          )
)

## Data Exploration Tab -------------------------------------------------------------------
Data_Exploration_tab <- tabPanel("Data Exploration",
                                 fluidPage(
                                   title = "Data Exploration",
                                   sidebarLayout(
                                     sidebarPanel(
                                       width = 3,
                                       tags$head(
                                         tags$style(HTML(
                                           '.selectize-input {
                                           max-height: 102px;
                                           overflow-y: auto;
                                           }'
                                         )
                                         )
                                       ),
                                       ### Sidebar --------------------------------------------
                                       
                                       #uiOutput("rendSubsetCol"),
                                       #uiOutput("rendSubsetCrit"),
                                       conditionalPanel(condition = "input.dataset == '2' | input.dataset == '4'",
                                                        uiOutput("rendGroupCol")
                                       ),
                                       conditionalPanel(condition = "input.dataset == '1'",
                                                        conditionalPanel(condition = "input.datasetheat != '5'",
                                                                         uiOutput("rendGroupColMulti")
                                                        ),
                                                        conditionalPanel(condition = "input.datasetheat == '5'",
                                                                         uiOutput("rendGroupColAvgHeat")
                                                        ),
                                       ),
                                       #### Heatmaps
                                       
                                       ##### MVG Heatmap
                                       
                                       conditionalPanel(condition = "input.dataset == '1'",
                                                        conditionalPanel(condition = "input.datasetheat == '1'",
                                                                         tabsetPanel(
                                                                           tabPanel("Data Input",
                                                                                    p(),
                                                                                    numericInput("NumFeatures", step = 1, label = "Number of Genes", value = 100),
                                                                                    hr(),
                                                                                    h4("Clustering Parameters"),
                                                                                    selectInput("VarianceMeasure", "Select Variance Measure",
                                                                                                choices = c("MAD","CV","VAR")),
                                                                                    fluidRow(
                                                                                      column(6,
                                                                                             checkboxInput("clustrowsMVG","Cluster Rows", value = T)
                                                                                      ),
                                                                                      column(6,
                                                                                             checkboxInput("clustcolsMVG","Cluster Columns", value = T)
                                                                                      )
                                                                                    ),
                                                                                    
                                                                                    uiOutput("ClusterMethodMVG"),
                                                                                    #numericInput("NumClusters", step = 1, label = "Number of Clusters (Cut Tree with ~k)", value = 2),
                                                                                    h4("Download Cluster Result:"),
                                                                                    downloadButton("downloadClusters", "Download .tsv"),
                                                                                    h4("Download Most Variable Genes List:"),
                                                                                    downloadButton("MVGdownload", "Download MVG .tsv"),
                                                                                    downloadButton("MVGdownloadgmt", "Download MVG .gmt")
                                                                           ),
                                                                           tabPanel("Figure Settings",
                                                                                    p(),
                                                                                    fluidRow(
                                                                                      column(6,
                                                                                             checkboxInput("ShowColNames1","Show Column Names", value = T),
                                                                                             numericInput("heatmapFont2.c", "Column Font Size:",
                                                                                                          min = 5, max = 75,
                                                                                                          value = 12, step = 1)
                                                                                      ),
                                                                                      column(6,
                                                                                             checkboxInput("ShowRowNames1","Show Row Names", value = T),
                                                                                             numericInput("heatmapFont2.r", "Row Font Size:",
                                                                                                          min = 5, max = 75,
                                                                                                          value = 9, step = 1)
                                                                                      )
                                                                                    ),
                                                                                    fluidRow(
                                                                                      column(6,
                                                                                             numericInput("heatmapHeight1","Download Height (in)",value = 20)
                                                                                      ),
                                                                                      column(6,
                                                                                             numericInput("heatmapWidth1","Download Width (in)",value = 15)
                                                                                      )
                                                                                    )
                                                                           )
                                                                         )
                                                                         
                                                        )
                                       ),
                                       
                                       ##### Custom Heatmap
                                       
                                       conditionalPanel(condition = "input.dataset == '1'",
                                                        conditionalPanel(condition = "input.datasetheat == '2'",
                                                                         tabsetPanel(
                                                                           tabPanel("Data Input",
                                                                                    p(),
                                                                                    selectInput("heatmapGeneSelec","Gene Selection:",
                                                                                                choices = NULL,
                                                                                                multiple = T, selected = 1),
                                                                                    textInput("userheatgenes", "Text Input of Gene List (space delimited):", value = ""),
                                                                                    selectizeInput("userheatsamp2", "Samples Selection:",
                                                                                                   choices = NULL, multiple = T, selected = 1),
                                                                                    h4("Clustering Parameters"),
                                                                                    fluidRow(
                                                                                      column(6,
                                                                                             checkboxInput("ClusRowOptCust","Cluster Rows", value = FALSE)
                                                                                      ),
                                                                                      column(6,
                                                                                             checkboxInput("ClusColOptCust","Cluster Columns", value = FALSE)
                                                                                      )
                                                                                    ),
                                                                                    uiOutput("rendClustMethodsCust"),
                                                                           ),
                                                                           tabPanel("Figure Settings",
                                                                                    p(),
                                                                                    fluidRow(
                                                                                      column(6,
                                                                                             checkboxInput("ShowRowNames2","Show Row Names", value = T),
                                                                                             numericInput("heatmapFont3.r", "Row Font Size:",
                                                                                                          min = 5, max = 75,
                                                                                                          value = 12, step = 1)
                                                                                      ),
                                                                                      column(6,
                                                                                             checkboxInput("ShowColNames2","Show Column Names", value = T),
                                                                                             numericInput("heatmapFont3.c", "Column Font Size:",
                                                                                                          min = 5, max = 75,
                                                                                                          value = 12, step = 1)
                                                                                      )
                                                                                    ),
                                                                                    fluidRow(
                                                                                      column(6,
                                                                                             numericInput("heatmapHeight2","Download Height (in)",value = 20)
                                                                                      ),
                                                                                      column(6,
                                                                                             numericInput("heatmapWidth2","Download Width (in)",value = 15)
                                                                                      )
                                                                                    )
                                                                           )
                                                                         )
                                                        )
                                       ),
                                       
                                       ##### DEG Heatmap
                                       
                                       conditionalPanel(condition = "input.dataset == '1'",
                                                        conditionalPanel(condition = "input.datasetheat == '3'",
                                                                         tabsetPanel(
                                                                           tabPanel("Data Input",
                                                                                    p(),
                                                                                    uiOutput("rendDEGcolHeat"),
                                                                                    uiOutput("rendDESeqDesignColRef_heat"),
                                                                                    radioButtons("volcanoCompChoice4","Comparison Method:",
                                                                                                 choices = c("Limma: Two groups","Limma: One group"), inline = T),
                                                                                    conditionalPanel(condition = "input.volcanoCompChoice4 == 'Limma: Two groups'",
                                                                                                     fluidRow(
                                                                                                       column(6,
                                                                                                              uiOutput("rendcomparisonA2_h")
                                                                                                       ),
                                                                                                       column (6,
                                                                                                               uiOutput("rendcomparisonB2_h"))
                                                                                                     )
                                                                                    ),
                                                                                    conditionalPanel(condition = "input.volcanoCompChoice4 == 'Limma: One group'",
                                                                                                     uiOutput("rendcomparisonA2_h_one")
                                                                                    ),
                                                                                    #conditionalPanel(condition = "input.RawCountQuantNorm.includes('Normalize Raw Counts') & input.volcanoCompChoice4 == 'DESeq2'",
                                                                                    #                 uiOutput("rendDESeqDesignCol_heat"),
                                                                                    #                 uiOutput("rendDESeqDesignColRef_heat")
                                                                                    #),
                                                                                    #fluidRow(
                                                                                    #  column(6,
                                                                                    #         uiOutput("rendcomparisonA2_h")
                                                                                    #  ),
                                                                                    #  column(6,
                                                                                    #         uiOutput("rendcomparisonB2_h")
                                                                                    #  )
                                                                                    #),
                                                                                    selectizeInput("DEGCovarSelectHeat","Select Covariates:", choices = NULL, selected = 1, multiple = T),
                                                                                    fluidRow(
                                                                                      column(6,
                                                                                             numericInput("fc_cutoff_h", "LogFC Threshold",
                                                                                                          min = 0, max = 5, step = 0.1, value = 0)
                                                                                      ),
                                                                                      column(6,
                                                                                             numericInput("p_cutoff_h", "P.Value Threshold:",
                                                                                                          min = 0, max = 10, step = 0.01, value = 0.05)
                                                                                      )
                                                                                    ),
                                                                                    verbatimTextOutput("GenesAboveCutoff1"),
                                                                                    numericInput("top_x_h", "Number of Top Hits to Show on Heatmap:",
                                                                                                 value = 100, min = 0),
                                                                                    hr(),
                                                                                    h4("Clustering Parameters"),
                                                                                    fluidRow(
                                                                                      column(6,
                                                                                             checkboxInput("ClusRowOptDEG","Cluster Rows", value = TRUE)
                                                                                      ),
                                                                                      column(6,
                                                                                             checkboxInput("ClusColOptDEG","Cluster Columns", value = TRUE)
                                                                                      )
                                                                                    ),
                                                                                    uiOutput("rendClustMethodsDEG")
                                                                           ),
                                                                           tabPanel("Figure Settings",
                                                                                    p(),
                                                                                    fluidRow(
                                                                                      column(6,
                                                                                             checkboxInput("ShowRowNames3","Show Row Names", value = T),
                                                                                             numericInput("heatmapFont3.r.deg", "Row Font Size:",
                                                                                                          min = 5, max = 75,
                                                                                                          value = 12, step = 1)
                                                                                      ),
                                                                                      column(6,
                                                                                             checkboxInput("ShowColNames3","Show Column Names", value = T),
                                                                                             numericInput("heatmapFont3.c.deg", "Column Font Size:",
                                                                                                          min = 5, max = 75,
                                                                                                          value = 12, step = 1)
                                                                                      )
                                                                                    ),
                                                                                    fluidRow(
                                                                                      column(6,
                                                                                             numericInput("heatmapHeight3","Download Height (in)",value = 20)
                                                                                      ),
                                                                                      column(6,
                                                                                             numericInput("heatmapWidth3","Download Width (in)",value = 15)
                                                                                      )
                                                                                    )
                                                                           )
                                                                         )
                                                                         
                                                        )
                                       ),
                                       
                                       ##### Average Expr Heatmap
                                       
                                       conditionalPanel(condition = "input.dataset == '1'",
                                                        conditionalPanel(condition = "input.datasetheat == '5'",
                                                                         tabsetPanel(
                                                                           tabPanel("Data Input",
                                                                                    p(),
                                                                                    uiOutput("rendSampCondSelectionCust"),
                                                                                    selectizeInput("avgheatmapGeneSelec","Gene Selection:",
                                                                                                   choices = NULL,
                                                                                                   multiple = T, selected = 1),
                                                                                    textInput("avguserheatgenes", "Text Input of Gene List (space delimited):", value = ""),
                                                                                    h4("Clustering Parameters"),
                                                                                    fluidRow(
                                                                                      column(6,
                                                                                             checkboxInput("AvgClusRowOptCust","Cluster Rows", value = FALSE)
                                                                                      ),
                                                                                      column(6,
                                                                                             checkboxInput("AvgClusColOptCust","Cluster Columns", value = FALSE)
                                                                                      )
                                                                                    ),
                                                                                    uiOutput("rendAvgClustMethodsCust")
                                                                           ),
                                                                           tabPanel("Figure Settings",
                                                                                    p(),
                                                                                    selectInput("ColorPalette5", "Select Color Palette:",
                                                                                                choices = c("Red/Blue" = "original",
                                                                                                            "OmniBlueRed" = "OmniBlueRed",
                                                                                                            "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                                                                            "Green/Black/Red" = "GreenBlackRed",
                                                                                                            "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                                                                            "Viridis" = "Viridis","Plasma" = "Plasma",
                                                                                                            "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")),
                                                                                    fluidRow(
                                                                                      column(6,
                                                                                             checkboxInput("ShowRowNames5","Show Row Names", value = T),
                                                                                             numericInput("heatmapFont3.r.deg.avg.cust", "Row Font Size:",
                                                                                                          min = 5, max = 75,
                                                                                                          value = 12, step = 1)
                                                                                      ),
                                                                                      column(6,
                                                                                             checkboxInput("ShowColNames5","Show Row Names", value = T),
                                                                                             numericInput("heatmapFont3.c.deg.avg.cust", "Column Font Size:",
                                                                                                          min = 5, max = 75,
                                                                                                          value = 12, step = 1)
                                                                                      )
                                                                                    )
                                                                           )
                                                                         )
                                                                         
                                                        )
                                       ),
                                       
                                       #### Gene scatter
                                       
                                       conditionalPanel(condition = "input.dataset == '2'",
                                                        tabsetPanel(
                                                          tabPanel("Data Input",
                                                                   p(),
                                                                   selectizeInput("scatterG1","Select Gene 1:",
                                                                                  choices = NULL,
                                                                                  multiple = F, selected = 1),
                                                                   selectizeInput("scatterG2","Select Gene 2:",
                                                                                  choices = NULL,
                                                                                  multiple = F, selected = 1),
                                                                   checkboxInput("logask","Log2 Transform Expression Data")
                                                          ),
                                                          tabPanel("Figure settings",
                                                                   p(),
                                                                   fluidRow(
                                                                     column(6,
                                                                            numericInput("GeneScatterTitleSize","Title Font Size",
                                                                                         value = 14, step = 1)
                                                                     ),
                                                                     column(6,
                                                                            numericInput("GeneScatterAxisSize","Axis Lable Font Size",
                                                                                         value = 12, step = 1)
                                                                     )
                                                                   )
                                                          )
                                                        ),
                                                        
                                       ),
                                       
                                       #### Boxplot
                                       
                                       conditionalPanel(condition = "input.dataset == '3'",
                                                        tabsetPanel(
                                                          tabPanel("Data Input",
                                                                   p(),
                                                                   uiOutput("rendBPgroupCriteria"),
                                                                   uiOutput("rendBPgroupSelection"),
                                                                   #selectInput("BPFeatureCategory","Select Feature Type:",
                                                                   #c("Matrix Features","Meta Features","Immune Deconvolution Features")),
                                                                   #            c("Matrix Features","Meta Features")),
                                                                   radioButtons("BPFeatureCategory","Select Feature Type:",
                                                                                choices = c("Matrix Features","Meta Features"), inline = T),
                                                                   uiOutput("rendImmuneDeconvMethods"),
                                                                   fluidRow(
                                                                     column(9, style = 'padding-right:2px;',
                                                                            selectizeInput("BPFeatSelection", label = "Select Feature:", choices = NULL,
                                                                                           multiple = F, selected = 1,width = "100%")
                                                                     ),
                                                                     column(3, style = 'padding-left:2px;',
                                                                            selectInput("BPlogOpt","Log:", choices = c("No Log","Log2","Log2+1","Log10","Log10+1"),
                                                                                        selected = "Log2+1")
                                                                     )
                                                                   ),
                                                                   fluidRow(
                                                                     column(5, style = 'padding-right:2px;',
                                                                            selectInput("BPplotstatComp","Stat Test Method:",
                                                                                        choices = c("none","wilcox.test","t.test","kruskal.test","anova"))
                                                                     ),
                                                                     column(4, style = 'padding-left:4px;padding-right:4px;',
                                                                            checkboxInput("BPplotsampledots","Include Dot Annotation", value = F)
                                                                     ),
                                                                     column(3, style = 'padding-left:4px;',
                                                                            numericInput("BPplotDotSize","Dot Size:", value = 2, step = 0.25)
                                                                     )
                                                                   ),
                                                                   fluidRow(
                                                                     column(5, style = 'padding-right:2px;',
                                                                            selectInput("BPplotXaxOrder","X-Axis Group Order",
                                                                                        choices = c("Ascending","Descending","Not Specificed"))
                                                                     ),
                                                                     column(4, style = 'padding-left:4px;padding-right:4px;',
                                                                            checkboxInput("BPremoveSingles","Remove groups 1 or less samples", value = T)
                                                                     ),
                                                                     column(3, style = 'padding-left:4px;',
                                                                            checkboxInput("BPflipBP","Flip Axis", value = F),
                                                                            radioButtons("BPorViolin",NULL,choices = c("Box Plot","Violin Plot"), selected = "Box Plot")
                                                                     )
                                                                   )
                                                          ),
                                                          #h4("Gene Selection:"),
                                                          #div(DT::dataTableOutput("GeneListTable2"), style = "font-size:10px"
                                                          #     )
                                                          # ),
                                                          tabPanel("Figure Settings",
                                                                   p(),
                                                                   fluidRow(
                                                                     column(6,
                                                                            selectInput("BPTheme","Theme:",
                                                                                        choices = c("Minimal" = "theme_minimal","Grey" = "theme_grey","BW" = "theme_bw",
                                                                                                    "Linedraw" = "theme_linedraw","Light" = "theme_light","Dark" = "theme_dark",
                                                                                                    "Classic" = "theme_classic","Void" = "theme_void","Test" = "theme_test"))
                                                                     ),
                                                                     column(6,
                                                                            selectInput("BPxAxisOrient","X-Axis Label Angle",
                                                                                        choices = c(90,45,0))
                                                                     )
                                                                   ),
                                                                   fluidRow(
                                                                     column(6,
                                                                            textInput("BPplotHeight","Plot Height:",value = "500px")
                                                                     ),
                                                                     column(6,
                                                                            textInput("BPplotWidth","Plot Width:",value = "100%")
                                                                     )
                                                                   ),
                                                                   fluidRow(
                                                                     column(4,
                                                                            numericInput("BPplot1XAxisSize","X-Axis Font:",
                                                                                         value = 16, step = 1)
                                                                     ),
                                                                     column(4,
                                                                            numericInput("BPplot1YAxisSize","Y-Axis Font:",
                                                                                         value = 16, step = 1)
                                                                     ),
                                                                     column(4,
                                                                            textInput("BPplot1YAxisLim","Y-Axis Limit:",
                                                                                      placeholder = "min,max")
                                                                     )
                                                                   ),
                                                                   h4("Figure Download Parameters"),
                                                                   fluidRow(
                                                                     column(4,
                                                                            numericInput("BPHeight","Height",value = 8)
                                                                     ),
                                                                     column(4,
                                                                            numericInput("BPWidth","Width",value = 10)
                                                                     ),
                                                                     column(4,
                                                                            selectInput("BPUnits","Units",choices = c("in","cm","mm","px"))
                                                                     )
                                                                   )
                                                                   #fluidRow(
                                                                   #  column(4,
                                                                   #         numericInput("boxplot1TitleSize","Title Font Size:",
                                                                   #                      value = 20, step = 1)
                                                                   #  ),
                                                                   #  column(4,
                                                                   #         numericInput("boxplot1AxisSize","Axis Font Size:",
                                                                   #                      value = 16, step = 1)
                                                                   #  ),
                                                                   #  column(4,
                                                                   #         numericInput("boxplotDot2", "Dot Size:",
                                                                   #                      min = 0, max = 5,
                                                                   #                      value = 0.5, step = .25)
                                                                   #  )
                                                                   #)
                                                          )
                                                        )
                                                        
                                       ),
                                       #### Barplot
                                       conditionalPanel(condition = "input.dataset == '4'",
                                                        tabsetPanel(
                                                          tabPanel("Data Parameters",
                                                                   p(),
                                                                   h4("Transform Data"),
                                                                   fluidRow(
                                                                     column(6, style = 'padding-right:4px;',
                                                                            checkboxInput("log2barplot","Log2(expr+1) Expression",value = T)
                                                                     ),
                                                                     column(6, style = 'padding-left:4px;',
                                                                            selectInput("errorbarplot","Error Bar Type",
                                                                                        choices = c("Standard Deviation","Standard Error","None"))
                                                                     )
                                                                   ),
                                                                   hr(),
                                                                   h4("Gene Selection:"),
                                                                   div(DT::dataTableOutput("GeneListTableBarPlot"), style = "font-size:10px")
                                                          ),
                                                          tabPanel("Figure Paramters",
                                                                   p(),
                                                                   textInput("barplotColoCodes","Color Code(s):",value = "", placeholder = "HEX or R Color Code(s) (Space Delim)"),
                                                                   hr(),
                                                                   fluidRow(
                                                                     column(6,
                                                                            checkboxInput("barplotsampledots","Include Dot Annotation", value = T)
                                                                     ),
                                                                     column(6,
                                                                            numericInput("barplotDotSize","Dot Size:", value = 0.5, step = 0.5)
                                                                     )
                                                                   ),
                                                                   hr(),
                                                                   fluidRow(
                                                                     column(6,
                                                                            textInput("barPlotYlim","Y-Axis Limits", value = "", placeholder = "min,max"),
                                                                            selectInput("barxAxisOrient","X-Axis Label Orientation",
                                                                                        choices = c(45,90,0))
                                                                     ),
                                                                     column(6,
                                                                            numericInput("barplotYbreaks","Y-Axis Breaks",value = "", step = 0.5),
                                                                            radioButtons("barplotXaxOrder","X-Axis Group Order",
                                                                                         choices = c("Ascending","Descending"))
                                                                     )
                                                                   ),
                                                                   hr(),
                                                                   fluidRow(
                                                                     column(6,
                                                                            numericInput("barplot1TitleSize","Title Font Size:",
                                                                                         value = 20, step = 1)
                                                                     ),
                                                                     column(6,
                                                                            numericInput("barplot1AxisSize","Axis Font Size:",
                                                                                         value = 16, step = 1)
                                                                     )
                                                                   )
                                                          )
                                                        )
                                       )
                                     ),
                                     ### Mainpanel -----------------------------------------
                                     mainPanel(
                                       tabsetPanel(
                                         id = "dataset",
                                         #### Heatmaps
                                         ##### MVG Heatmaps
                                         tabPanel("Heatmaps",
                                                  tabsetPanel(
                                                    id = "datasetheat",
                                                    tabPanel("Unsupervised Clustering Heatmap",
                                                             jqui_resizable(plotOutput("heatmap1", width = "100%", height = "1000px")),
                                                             p(),
                                                             fluidRow(
                                                               downloadButton("dnldPlotSVG_heat1","Download as SVG"),
                                                               #downloadButton("dnldPlotPDF_heat1","Download as PDF")
                                                             ),
                                                             value = 1),
                                                    ##### Custom Heatmaps
                                                    tabPanel("Custom Heatmap",
                                                             withSpinner(jqui_resizable(plotOutput("heatmap2", width = "100%", height = "1000px")), type = 6),
                                                             p(),
                                                             fluidRow(
                                                               downloadButton("dnldPlotSVG_heat2","Download as SVG"),
                                                               #downloadButton("dnldPlotPDF_heat2","Download as PDF")
                                                             ),
                                                             value = 2),
                                                    ##### DEG Heatmaps
                                                    tabPanel("DEG Heatmap",
                                                             withSpinner(jqui_resizable(plotOutput("heatmap3", width = "100%", height = "1000px")), type = 6),
                                                             p(),
                                                             fluidRow(
                                                               downloadButton("dnldPlotSVG_heat3","Download as SVG"),
                                                               #downloadButton("dnldPlotPDF_heat3","Download as PDF")
                                                             ),
                                                             value = 3),
                                                    ##### Average Expr Heatmaps
                                                    tabPanel("Custom Genes Average Expression Heatmap",
                                                             jqui_resizable(plotOutput("avgheatmap1Cust", width = "100%", height = "1000px")),
                                                             #withSpinner(jqui_resizable(plotOutput("avgheatmap1Cust", width = "100%", height = "1000px")), type = 6),
                                                             p(),
                                                             fluidRow(
                                                               downloadButton("dnldPlotSVG_heat5","Download as SVG"),
                                                               downloadButton("dnldPlotPDF_heat5","Download as PDF")
                                                             ),
                                                             value = 5)
                                                  ),
                                                  value = 1
                                         ),
                                         #### Scatter Plot
                                         tabPanel("Gene Scatter Plot",
                                                  p(),
                                                  withSpinner(jqui_resizable(plotlyOutput("geneScatter0", height = "500px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldPlotSVG_scatter","Download as SVG"),
                                                    downloadButton("dnldPlotPDF_scatter","Download as PDF")
                                                  ),
                                                  p(),
                                                  DT::dataTableOutput("geneScatterTable"),
                                                  downloadButton("geneScatterDownload", "Download Non-log2 Transformed .tsv"),
                                                  value = 2),
                                         #### Box Plot
                                         tabPanel("Box Plot",
                                                  p(),
                                                  #withSpinner(jqui_resizable(plotOutput('boxplot3', width = "100%", height = "600px")), type = 6),
                                                  #jqui_resizable(plotOutput('boxplot3', width = "100%", height = "600px")),
                                                  jqui_resizable(plotOutput('uncorrected_Box_plot', width = "100%", height = "600px")),
                                                  fluidRow(
                                                    downloadButton("dnldPlotSVG_exprBox","Download as SVG"),
                                                    downloadButton("dnldPlotPDF_exprBox","Download as PDF")
                                                  ),
                                                  p(),
                                                  DT::dataTableOutput("DataExplor_Box_plot_df"),
                                                  downloadButton("dnldDataExplor_Box_plot_df","Download table"),
                                                  value = 3),
                                         #### Bar Plot
                                         tabPanel("Bar Plot",
                                                  p(),
                                                  jqui_resizable(plotOutput('barplot', width = "100%", height = "600px")),
                                                  fluidRow(
                                                    downloadButton("dnldPlotSVG_exprBar","Download as SVG"),
                                                    downloadButton("dnldPlotPDF_exprBar","Download as PDF")
                                                  ),
                                                  value = 4)
                                         
                                       )
                                     )
                                   )
                                 )
)

## DGE Tab -----------------------------------------
DGE_tab <- tabPanel("Differential Expression Analysis",
                    fluidPage(
                      title = "Differential Expression Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          width = 3,
                          
                          #uiOutput("rendDEGSubsetCol"),
                          #uiOutput("rendDEGSubsetCrit"),
                          #uiOutput("rendDEGGroupCol"),
                          
                          ### Sidebar -----------------------------------------
                          #### Volcano and MA plot
                          
                          conditionalPanel(condition = "input.datasettwo == '1' || input.datasettwo == '2'",
                                           tabsetPanel(
                                             tabPanel("Data Input",
                                                      h4("Condition Selection:"),
                                                      uiOutput("rendDEGGroupCol"),
                                                      uiOutput("rendDESeqDesignColRef_vol"),
                                                      radioButtons("volcanoCompChoice","Comparison Method:",
                                                                   choices = c("Limma: Two groups","Limma: One group"), inline = T),
                                                      conditionalPanel(condition = "input.volcanoCompChoice == 'Limma: Two groups'",
                                                                       fluidRow(
                                                                         column(6,
                                                                                uiOutput("rendcomparisonA2")
                                                                         ),
                                                                         column (6,
                                                                                 uiOutput("rendcomparisonB2"))
                                                                       )
                                                      ),
                                                      conditionalPanel(condition = "input.volcanoCompChoice == 'Limma: One group'",
                                                                       uiOutput("rendcomparisonA2_one")
                                                      ),
                                                      #conditionalPanel(condition = "input.RawCountQuantNorm.includes('Normalize Raw Counts') & input.volcanoCompChoice == 'DESeq2'",
                                                      #                 uiOutput("rendDESeqDesignCol_vol"),
                                                      #                 uiOutput("rendDESeqDesignColRef_vol")
                                                      #),
                                                      selectizeInput("DEGCovarSelect","Select Covariates:", choices = NULL, selected = 1, multiple = T),
                                                      h4("Threshold Parameters"),
                                                      fluidRow(
                                                        column(6,
                                                               numericInput("fc_cutoff", "LogFC Threshold",
                                                                            min = 0, max = 5, step = 0.1, value = 0)
                                                        ),
                                                        column(6,
                                                               numericInput("p_cutoff", "P.Value Threshold:",
                                                                            min = 0, max = 10, step = 0.1, value = 0.05)
                                                        )
                                                      ),
                                                      conditionalPanel(condition = "input.datasettwo == '1'" ,
                                                                       radioButtons("PorAdjPval","Y-Axis Scale:",choices = c("-log10(P.value)","-log10(Adjusted P.value)"),
                                                                                    inline = T)
                                                      ),
                                                      h4("Gene Selection Parameters:"),
                                                      numericInput("top_x", "Number of Top Hits:", value = 10,
                                                                   min = 0),
                                                      selectizeInput("userGeneSelec", "User Selected Hits:",
                                                                     choices = NULL, multiple = T, selected = 1),
                                                      textInput("userGeneSelec2", "Text Input of Gene List (space or tab delimited):", value = "")
                                             ),
                                             tabPanel("Figure Settings",
                                                      h4("Font Sizes"),
                                                      fluidRow(
                                                        column(4,
                                                               numericInput("VolMAAxisSize","Axis Lable",
                                                                            value = 18, step = 1)
                                                        ),
                                                        column(4,
                                                               numericInput("VolMATickSize","Axis Tick",
                                                                            value = 16, step = 1)
                                                        ),
                                                        column(4,
                                                               numericInput("VolMAAnnoSize","Annotation Text",
                                                                            value = 6, step = 1)
                                                        )
                                                      ),
                                                      h4("Download Parameters"),
                                                      fluidRow(
                                                        column(4,
                                                               numericInput("VolMAdnldHeight","Height",
                                                                            value = 8, step = 1, min = 1)
                                                        ),
                                                        column(4,
                                                               numericInput("VolMAdnldWidth","Width",
                                                                            value = 10, step = 1, min = 1)
                                                        ),
                                                        column(4,
                                                               selectInput("VolMAdnldSizeUnits","Units",
                                                                           choices = c("in","cm","mm","px"))
                                                        )
                                                      )
                                             )
                                           )
                                           
                          ),
                          
                          #### Boxplot
                          
                          conditionalPanel(condition = "input.datasettwo == '3'",
                                           tabsetPanel(
                                             tabPanel("Data Input",
                                                      p(),
                                                      uiOutput("rendBPgroupCriteria2"),
                                                      uiOutput("rendBPgroupSelection2"),
                                                      radioButtons("BPFeatureCategory2","Select Feature Type:",
                                                                   choices = c("Matrix Features","Meta Features"), inline = T),
                                                      uiOutput("rendImmuneDeconvMethods2"),
                                                      fluidRow(
                                                        column(9, style = 'padding-right:2px;',
                                                               selectizeInput("BPFeatSelection2", label = "Select Feature:", choices = NULL,
                                                                              multiple = F, selected = 1,width = "100%")
                                                        ),
                                                        column(3, style = 'padding-left:2px;',
                                                               selectInput("BPlogOpt2","Log:", choices = c("No Log","Log2","Log2+1","Log10","Log10+1"),
                                                                           selected = "Log2+1")
                                                        )
                                                      ),
                                                      fluidRow(
                                                        column(5, style = 'padding-right:2px;',
                                                               selectInput("BPplotstatComp2","Stat Test Method:",
                                                                           choices = c("none","wilcox.test","t.test","kruskal.test","anova"))
                                                        ),
                                                        column(4, style = 'padding-left:4px;padding-right:4px;',
                                                               checkboxInput("BPplotsampledots2","Include Dot Annotation", value = F)
                                                        ),
                                                        column(3, style = 'padding-left:4px;',
                                                               numericInput("BPplotDotSize2","Dot Size:", value = 2, step = 0.25)
                                                        )
                                                      ),
                                                      fluidRow(
                                                        column(5, style = 'padding-right:2px;',
                                                               selectInput("BPplotXaxOrder2","X-Axis Group Order",
                                                                           choices = c("Ascending","Descending","Not Specificed"))
                                                        ),
                                                        column(4, style = 'padding-left:4px;padding-right:4px;',
                                                               checkboxInput("BPremoveSingles2","Remove groups 1 or less samples", value = T)
                                                        ),
                                                        column(3, style = 'padding-left:4px;',
                                                               checkboxInput("BPflipBP2","Flip Axis", value = F),
                                                               radioButtons("BPorViolin2",NULL,choices = c("Box Plot","Violin Plot"), selected = "Box Plot")
                                                        )
                                                      )
                                                      #h4("Gene Selection:"),
                                                      #div(DT::dataTableOutput("GeneListTable"), style = "font-size:10px")
                                             ),
                                             tabPanel("Figure Settings",
                                                      p(),
                                                      fluidRow(
                                                        column(6,
                                                               selectInput("BPTheme2","Theme:",
                                                                           choices = c("Minimal" = "theme_minimal","Grey" = "theme_grey","BW" = "theme_bw",
                                                                                       "Linedraw" = "theme_linedraw","Light" = "theme_light","Dark" = "theme_dark",
                                                                                       "Classic" = "theme_classic","Void" = "theme_void","Test" = "theme_test"))
                                                        ),
                                                        column(6,
                                                               selectInput("BPxAxisOrient2","X-Axis Label Angle",
                                                                           choices = c(90,45,0))
                                                        )
                                                      ),
                                                      fluidRow(
                                                        column(6,
                                                               textInput("BPplotHeight2","Plot Height:",value = "500px")
                                                        ),
                                                        column(6,
                                                               textInput("BPplotWidth2","Plot Width:",value = "100%")
                                                        )
                                                      ),
                                                      fluidRow(
                                                        column(4,
                                                               numericInput("BPplot1XAxisSize2","X-Axis Font:",
                                                                            value = 16, step = 1)
                                                        ),
                                                        column(4,
                                                               numericInput("BPplot1YAxisSize2","Y-Axis Font:",
                                                                            value = 16, step = 1)
                                                        ),
                                                        column(4,
                                                               textInput("BPplot1YAxisLim2","Y-Axis Limit:",
                                                                         placeholder = "min,max")
                                                        )
                                                      ),
                                                      h4("Figure Download Parameters"),
                                                      fluidRow(
                                                        column(4,
                                                               numericInput("BPHeight2","Height",value = 8)
                                                        ),
                                                        column(4,
                                                               numericInput("BPWidth2","Width",value = 10)
                                                        ),
                                                        column(4,
                                                               selectInput("BPUnits2","Units",choices = c("in","cm","mm","px"))
                                                        )
                                                      )
                                                      #fluidRow(
                                                      #  column(4,
                                                      #         numericInput("boxplot2TitleSize","Title Font Size",
                                                      #                      value = 20, step = 1)
                                                      #  ),
                                                      #  column(4,
                                                      #         numericInput("boxplot2AxisSize","Axis Font Size",
                                                      #                      value = 16, step = 1)
                                                      #  ),
                                                      #  column(4,
                                                      #         numericInput("boxplotDot", "Dot Size:",
                                                      #                      min = 0, max = 5,
                                                      #                      value = 0.5, step = .25)
                                                      #  )
                                                      #)
                                             )
                                           )
                          ),
                          
                          #### Avg Expr Scatter
                          
                          conditionalPanel(condition = "input.datasettwo == '4'",
                                           tabsetPanel(
                                             tabPanel("Data Input",
                                                      p(),
                                                      uiOutput("rendAvgGroupCol"),
                                                      fluidRow(
                                                        column(6,
                                                               uiOutput("rendcomparisonA2.avg")
                                                        ),
                                                        column(6,
                                                               uiOutput("rendcomparisonB2.avg")
                                                        )
                                                      ),
                                                      checkboxInput("AvgExpLogFC","Log2 Transform Expression Data", value = TRUE),
                                                      h4("Gene Annotation Selection"),
                                                      selectizeInput("scatterGeneSelec", "Select Genes to Annotate in Plot:",
                                                                     choices = NULL, multiple = T, selected = 1),
                                                      textInput("gsSelection2", "Text Input of Genes (space or tab delimited):", value = ""),
                                                      uiOutput("hover_info3")
                                             ),
                                             tabPanel("Figure Settings",
                                                      p(),
                                                      fluidRow(
                                                        column(6,
                                                               numericInput("AvgExprScatterTitleSize","Title Font Size",
                                                                            value = 18, step = 1)
                                                        ),
                                                        column(6,
                                                               numericInput("AvgExprScatterAxisSize","Axis Lable Font Size",
                                                                            value = 16, step = 1)
                                                        )
                                                      )
                                             )
                                           )
                          ),
                          
                          
                          #### DEG Table
                          
                          conditionalPanel(condition = "input.datasettwo == '6'",
                                           h4("Condition Selection:"),
                                           uiOutput("rendDEGtabGroupCol"),
                                           radioButtons("volcanoCompChoice2","Comparison Method:",
                                                        choices = c("Limma: Two groups","Limma: One group"), inline = T),
                                           conditionalPanel(condition = "input.volcanoCompChoice2 == 'Limma: Two groups'",
                                                            fluidRow(
                                                              column(6,
                                                                     uiOutput("rendcomparisonA2.DEG")
                                                              ),
                                                              column (6,
                                                                      uiOutput("rendcomparisonB2.DEG"))
                                                            )
                                           ),
                                           conditionalPanel(condition = "input.volcanoCompChoice2 == 'Limma: One group'",
                                                            uiOutput("rendcomparisonA2.DEG_one")
                                           ),
                                           conditionalPanel(condition = "input.RawCountQuantNorm.includes('Normalize Raw Counts') & input.volcanoCompChoice2 == 'DESeq2'",
                                                            uiOutput("rendDESeqDesignCol_tab"),
                                                            uiOutput("rendDESeqDesignColRef_tab")
                                           ),
                                           selectizeInput("DEGCovarSelectTab","Select Covariates:", choices = NULL, selected = 1, multiple = T),
                                           h4("Download DEG Table as GMT File"),
                                           textInput("DEGfileName", "File Name for Download:",value = "DEGgeneSet"),
                                           textInput("DEGGeneSetName1", "Upregulated Gene Set Name:",value = NULL),
                                           conditionalPanel(condition = "input.UpDnChoice == 'UpAndDown_Regulated'",
                                                            textInput("DEGGeneSetName2", "Downregulated Gene Set Name:",value = NULL)
                                           ),
                                           fluidRow(
                                             column(6,
                                                    numericInput("fc_cutoff2", "LogFC Threshold",
                                                                 min = 0, max = 5, step = 0.1, value = 0)
                                             ),
                                             column(6,
                                                    numericInput("p_cutoff2", "P.Value Cutoff:",
                                                                 min = 0, max = 10, step = 0.1, value = 0.05)
                                             )
                                           ),
                                           selectInput("UpDnChoice","Up-regulated or Down-regulated:",
                                                       choices = c("UpAndDown_Regulated","Up_Regulated","Down_Regulated")),
                                           numericInput("top_x2", "Number of Top Hits:", value = 100),
                                           downloadButton("DEGgmtDownload", "DEG Geneset .gmt"),
                                           downloadButton("DEGtsvDownload", "DEG Geneset .tsv")
                          )
                        ),
                        ### Main panel -----------------------------------------
                        mainPanel(
                          tabsetPanel(
                            id = "datasettwo",
                            #### Volcano plot
                            tabPanel("Volcano Plot",
                                     p(),
                                     fluidRow(
                                       column(8,
                                              verbatimTextOutput("VolGroupsText"),
                                              jqui_resizable(plotOutput('Volcano3', width = "800px", height = "550px",
                                                                        hover = hoverOpts("plot_hover", delay = 10, delayType = "debounce")))
                                       ),
                                       column(4,
                                              uiOutput("hover_info")
                                       )
                                     ),
                                     fluidRow(
                                       downloadButton("dnldPlotSVG_vol","Download as SVG"),
                                       downloadButton("dnldPlotPDF_vol","Download as PDF")
                                     ),
                                     value = 1),
                            #### MA plot
                            tabPanel("MA Plot",
                                     p(),
                                     fluidRow(
                                       column(8,
                                              verbatimTextOutput("MAGroupsText"),
                                              jqui_resizable(plotOutput('MAPlot1', width = "800px", height = "550px",
                                                                        hover = hoverOpts("plot_hover2", delay = 10, delayType = "debounce")))
                                       ),
                                       column(4,
                                              uiOutput("hover_info2")
                                       )
                                     ),
                                     #verbatimTextOutput("MAGroupsText"),
                                     fluidRow(
                                       downloadButton("dnldPlotSVG_MA","Download as SVG"),
                                       downloadButton("dnldPlotPDF_MA","Download as PDF")
                                     ),
                                     #uiOutput("hover_info2"),
                                     value = 2),
                            #### Box plot
                            tabPanel("Box Plots",
                                     p(),
                                     #withSpinner(jqui_resizable(plotOutput('boxplot1', width = "100%", height = "600px")), type = 6),
                                     #jqui_resizable(plotOutput('boxplot1', width = "100%", height = "600px")),
                                     jqui_resizable(plotOutput('uncorrected_Box_plot2', width = "100%", height = "600px")),
                                     fluidRow(
                                       downloadButton("dnldPlotSVG_exprBox2","Download as SVG"),
                                       downloadButton("dnldPlotPDF_exprBox2","Download as PDF")
                                     ),
                                     p(),
                                     DT::dataTableOutput("DataExplor_Box_plot_df2"),
                                     downloadButton("dnldDataExplor_Box_plot_df2","Download table"),
                                     value = 3),
                            #### Avg Expr Scatter plot
                            tabPanel("Average Expression Scatter Plot",
                                     p(),
                                     jqui_resizable(plotOutput("AvggeneScatter2", height = "500px",
                                                               hover = hoverOpts("plot_hover3", delay = 10, delayType = "debounce"))),
                                     #withSpinner(jqui_resizable(plotOutput("AvggeneScatter2", height = "500px")), type = 6),
                                     #withSpinner(jqui_resizable(plotOutput("AvggeneScatter2", height = "500px",
                                     #                                      hover = hoverOpts("plot_hover", delay = 10, delayType = "debounce"))), type = 6),
                                     fluidRow(
                                       downloadButton("dnldPlotSVG_AvgScatter","Download as SVG"),
                                       downloadButton("dnldPlotPDF_AvgScatter","Download as PDF")
                                     ),
                                     p(),
                                     DT::dataTableOutput("AvggeneScatterTable"),
                                     downloadButton("AvggeneScatterDownload", "Download Non-log2 Transformed .tsv"),
                                     value = 4),
                            #### DEG Table
                            tabPanel("DEG Table",
                                     p(),
                                     verbatimTextOutput("degtext"),
                                     p(),
                                     downloadButton("DEGtableDownload", "Download DEG table .txt"),
                                     p(),
                                     div(DT::dataTableOutput("DEGtable1"), style = "font-size:12px"),
                                     value = 6)
                          )
                        )
                      )
                    )
)

## GSEA Tab -----------------------------------------
GSEA_tab <- tabPanel("GSEA Analysis",
                     fluidPage(
                       title = "GSEA Analysis",
                       sidebarLayout(
                         sidebarPanel(
                           width = 3,
                           
                           #uiOutput("rendGSEASubsetCol"),
                           #uiOutput("rendGSEASubsetCrit"),
                           
                           ### Sidebar -----------------------------------------
                           tabsetPanel(id = "GSEA",
                                       tabPanel("Data Input",
                                                p(),
                                                conditionalPanel(condition = "input.datasetthree != '6'",
                                                                 h4("Condition Selection:"),
                                                                 uiOutput("rendGSEAGroupCol")
                                                ),
                                                conditionalPanel(condition = "input.datasetthree == '2'",
                                                                 uiOutput("rendGSEAGroupColMulti")
                                                ),
                                                conditionalPanel(condition = "input.datasetthree == '6'",
                                                                 conditionalPanel(condition = "input.datasetssheat == '2' | input.datasetssheat == '1'",
                                                                                  uiOutput("rendssGSEAGroupColMulti")
                                                                 ),
                                                                 conditionalPanel(condition = "input.datasetssheat == '3'",
                                                                                  uiOutput("rendssGSEAGroupColAvg")
                                                                 )
                                                ),
                                                conditionalPanel(condition = "input.datasetthree == '1' | input.datasetthree == '2' | input.datasetthree == '4'",
                                                                 radioButtons("volcanoCompChoice3","Comparison Method:",
                                                                              choices = c("Two groups","One group"), inline = T),
                                                                 conditionalPanel(condition = "input.volcanoCompChoice3 == 'Two groups'",
                                                                                  fluidRow(
                                                                                    column(6,
                                                                                           uiOutput("rendcomparisonA")
                                                                                    ),
                                                                                    column (6,
                                                                                            uiOutput("rendcomparisonB")
                                                                                    )
                                                                                  )
                                                                 ),
                                                                 conditionalPanel(condition = "input.volcanoCompChoice3 == 'One group'",
                                                                                  uiOutput("rendcomparisonA_one")
                                                                 ),
                                                                 #fluidRow(
                                                                 #  column(6,
                                                                 #         uiOutput("rendcomparisonA")
                                                                 #  ),
                                                                 #  column(6,
                                                                 #         uiOutput("rendcomparisonB")
                                                                 #  )
                                                                 #)
                                                ),
                                                conditionalPanel(condition = "input.datasetthree == '4'",
                                                                 h4("GSEA Threshold Parameter:"),
                                                                 numericInput("userPval", "Pvalue Cutoff", value = 1.0, width = '100%')
                                                ),
                                                conditionalPanel(condition = "input.datasetthree == '5' | input.datasetthree =='6'",
                                                                 p(),
                                                                 selectInput("ssGSEAtype","Choose ssGSEA Method",
                                                                             choices = c("ssgsea","gsva","zscore","plage"))
                                                ),
                                                conditionalPanel(condition = "input.datasetthree == '5'",
                                                                 selectInput("boxplotcompare", "Boxplot Stat Compare Method:",
                                                                             choices = c("none","wilcox.test","t.test","kruskal.test","anova"))
                                                ),
                                                conditionalPanel(condition = "input.datasetthree == '1' | input.datasetthree == '2' | input.datasetthree == '4' | input.datasetthree == '5'",
                                                                 tabsetPanel(
                                                                   id = "tables",
                                                                   tabPanel("Gene Sets",
                                                                            p(),
                                                                            selectInput("GeneSetCat_Select","Select Geneset Category:", choices = NULL, selected = 1),
                                                                            radioButtons("GeneSetMultiOpt",NULL, choices = c("Select single gene set","Select multiple gene sets"), inline = T),
                                                                            div(DT::dataTableOutput("GeneSetTable"), style = "font-size:10px"),
                                                                            value = 1
                                                                   ),
                                                                   #tabPanel("MSigDB Gene Sets",
                                                                   #         p(),
                                                                   #         div(DT::dataTableOutput("msigdbTable"), style = "font-size:10px;"),
                                                                   #         value = 1),
                                                                   #tabPanel(paste("Cell Marker"),
                                                                   #         p(),
                                                                   #         div(DT::dataTableOutput("tab2table"), style = "font-size:10px;"),
                                                                   #         value = 3),
                                                                   tabPanel("Use your own gene set",
                                                                            p(),
                                                                            fileInput("userGeneSet","Gene Set Upload", accept = c(".gmt",".tsv",".txt",".RData")),
                                                                            div(DT::dataTableOutput("userGeneSetTable"), style = "font-size:10px"),
                                                                            #fluidRow(
                                                                            #  column(9,
                                                                            #         fileInput("userGeneSet","Gene Set Upload", accept = c(".gmt",".tsv",".txt",".RData"))
                                                                            #         #uiOutput("user.gmt")
                                                                            #  ),
                                                                            #  column(3,
                                                                            #         checkboxInput("UserGSheaderCheck","Header",value = T)
                                                                            #  )
                                                                            #),
                                                                            #uiOutput("user.RDataButton"),
                                                                            #uiOutput("RDataMessage"),
                                                                            #uiOutput("user.GStable"),
                                                                            value = 5)
                                                                 )
                                                ),
                                                conditionalPanel(condition = "input.datasetthree =='6'",
                                                                 p(),
                                                                 h4("Gene Set Selection"),
                                                                 #uiOutput("rendssgseaHeatGS"),
                                                                 selectInput("GeneSetCat_Select_heat","Select Geneset Category:", choices = NULL, selected = 1, multiple = T),
                                                                 selectizeInput("ssgseaHeatGS","Select Gene Sets:", choices = NULL, selected = 1, multiple = T),
                                                                 h4("Sample Selection"),
                                                                 selectizeInput("userheatsampSS", "",choices = NULL, multiple = T, selected = 1),
                                                )
                                       ),
                                       tabPanel("Figure Settings",
                                                p(),
                                                conditionalPanel(condition = "input.datasetthree == '2'",
                                                                 fluidRow(
                                                                   column(6,
                                                                          numericInput("GSEAheatmapHeight1","Download Height (in)",value = 20)
                                                                   ),
                                                                   column(6,
                                                                          numericInput("GSEAheatmapWidth1","Download Width (in)",value = 15)
                                                                   )
                                                                 )
                                                ),
                                                conditionalPanel(condition = "input.datasetthree == '6'",
                                                                 conditionalPanel("input.datasetssheat == '3'",
                                                                                  selectInput("ColorPalette_gseaHeat", "Select Color Palette:",
                                                                                              choices = c("Red/Blue" = "original",
                                                                                                          "OmniBlueRed" = "OmniBlueRed",
                                                                                                          "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                                                                          "Green/Black/Red" = "GreenBlackRed",
                                                                                                          "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                                                                          "Viridis" = "Viridis","Plasma" = "Plasma",
                                                                                                          "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens"))
                                                                 ),
                                                                 fluidRow(
                                                                   column(6,
                                                                          checkboxInput("ShowColNamesSSheat","Show Heatmap Column Names", value = T),
                                                                          numericInput("heatmapFont1.c", "Heatmap Column Font Size:",
                                                                                       min = 5, max = 75,
                                                                                       value = 12, step = 1),
                                                                          checkboxInput("clustcolsSSheat","Cluster Heatmap Columns", value = F)
                                                                   ),
                                                                   column(6,
                                                                          checkboxInput("ShowRowNames1SSheat","Show Heatmap Row Names", value = T),
                                                                          numericInput("heatmapFont1.r", "Heatmap Row Font Size:",
                                                                                       min = 5, max = 75,
                                                                                       value = 10, step = 1),
                                                                          checkboxInput("clustrowsSSheat","Cluster Heatmap Rows", value = F)
                                                                   )
                                                                 ),
                                                                 uiOutput("rendClusterMethodSSheat"),
                                                                 fluidRow(
                                                                   column(6,
                                                                          numericInput("ssgseaheatmapHeight","Download Height (in)",value = 20)
                                                                   ),
                                                                   column(6,
                                                                          numericInput("ssgseaheatmapWidth","Download Width (in)",value = 15)
                                                                   )
                                                                 )
                                                ),
                                                conditionalPanel(condition = "input.datasetthree == '5'",
                                                                 fluidRow(
                                                                   column(4,
                                                                          numericInput("gseaBoxTitleSize","Title Font Size",
                                                                                       value = 20, step = 1)
                                                                   ),
                                                                   column(4,
                                                                          numericInput("gseaBoxAxisSize","Axis Lable Font Size",
                                                                                       value = 16, step = 1)
                                                                   ),
                                                                   column(4,
                                                                          numericInput("boxplotDotss", "Dot Size:",
                                                                                       min = 0, max = 5,
                                                                                       value = 0.75, step = .25)
                                                                   )
                                                                 )
                                                )
                                                
                                       )
                           )
                         ),
                         ### Main Panel -----------------------------------------
                         mainPanel(
                           span(textOutput("GroupSelectionError"), style="color:red"),
                           tabsetPanel(
                             id = "datasetthree",
                             #### Enrichment Plot
                             tabPanel("Enrichment Plot",
                                      h3("GSEA Enrichment Plot"),
                                      verbatimTextOutput("NESandPval"),
                                      fluidRow(
                                        column(6,
                                               h3(""),
                                               withSpinner(plotOutput("enrichplot0", width = "500px", height = "450px"), type = 6),
                                               fluidRow(
                                                 downloadButton("dnldPlotSVG_gsea","Download as SVG"),
                                                 downloadButton("dnldPlotPDF_gsea","Download as PDF")
                                               )
                                        ),
                                        column(6,
                                               h3("Leading Edge Genes (~Signal2Noise Ranking)"),
                                               div(DT::dataTableOutput("LeadingEdgeGenes"), style = "font-size:12px; height:400px; overflow-y: scroll"),
                                               p(),
                                               downloadButton("LEGdownload", "Download .tsv")
                                        )
                                      ),
                                      value = 1),
                             #### GSEA Heatmap
                             tabPanel("GSEA Heatmap",
                                      withSpinner(jqui_resizable(plotOutput("heatmap0", width = "100%", height = "1000px")), type = 6),
                                      fluidRow(
                                        downloadButton("dnldPlotSVG_gseaHeat","Download as SVG")
                                      ),
                                      value = 2),
                             #### Enrich Sig Table
                             tabPanel("Generate Enriched Signatures Table",
                                      p(),
                                      p("Please note this may take several minutes depending on size and quantity of gene sets in GMT file."),
                                      #uiOutput("genESTbutton"),
                                      actionButton("GenerateEST", "Generate enriched signature table for selected gene set"),
                                      p(),
                                      withSpinner(DT::dataTableOutput("enrich_sig_table_gen"), type = 6),
                                      downloadButton("enrich_sig_download.u","Download .tsv"),
                                      value = 4),
                             #### ssGSEA Box Plot
                             tabPanel("ssGSEA Boxplots",
                                      p(),
                                      jqui_resizable(plotOutput('boxplot2', width = "100%", height = "500px")),
                                      fluidRow(
                                        downloadButton("dnldPlotSVG_gseaBox","Download as SVG"),
                                        downloadButton("dnldPlotPDF_gseaBox","Download as PDF")
                                      ),
                                      p(),
                                      DT::dataTableOutput("ssGSEAtable"),
                                      downloadButton("ssGSEAdownload", "Download .tsv"),
                                      value = 5),
                             #### ssGSEA Heatmaps
                             tabPanel("ssGSEA Heatmap",
                                      tabsetPanel(
                                        id = "datasetssheat",
                                        tabPanel("ssGSEA zScore Heatmap",
                                                 p(),
                                                 withSpinner(jqui_resizable(plotOutput('ssgseaheatmap', width = "100%", height = "800px")), type = 6),
                                                 fluidRow(
                                                   downloadButton("dnldPlotSVG_ssgseaHeat","Download as SVG")
                                                 ),
                                                 p(),
                                                 withSpinner(DT::dataTableOutput("ssgseaheatmap_df"), type = 6),
                                                 downloadButton("dnldssgseaheatmap_df","Download Table"),
                                                 value = 1
                                        ),
                                        tabPanel("ssGSEA Raw Difference Heatmap",
                                                 p(),
                                                 withSpinner(jqui_resizable(plotOutput('ssgseaheatmap2', width = "100%", height = "800px")), type = 6),
                                                 fluidRow(
                                                   downloadButton("dnldPlotSVG_ssgseaHeat2","Download as SVG")
                                                 ),
                                                 p(),
                                                 withSpinner(DT::dataTableOutput("ssgseaheatmap2_df"), type = 6),
                                                 downloadButton("dnldssgseaheatmap2_df","Download Table"),
                                                 value = 2
                                        ),
                                        tabPanel("Average Raw Difference ssGSEA Heatmap",
                                                 p(),
                                                 withSpinner(jqui_resizable(plotOutput('ssgseaheatmap3', width = "100%", height = "800px")), type = 6),
                                                 fluidRow(
                                                   downloadButton("dnldPlotSVG_ssgseaHeat3","Download as SVG"),
                                                   downloadButton("dnldPlotPDF_ssgseaHeat3","Download as PDF")
                                                 ),
                                                 p(),
                                                 withSpinner(DT::dataTableOutput("ssgseaheatmap3_df"), type = 6),
                                                 downloadButton("dnldssgseaheatmap3_df","Download Table"),
                                                 value = 3
                                        )
                                      ),
                                      value = 6)
                           )
                         )
                       )
                     )
)


# Render UI ---------------------------------------------------------------------------


if (Password_Protected) {
  ui <- navbarPage(paste("{ EASY Expression Analysis }", sep=" "),
                   id = "tabs",
                   collapsible = TRUE,
                   login_tab)
} else {
  ui <- navbarPage(paste("{ EASY Expression Analysis }", sep=" "),
                   id = "tabs",
                   collapsible = TRUE,
                   DataInput_tab,
                   Data_Exploration_tab,
                   DGE_tab,
                   GSEA_tab,
                   About_tab)
}


# Server -----------------------------------------------------------------------

server <- function(input, output, session) {
  
  
  
  
  # Password Protection --------------------------------------------------------
  # hack to add the logout button to the navbar on app launch
  if (Password_Protected) {
    insertUI(
      selector = ".navbar .container-fluid .navbar-collapse",
      ui = tags$ul(
        class="nav navbar-nav navbar-right",
        tags$li(
          div(
            style = "padding: 10px; padding-top: 8px; padding-bottom: 0;",
            #shinyauthr::logoutUI("logout")
            logoutUI("logout")
          )
        )
      )
    )
  }
  
  # call the shinyauthr login and logout server modules
  if (Password_Protected) {
    #credentials <- shinyauthr::loginServer(
    credentials <- loginServer(
      id = "login",
      data = user_base,
      user_col = "user",
      pwd_col = "password",
      sodium_hashed = FALSE,
      #sodium_hashed = TRUE,
      reload_on_logout = TRUE,
      log_out = reactive(logout_init())
    )
  } else {
    credentials <- reactive({
      list(user_auth = TRUE)
    })
  }
  
  #logout_init <- shinyauthr::logoutServer(
  if (Password_Protected) {
    logout_init <- logoutServer(
      id = "logout",
      active = reactive(credentials()$user_auth)
    )
  }
  
  observeEvent(credentials()$user_auth, {
    # if user logs in successfully
    if (Password_Protected) {
      if (credentials()$user_auth) {
        # remove the login tab
        removeTab("tabs", "login")
        # add home tab
        appendTab("tabs", DataInput_tab, select = TRUE)
        appendTab("tabs", Data_Exploration_tab, select = FALSE)
        appendTab("tabs", DGE_tab, select = FALSE)
        appendTab("tabs", GSEA_tab, select = FALSE)
        appendTab("tabs", About_tab, select = FALSE)
      }
    }
    
    if (credentials()$user_auth) {
      
      # Reactive Val Start ----------------------------------------------------------
      ProjectName_react <- reactiveVal(ProjectName)
      ExpressionMatrix_file_react <- reactiveVal(expression_file)
      MetaData_file_react <- reactiveVal(meta_file)
      MM_react <- reactiveVal(mouse)
      Subsetting_Feature_react <- reactiveVal(Subsetting_Feature)
      SubCrit_reactVal <- reactiveVal(1)
      Feature_Selected_react <- reactiveVal(Feature_Selected)
      SubCol_reactVal <- reactiveVal(Subsetting_Feature_react())
      SubCrit_reactVal <- reactiveVal(SubCrit_reactVal())
      metacol_reactVal <- reactiveVal(Feature_Selected_react())
      
      gse_id_react <- reactiveVal()
      recount_id_react <- reactiveVal()
      recount_obj <- reactiveVal()
      recount3_proj <- reactiveVal()
      user_upload <- reactiveVal()
      expr_input <- reactiveVal()
      meta_input <- reactiveVal()
      expr_raw <- reactiveVal()
      A_raw <- reactiveVal()
      ImmDeconv_react <- reactiveVal()
      geneList_raw <- reactiveVal()
      Gene_raw <- reactiveVal()
      FileCheckAlerts_react <- reactiveVal()
      
      # Data Input Tab ---------------------------------------------------------
      
      observe({
        if (isTruthy(Subsetting_Feature_react())) {
          SubCol_reactVal(Subsetting_Feature_react())
        }
      })
      observe({
        if (isTruthy(Feature_Selected_react())) {
          metacol_reactVal(Feature_Selected_react())
        }
      })
      
      output$rendExprFileInput <- renderUI({
        refresh <- input$UseExpData
        if (MM_react()) {
          HM_select <- "Mouse"
        } else {
          HM_select <- "Human"
        }
        fluidRow(
          column(8,
                 fileInput("ExprFileInput","Expression Matrix")
          ),
          column(4, style = "margin-top:15px",
                 shiny::radioButtons("HumanOrMouse",NULL,c("Human","Mouse","Mouse to Human Conv"), selected = HM_select)
          )
        )
      })
      output$rendRawCountNorm <- renderUI({
        if ("Normalize Raw Counts" %in% input$RawCountQuantNorm) {
          if (isTruthy(recount_id_react()))  {
            selectInput("RawCountNorm","Normalize Raw Counts:",
                        c("No Normalization" = "none",
                          "DESeq2" = "DESeq2",
                          "TPM" = "TPM",
                          "RPKM" = "RPKM",
                          "TMM" = "TMM",
                          "Upper Quartile" = "upperquartile"),
                        selected = "TPM")
          } else {
            selectInput("RawCountNorm","Normalize Raw Counts:",
                        c("No Normalization" = "none",
                          "DESeq2" = "DESeq2",
                          "TMM" = "TMM",
                          "Upper Quartile" = "upperquartile"),
                        selected = "TMM")
          }
          
        }
      })
      output$rendDESeqDesignCol <- renderUI({
        req(input$RawCountNorm)
        if ("Normalize Raw Counts" %in% input$RawCountQuantNorm) {
          if (input$RawCountNorm == "DESeq2") {
            meta <- meta_react()
            selectInput("DESeqDesignCol","Design Column:",choices = c("Select column",colnames(meta)[-1]))
          }
        }
        
      })
      output$rendDESeqDesignColRef <- renderUI({
        
        req(input$DESeqDesignCol)
        if ("Normalize Raw Counts" %in% input$RawCountQuantNorm) {
          if (input$RawCountNorm == "DESeq2") {
            if (input$DESeqDesignCol != "Select column") {
              meta <- meta_react()
              selectInput("DESeqDesignColRef","Design Reference:",choices = unique(meta[,input$DESeqDesignCol]))
            }
          }
        }
        
        
      })
      observe({
        if ("Normalize Raw Counts" %in% input$RawCountQuantNorm) {
          updateRadioButtons(session,"volcanoCompChoice", choices = c("Limma: Two groups","Limma: One group","DESeq2"), inline = T)
          updateRadioButtons(session,"volcanoCompChoice2", choices = c("Limma: Two groups","Limma: One group","DESeq2"), inline = T)
          updateRadioButtons(session,"volcanoCompChoice4", choices = c("Limma: Two groups","Limma: One group","DESeq2"), inline = T)
        }
        if (!"Normalize Raw Counts" %in% input$RawCountQuantNorm) {
          updateRadioButtons(session,"volcanoCompChoice", choices = c("Limma: Two groups","Limma: One group"), inline = T)
          updateRadioButtons(session,"volcanoCompChoice2", choices = c("Limma: Two groups","Limma: One group"), inline = T)
          updateRadioButtons(session,"volcanoCompChoice4", choices = c("Limma: Two groups","Limma: One group"), inline = T)
        }
      })
      
      
      #output$rendDESeqDesignCol_heat <- renderUI({
      #  req(input$RawCountNorm)
      #  if ("Normalize Raw Counts" %in% input$RawCountQuantNorm) {
      #    meta <- meta_react()
      #    selectInput("DESeqDesignCol_heat","Design Column:",choices = c("Select column",colnames(meta)[-1]))
      #  }
      #})
      output$rendDESeqDesignColRef_heat <- renderUI({
        req(metacol_reactVal())
        if ("Normalize Raw Counts" %in% input$RawCountQuantNorm) {
          if (input$volcanoCompChoice4 == "DESeq2") {
            meta <- meta_react()
            selectInput("DESeqDesignColRef_heat","Design Reference:",choices = unique(meta[,metacol_reactVal()]))
          }
        }
      })
      #output$rendDESeqDesignCol_vol <- renderUI({
      #  req(input$RawCountNorm)
      #  if ("Normalize Raw Counts" %in% input$RawCountQuantNorm) {
      #    meta <- meta_react()
      #    selectInput("DESeqDesignCol_vol","Design Column:",choices = c("Select column",colnames(meta)[-1]))
      #  }
      #})
      output$rendDESeqDesignColRef_vol <- renderUI({
        req(metacol_reactVal())
        if ("Normalize Raw Counts" %in% input$RawCountQuantNorm) {
          if (input$volcanoCompChoice == "DESeq2") {
            meta <- meta_react()
            selectInput("DESeqDesignColRef_vol","Design Reference:",choices = unique(meta[,metacol_reactVal()]))
          }
        }
      })
      #output$rendDESeqDesignCol_tab <- renderUI({
      #  req(input$RawCountNorm)
      #  if ("Normalize Raw Counts" %in% input$RawCountQuantNorm) {
      #    meta <- meta_react()
      #    selectInput("DESeqDesignCol_tab","Design Column:",choices = c("Select column",colnames(meta)[-1]))
      #  }
      #})
      output$rendDESeqDesignColRef_tab <- renderUI({
        req(metacol_reactVal())
        if ("Normalize Raw Counts" %in% input$RawCountQuantNorm) {
          meta <- meta_react()
          selectInput("DESeqDesignColRef_tab","Design Reference:",choices = unique(meta[,metacol_reactVal()]))
        }
      })
      output$rendClinFileInput <- renderUI({
        refresh <- input$UseExpData
        fileInput("ClinFileInput","Meta Data")
      })
      output$rendExprFilePrevHeader <- renderUI({
        req(ExpressionMatrix_file_react())
        h3("Expression File Preview")
      })
      output$rendExprHead <- renderUI({
        req(ExpressionMatrix_file_react())
        radioButtons("ExprHead",NULL, choices = c("View table head","View entire table"), inline = T)
      })
      output$rendClinFilePrevHeader <- renderUI({
        req(MetaData_file_react())
        h3("Meta File Preview")
      })
      output$rendClinHead <- renderUI({
        req(MetaData_file_react())
        radioButtons("ClinHead",NULL, choices = c("View table head","View entire table"), inline = T)
      })
      
      output$rendMetaColOfInterest <- renderUI({
        req(meta_input())
        selectizeInput("MetaColOfInterest","Meta column of interest:",choices = NULL, selected = "")
      })
      
      metaCol_ofInt_Selected <- reactive(input$MetaColOfInterest)
      observeEvent(meta_input(),{
        if (isTruthy(metaCol_ofInt_Selected())) {
          selectedFeat <- metaCol_ofInt_Selected()
        } else { selectedFeat <- "" }
        updateSelectizeInput(session,"MetaColOfInterest",choices = colnames(meta_input()), server = T,
                             selected = selectedFeat)
      })
      observe({
        if (isTruthy(input$MetaColOfInterest)) {
          meta <- meta_input()
          meta <- meta %>% relocate(any_of(input$MetaColOfInterest), .after = 1) %>% as.data.frame()
          meta_input(meta)
          Feature_Selected_react(input$MetaColOfInterest)
          metacol_reactVal(input$MetaColOfInterest)
        }
      })
      
      observe({
        req(meta_input())
        updateSelectizeInput(session,"ColToSplit",choices = colnames(meta_input()), server = T, selected = "")
        updateTextInput(session,"SplitDelim",value = "")
        updateTextInput(session,"SplitNewColNames",value = "")
      })
      observe({
        req(meta_input())
        updateSelectizeInput(session,"ColsToMerge",choices = colnames(meta_input()), server = T)
        updateTextInput(session,"MergeDelim",value = "")
        updateTextInput(session,"MergeNewColName",value = "")
      })
      
      # Observe URL input ----------------------------------------------------------
      
      # Observe URL Inputs
      observe({
        query <- parseQueryString(session$clientData$url_search)
        print(query)
        if (!is.null(query[['expr']]) && !is.null(query[['meta']])) {
          ProjectName_react(query[["proj"]])
          ExpressionMatrix_file_react(query[['expr']])
          MetaData_file_react(query[['meta']])
          if ('mouse' %in% names(query)) {
            if (query[['mouse']] == "") {
              MM_react(TRUE)
            } else {
              MM_react(query[['mouse']])
            }
          } else {
            MM_react(FALSE)
          }
          Subsetting_Feature_react(query[['subset']])
          Feature_Selected_react(query[['feature']])
        }
        if (!is.null(query[['gseid']]) && is.null(query[['expr']]) && is.null(query[['meta']]) && is.null(query[['recount']])) {
          gse_id_react(query[['gseid']])
          print(gse_id_react)
          ExpressionMatrix_file_react(query[['gseid']])
          MetaData_file_react(query[['gseid']])
          Subsetting_Feature_react(query[['subset']])
          Feature_Selected_react(query[['feature']])
        }
        if (is.null(query[['gseid']]) && is.null(query[['expr']]) && is.null(query[['meta']]) && !is.null(query[['recount']])) {
          recount_id_react(query[['recount']])
          print(recount_id_react)
          ExpressionMatrix_file_react(query[['recount']])
          MetaData_file_react(query[['recount']])
          Subsetting_Feature_react(query[['subset']])
          Feature_Selected_react(query[['feature']])
        }
      })
      
      ### Example Data ---------------------------------------------------------
      observeEvent(input$UseExpData, {
        ProjectName_react("TCGA_CHOL")
        ExpressionMatrix_file_react(ExampleExpr_File)
        MetaData_file_react(ExampleClin_File)
        metacol_reactVal("ajcc_pathologic_tumor_stage")
        Feature_Selected_react("ajcc_pathologic_tumor_stage")
        MM_react(FALSE)
        Subsetting_Feature_react("Select All Samples")
        SubCol_reactVal("Select All Samples")
        updateTextInput(session,"UserProjectName",value = "TCGA_CHOL")
      })
      
      observe({
        meta <- meta_react()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
          metacol_reactVal(metacol)
        }
        
      })
      
      ### User Upload ----------------------------------------------------------
      observe({
        ProjectName_react(input$UserProjectName)
        if (isTruthy(input$ExprFileInput$datapath)) {
          ExpressionMatrix_file_react(input$ExprFileInput$datapath)
        }
        if (isTruthy(input$ClinFileInput$datapath)) {
          MetaData_file_react(input$ClinFileInput$datapath)
        }
        if (isTruthy(input$HumanOrMouse)) {
          if (input$HumanOrMouse == "Mouse") {
            MM_react(TRUE)
          } else {
            MM_react(FALSE)
          }
        }
      })
      
      # Load in URL or given file
      observe({
        
        if (isTruthy(ExpressionMatrix_file_react()) & isTruthy(MetaData_file_react()) && !isTruthy(gse_id_react()) && !isTruthy(recount_id_react())) {
          # Expression file
          withProgress(message = "Processing", value = 0, {
            incProgress(0.5, detail = "Loading Expression Data")
            if (file.exists(ExpressionMatrix_file_react())) {
              if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("txt","tsv","zip","gz","TXT","TSV","ZIP","GZ")) {
                expr <- as.data.frame(read_delim(ExpressionMatrix_file_react(), delim = '\t', col_names = T))
              } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("csv","CSV")) {
                expr <- as.data.frame(read_delim(ExpressionMatrix_file_react(), delim = ',', col_names = T))
              } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("RData","rdata")) {
                expr <- loadRData(ExpressionMatrix_file_react())
              } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("rds","RDS")) {
                expr <- readRDS(ExpressionMatrix_file_react())
              }
            } else {
              if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("txt","tsv","TXT","TSV")) {
                expr <- as.data.frame(read_delim(url(ExpressionMatrix_file_react()), delim = '\t', col_names = T))
              } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("zip","gz","ZIP","GZ")) {
                expr <- as.data.frame(read_delim(getZip(ExpressionMatrix_file_react()), delim = '\t', col_names = T))
              } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("csv","CSV")) {
                expr <- as.data.frame(read_delim(url(ExpressionMatrix_file_react()), delim = ',', col_names = T))
              } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("RData","rdata")) {
                expr <- loadRData(url(ExpressionMatrix_file_react()))
              } else if (tools::file_ext(ExpressionMatrix_file_react()) %in% c("rds","RDS")) {
                expr <- readRDS(url(ExpressionMatrix_file_react()))
              }
            }
            incProgress(0.5, detail = "Complete!")
          })
          withProgress(message = "Processing", value = 0, {
            incProgress(0.5, detail = "Loading Meta Data")
            if (file.exists(MetaData_file_react())) {
              if (tools::file_ext(MetaData_file_react()) %in% c("txt","tsv","zip","gz","TXT","TSV","ZIP","GZ")) {
                meta <- as.data.frame(read_delim(MetaData_file_react(), delim = '\t', col_names = T))
              } else if (tools::file_ext(MetaData_file_react()) %in% c("csv","CSV")) {
                meta <- as.data.frame(read_delim(MetaData_file_react(), delim = ',', col_names = T))
              } else if (tools::file_ext(MetaData_file_react()) %in% c("RData","rdata")) {
                meta <- loadRData(MetaData_file_react())
              } else if (tools::file_ext(MetaData_file_react()) %in% c("rds","RDS")) {
                meta <- readRDS(MetaData_file_react())
              }
            } else {
              if (tools::file_ext(MetaData_file_react()) %in% c("txt","tsv","TXT","TSV")) {
                meta <- as.data.frame(read_delim(url(MetaData_file_react()), delim = '\t', col_names = T))
              } else if (tools::file_ext(MetaData_file_react()) %in% c("zip","gz","ZIP","GZ")) {
                meta <- as.data.frame(read_delim(getZip(MetaData_file_react()), delim = '\t', col_names = T))
                #meta <- as.data.frame(read_delim(getZip(MetaData_file_react()), delim = '\t', col_names = T))
              } else if (tools::file_ext(MetaData_file_react()) %in% c("csv","CSV")) {
                meta <- as.data.frame(read_delim(url(MetaData_file_react()), delim = ',', col_names = T))
              } else if (tools::file_ext(MetaData_file_react()) %in% c("RData","rdata")) {
                meta <- loadRData(url(MetaData_file_react()))
              } else if (tools::file_ext(MetaData_file_react()) %in% c("rds","RDS")) {
                meta <- readRDS(url(MetaData_file_react()))
              }
            }
            incProgress(0.5, detail = "Complete!")
          })
          user_upload(list(expr = expr,
                           meta = meta))
        } else if (isTruthy(ExpressionMatrix_file_react()) & isTruthy(MetaData_file_react()) & isTruthy(gse_id_react())) {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.5, detail = "Loading GSE Data")
            gse_id = as.character(gse_id_react())
            gset <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE)
            gse <- getGEO(gse_id, GSEMatrix = FALSE, getGPL = TRUE)
            user_upload(list(gset = gset,
                             gse = gse))
            incProgress(0.5, detail = "Complete!")
          })
        } else if (isTruthy(ExpressionMatrix_file_react()) & isTruthy(MetaData_file_react()) & isTruthy(recount_id_react())) {
          if (RecountAvail) {
            withProgress(message = "Processing", value = 0, {
              incProgress(0.5, detail = "Loading Recount Data")
              recount3_Projects <- available_projects()
              recount3_proj(recount3_Projects)
              RecountID <- recount_id_react()
              if (RecountID %in% recount3_Projects$project) {
                Proj_home <- recount3_Projects[which(recount3_Projects$project == RecountID),"project_home"]
                Proj_org <- recount3_Projects[which(recount3_Projects$project == RecountID),"organism"]
                Proj_source <- toupper(recount3_Projects[which(recount3_Projects$project == RecountID),"file_source"])
                File_prefix <- paste0(Proj_source,"_",RecountID)
                recount_rse <- recount3::create_rse_manual(
                  project = RecountID,
                  project_home = Proj_home,
                  organism = Proj_org,
                  annotation = "gencode_v26",
                  type = "gene",
                  verbose = TRUE
                )
                if (Proj_org == "mouse") {
                  MM_react(TRUE)
                } else {
                  MM_react(FALSE)
                }
                user_upload(recount_rse)
              }
              incProgress(0.5, detail = "Complete!")
            })
          }
        }
      })
      
      observe({
        if (isTruthy(recount_id_react())) {
          updateCheckboxGroupInput(session,"RawCountQuantNorm",choices = c("Input data is log-transformed","Normalize Raw Counts","Quantile Normalization","Filter Matrix"),
                                   selected = "Normalize Raw Counts")
        }
      })
      
      observe({
        req(input$HumanOrMouse)
        if (input$HumanOrMouse == "Mouse to Human Conv"){
          MM_react(FALSE)
        }
      })
      
      observe({
        FileCheckAlerts_list <- c()
        req(user_upload())
        withProgress(message = "Processing", value = 0, {
          incProgress(0.5, detail = "Formatting Input Data")
          if (isTruthy(ExpressionMatrix_file_react()) & isTruthy(MetaData_file_react()) && !isTruthy(gse_id_react()) && !isTruthy(recount_id_react())) {
            req(user_upload())
            expr <- user_upload()$expr
            meta <- user_upload()$meta
            expr_raw(expr)
          }
          if (isTruthy(ExpressionMatrix_file_react()) & isTruthy(MetaData_file_react()) & isTruthy(gse_id_react())) {
            req(user_upload())
            gset <- user_upload()$gset
            gse <- user_upload()$gse
            
            gse_id = as.character(gse_id_react())
            platform <- names(GPLList(gse))
            
            if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
            gset <- gset[[idx]]
            expr <- exprs(gset)
            df <- as.data.frame(expr)
            df <- arrange(df)
            df$ID <- rownames(df)
            df <- df %>% relocate(ID)
            df$ID <- trimws(df$ID)
            df$ID <- as.character(df$ID)
            full_table <- read_delim(paste("https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/GPL_Files/", platform, ".txt", sep=""), delim = '\t', col_names = T)
            full_table[1] <- NULL
            colnames(full_table)[1] = "ID"
            colnames(full_table)[2] = "Symbol"
            full_table <- arrange(full_table)
            full_table = data.frame(full_table)
            full_table$ID <- as.character(full_table$ID)
            # convert id to symbol for expr matrix
            expr <- inner_join(df, full_table, by = "ID")
            expr <- relocate(expr, "Symbol", .before = 1)
            expr$ID <- NULL
            expr <- expr[complete.cases(expr$Symbol), ]
            expr$Symbol <- trimws(expr$Symbol)
            expr_raw(expr)
            
            meta <- pData(gset)
            org <- meta[1,"organism_ch1"]
            meta <- relocate(meta, geo_accession, .before = title)
            # cleaning the meta data
            # Select columns that start with "characteristic"
            char_cols <- grep("^characteristic", colnames(meta), value = TRUE)
            extr_meta <- meta[, c("geo_accession", char_cols)]
            
            original_df <- extr_meta
            
            # Initialize an empty data frame to store results
            new_df <- data.frame(sample = character(0), title = character(0), value = character(0), stringsAsFactors = FALSE)
            
            # Iterate through rows and columns to extract values and create new rows
            for (i in seq_along(original_df$geo_accession)) {
              row_values <- as.character(original_df[i, -1])  # Convert to character
              row_values <- row_values[which(!row_values == "")]
              
              for (j in seq_along(row_values)) {
                split_value <- trimws(strsplit(row_values[j], ":")[[1]])
                title <- split_value[1]
                value <- paste(split_value[-1], collapse = ":")
                
                new_row <- data.frame(geo_accession = original_df$geo_accession[i], title = title, value = value, stringsAsFactors = FALSE)
                new_df <- rbind(new_df, new_row)
              }
            }
            
            final_extr_meta <- tidyr::pivot_wider(new_df, names_from = title, values_from = value)
            final_extr_meta <- as.data.frame(final_extr_meta)
            meta_merged <- dplyr::full_join(final_extr_meta, meta[,which(!colnames(meta) %in% char_cols)], by = "geo_accession")
            meta <- meta_merged
            if (tolower(org) == "mus musculus") {
              MM_react(TRUE)
            } else (
              MM_react(FALSE)
            )
          }
          if (isTruthy(ExpressionMatrix_file_react()) & isTruthy(MetaData_file_react()) & isTruthy(recount_id_react())) {
            if (RecountAvail) {
              req(user_upload())
              RecountID <- recount_id_react()
              if (RecountID %in% recount3_proj()$project) {
                Proj_home <- recount3_proj()[which(recount3_proj()$project == RecountID),"project_home"]
                Proj_org <- recount3_proj()[which(recount3_proj()$project == RecountID),"organism"]
                Proj_source <- toupper(recount3_proj()[which(recount3_proj()$project == RecountID),"file_source"])
                File_prefix <- paste0(Proj_source,"_",RecountID)
                
                recount_rse <- user_upload()
                expr <- assay(recount_rse)
                expr_raw(expr)
                if ("Normalize Raw Counts" %in% input$RawCountQuantNorm) {
                  if (isTruthy(input$RawCountNorm)) {
                    if (input$RawCountNorm %in% c("TPM","RPKM")) {
                      if (input$RawCountNorm == "TPM") {
                        recount_rse <- user_upload()
                        assays(recount_rse)$counts <- transform_counts(recount_rse)
                        expr <- as.data.frame(recount::getTPM(recount_rse))
                        message <- paste0("Raw counts normalized by TPM")
                        FileCheckAlerts_list <- c(FileCheckAlerts_list,message)
                      } else if (input$RawCountNorm == "RPKM") {
                        recount_rse <- user_upload()
                        assays(recount_rse)$counts <- transform_counts(recount_rse)
                        expr <- as.data.frame(recount::getRPKM(recount_rse))
                        message <- paste0("Raw counts normalized by RPKM")
                        FileCheckAlerts_list <- c(FileCheckAlerts_list,message)
                      }  else if (input$RawCountNorm == "none") {
                        recount_rse <- user_upload()
                        expr <- assay(recount_rse)
                      }
                    }
                  }
                }
                
                expr_anno <- as.data.frame(recount_rse@rowRanges@elementMetadata@listData)
                expr <- merge(expr_anno[,c(5,7)],expr,by.x = "gene_id", by.y = 0,all.y = T)
                expr <- expr[,-1]
                
                
                meta <- as.data.frame(t(do.call(rbind,recount_rse@colData@listData)))
                
                if (Proj_source == "TCGA") {
                  names(expr)[-1][match(meta[,"external_id"], names(expr)[-1])] = meta[,"tcga.tcga_barcode"]
                  meta <- meta %>% relocate(tcga.tcga_barcode)
                } else {
                  meta <- meta %>% relocate("external_id")
                }
              }
            }
          }
          
          SpecDetect <- detect_species(expr[,1])
          if (SpecDetect == "mouse") {
            updateRadioButtons(session,"HumanOrMouse",NULL,c("Human","Mouse","Mouse to Human Conv"), selected = "Mouse")
            MM_react(TRUE)
            FileCheckAlerts_list <- c(FileCheckAlerts_list,paste0("Mouse gene symbols detected."))
          }
          
          colnames(expr)[1] <- "Gene"
          expr$Gene <- date_to_gene(expr$Gene)
          if (isTruthy(MM_react())) {
            if (as.logical(MM_react())) {
              expr$Gene <- date_to_gene(expr$Gene, mouse = T)
            }
          }
          
          
          if (isTruthy(input$HumanOrMouse)) {
            if (input$HumanOrMouse == "Mouse to Human Conv") {
              conv_df <- as.data.frame(fread(MM_HS_Conversion_File))
              colnames(conv_df)[1] <- colnames(expr)[1]
              expr_new <- merge(conv_df,expr)
              message1 <- "Mouse gene symbols converted to human"
              FileCheckAlerts_list <- c(FileCheckAlerts_list,message1)
              if (nrow(expr) == nrow(expr_new)) {
                message2 <- "All genes converted successfully"
                FileCheckAlerts_list <- c(FileCheckAlerts_list,message2)
              } else {
                message2 <- paste0(abs(nrow(expr)-nrow(expr_new))," genes unable to be converted.")
                FileCheckAlerts_list <- c(FileCheckAlerts_list,message2)
              }
              expr <- expr_new
            }
            
          }
          
          
          
          if (length(expr[sapply(expr, is.infinite)]) > 0) {
            message <- paste0(length(expr[sapply(expr, is.infinite)]), " infinite values found. Replaced with NA." )
            FileCheckAlerts_list <- c(FileCheckAlerts_list,message)
            expr[sapply(expr, is.infinite)] <- NA
          }
          
          # Remove Expression with NA
          exprRows <- nrow(expr)
          expr <- expr %>%
            drop_na()
          if (exprRows > nrow(expr)) {
            message <- paste0(exprRows-nrow(expr), " features with NA values found. Features removed." )
            FileCheckAlerts_list <- c(FileCheckAlerts_list,message)
          }
          # Check that expression data is numeric
          isChar <- unname(which(sapply(expr, function(x) is.character(x))))
          isChar <-  isChar[-1]
          if (length(isChar) > 0) {
            expr[isChar] <- sapply(expr[isChar],as.numeric)
          }
          # Remove Duplicate genes
          expr_dup <- expr[which(expr[,1] %in% expr[,1][duplicated(expr[,1])]),]
          expr_nondup <- expr[which(!expr[,1] %in% expr[,1][duplicated(expr[,1])]),]
          if (nrow(expr_dup) > 0) {
            expr_dup <- expr_dup %>%
              group_by(Gene) %>%
              summarise_all(max)
          }
          expr <- rbind(expr_dup,expr_nondup)
          if (nrow(expr_dup) > 0) {
            message <- paste0(length(unique(expr_dup[,1])), " duplicate features found. Features reduced to those with maximum value." )
            FileCheckAlerts_list <- c(FileCheckAlerts_list,message)
          }
          expr <- as.data.frame(expr)
          # Make rownames <- genenames
          row.names(expr) <- expr[,1]
          expr <- expr[,-1]
          expr = expr[order(row.names(expr)), ]
          expr_col <- colnames(expr)
          
          if ("Input data is log-transformed" %in% input$RawCountQuantNorm) {
            expr <- 2^as.matrix(expr)
            message <- paste0("Features exponentiated.")
            FileCheckAlerts_list <- c(FileCheckAlerts_list,message)
          }
          if ("Normalize Raw Counts" %in% input$RawCountQuantNorm) {
            if (isTruthy(input$RawCountNorm)) {
              if (input$RawCountNorm %in% c("TMM","upperquartile")) {
                mat <- as.matrix(expr)
                mat_dgeList <- DGEList(counts = as.matrix(mat))
                mat_dgeList_Norm <- edgeR::calcNormFactors(mat_dgeList, method = input$RawCountNorm)
                expr <- edgeR::cpm(mat_dgeList_Norm)
                message <- paste0(input$RawCountNorm," normalization preformed")
                FileCheckAlerts_list <- c(FileCheckAlerts_list,message)
              } else if (input$RawCountNorm == "DESeq2") {
                if (isTruthy(input$DESeqDesignCol)) {
                  desCol <- input$DESeqDesignCol
                  mat <- as.matrix(expr)
                  sortedIndices <- match(colnames(mat), meta[,1])
                  meta_deseq <- meta[sortedIndices,,drop=FALSE]
                  rownames(meta_deseq) <- meta_deseq[,1]
                  if (desCol != "Select column") {
                    desRef <- input$DESeqDesignColRef
                    if (isTruthy(desRef)) {
                      meta_deseq[,desCol] <- relevel(factor(meta_deseq[,desCol]), ref = desRef)
                      colnames(meta_deseq) <- gsub(" ","_",colnames(meta_deseq))
                      colnames(meta_deseq) <- gsub("[[:punct:]]","_",colnames(meta_deseq))
                      desCol <- gsub(" ","_",desCol)
                      desCol <- gsub("[[:punct:]]","_",desCol)
                      dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat,
                                                            colData = meta_deseq,
                                                            design = as.formula(paste0("~ ",desCol)))
                      dds <- DESeq(dds)
                      res <- results(dds)
                      resOrdered <- res[order(res$padj),]
                      normalizedCounts <- counts(dds, normalized=TRUE)
                      expr <- as.matrix(normalizedCounts)
                      message <- paste0("DESeq2 Normalization performed using ",desCol," as the design column and ",desRef," as the reference")
                      FileCheckAlerts_list <- c(FileCheckAlerts_list,message)
                    }
                  }
                }
              }
            }
          }
          if ("Quantile Normalization" %in% input$RawCountQuantNorm) {
            expr_rows <- rownames(expr)
            expr_cols <- colnames(expr)
            expr <- preprocessCore::normalize.quantiles(as.matrix(expr))
            rownames(expr) <- expr_rows
            colnames(expr) <- expr_cols
            message <- paste0("Quantile Normalization preformed")
            FileCheckAlerts_list <- c(FileCheckAlerts_list,message)
          }
          if ("Filter Matrix" %in% input$RawCountQuantNorm & isTruthy(input$FilterNum) & isTruthy(input$FilterProp)) {
            FilterProp <- input$FilterProp/100
            expr_pass <- apply(expr,1,function(x) ExprFilter2(x, input$FilterNum, FilterProp))
            expr <- expr[which(expr_pass == TRUE),]
            message <- paste0("Features with an expression less than ",input$FilterNum," in at least ",
                              input$FilterProp,"% of Samples filtered out. ",
                              length(which(expr_pass == FALSE))," features removed.")
            FileCheckAlerts_list <- c(FileCheckAlerts_list,message)
          }
          
          #gene list file from expression data
          Gene <- rownames(expr)
          geneList <- as.data.frame(Gene)
          geneList_raw(geneList);
          Gene_raw(Gene)
          
          #meta[,1] <- gsub("[_.-]", "_", meta[,1])
          colnames(meta)[1] <- "SampleName"
          
          #for heatmap sample selection
          sampsames <- intersect(colnames(expr),meta[,1])
          
          if (length(sampsames) != ncol(expr) | length(sampsames) != nrow(meta)) {
            message <- paste0("Mistmatching or missing sample names found between feature matrix and meta data. Reduced to only similar samples (N=",length(sampsames),")")
            FileCheckAlerts_list <- c(FileCheckAlerts_list,message)
          }
          
          #ensure expression samples and meta are exact
          expr <- expr[,sampsames]
          meta <- meta[which(meta[,1] %in% sampsames),]
          
          #if (immudecon_check == TRUE & "Perform Immune Deconvolution" %in% input$RawCountQuantNorm) {
          #  withProgress(message = "Processing", value = 0, {
          #    incProgress(0.5, detail = "Performing Immune Deconvolution")
          #    mcp_counter_decon <- as.data.frame(deconvolute(expr, "mcp_counter"))
          #    rownames(mcp_counter_decon) <- paste0(mcp_counter_decon[,1],"_MCP_Counter_Immunedeconv")
          #    estimate_decon <- as.data.frame(deconvolute(expr, "estimate"))
          #    rownames(estimate_decon) <- paste0(estimate_decon[,1],"_Estimate_Immunedeconv")
          #    imm_deconv <- rbind(mcp_counter_decon,estimate_decon)[,-1]
          #    ImmDeconv_react(imm_deconv)
          #    imm_deconv <- as.data.frame(t(imm_deconv))
          #    meta <- merge(meta,imm_deconv, by.x = colnames(meta)[1], by.y = 0, all.x = T)
          #    incProgress(0.5, detail = "Complete")
          #  })
          #}
          
          
          expr_input(expr)
          meta_input(meta)
          
          
          
          if (isTruthy(MM_react())) {
            if (as.logical(MM_react())) {
              CTKgenes <- CTKgenes_mm
              updateRadioButtons(session,"HumanOrMouse",selected = "Mouse")
            } else {
              CTKgenes <- CTKgenes_hs
            }
          } else {
            CTKgenes <- CTKgenes_hs
          }
          
          CTKgenes <- CTKgenes[which(CTKgenes %in% Gene)]
          updateSelectizeInput(session = session, inputId = "heatmapGeneSelec",
                               choices = Gene,selected = CTKgenes, server = T)
          updateSelectizeInput(session = session, inputId = "userheatsamp2",
                               choices = sampsames,selected = sampsames, server = T)
          updateSelectizeInput(session = session, inputId = "avgheatmapGeneSelec",
                               choices = Gene,selected = CTKgenes, server = T)
          updateSelectizeInput(session = session, inputId = "scatterG1",
                               choices = Gene,selected = Gene[1], server = T)
          updateSelectizeInput(session = session, inputId = "scatterG2",
                               choices = Gene,selected = Gene[2], server = T)
          updateSelectizeInput(session = session, inputId = "userGeneSelec",
                               choices = sort(as.vector(geneList[,1])),selected = NULL, server = T)
          updateSelectizeInput(session = session, inputId = "scatterGeneSelec",
                               choices = rownames(expr),selected = NULL, server = T)
          updateSelectizeInput(session = session, inputId = "userheatsampSS",
                               choices = sampsames,selected = sampsames, server = T)
          
          FileCheckAlerts_react(FileCheckAlerts_list)
          incProgress(0.5, detail = "Complete!")
        })
      })
      
      
      observeEvent(input$SplitColumn, {
        req(meta_input())
        meta <- meta_input()
        meta_col_names_start <- colnames(meta)
        Split_col_name <- input$ColToSplit
        Split_by <- input$SplitDelim
        New_col_names_in <- input$SplitNewColNames
        New_col_names <- strsplit(New_col_names_in,",")[[1]]
        if (isTruthy(Split_col_name) & isTruthy(Split_by) & length(New_col_names) > 1) {
          meta <- meta %>%
            separate_wider_delim(!!sym(Split_col_name),Split_by,names = New_col_names,cols_remove = FALSE,
                                 too_few = "align_start", too_many = "merge", names_repair = "unique")
          meta_col_names_end <- colnames(meta)
          meta <- meta %>%
            relocate(any_of(c(Split_col_name,setdiff(meta_col_names_end,meta_col_names_start))), .after = 1) %>%
            as.data.frame()
          meta_input(meta)
        }
      })
      observeEvent(input$MergeColumns, {
        req(meta_input())
        meta <- meta_input()
        Merg_col_name <- input$ColsToMerge
        merge_with <- input$MergeDelim
        New_col_name <- input$MergeNewColName
        if (length(Merg_col_name) > 1 & isTruthy(merge_with) & isTruthy(New_col_name)) {
          meta[,New_col_name] <- apply(meta[,Merg_col_name],1,function(x) {paste(x,collapse = merge_with)})
          meta <- meta %>% relocate(any_of(c(Merg_col_name,New_col_name)),.after = 1) %>%
            as.data.frame()
          meta_input(meta)
        }
      })
      
      
      ## Data Preview ----------------------------------------------------------
      
      observe({
        if (immudecon_check) {
          if ("Perform Immune Deconvolution" %in% input$RawCountQuantNorm) {
            req(expr_input())
            expr <- expr_input()
            withProgress(message = "Processing", value = 0, {
              incProgress(0.5, detail = "Performing Immune Deconvolution")
              mcp_counter_decon <- as.data.frame(deconvolute(expr, "mcp_counter"))
              rownames(mcp_counter_decon) <- paste0(mcp_counter_decon[,1],"_MCP_Counter_Immunedeconv")
              estimate_decon <- as.data.frame(deconvolute(expr, "estimate"))
              rownames(estimate_decon) <- paste0(estimate_decon[,1],"_Estimate_Immunedeconv")
              imm_deconv <- rbind(mcp_counter_decon,estimate_decon)[,-1]
              ImmDeconv_react(imm_deconv)
              incProgress(0.5, detail = "Complete")
            })
          }
        }
      })
      
      observe({
        if (immudecon_check) {
          req(input$ImmDeconvAnalysis)
          if (input$ImmDeconvAnalysis == "Immune Deconvolution") {
            req(ImmDeconv_react())
            expr <- ImmDeconv_react()
            Gene <- rownames(expr)
            geneList <- as.data.frame(Gene)
            geneList_raw(geneList);
            Gene_raw(Gene)
            CTKgenes <- CTKgenes[which(CTKgenes %in% Gene)]
            updateSelectizeInput(session = session, inputId = "heatmapGeneSelec",
                                 choices = Gene,selected = CTKgenes, server = T)
            updateSelectizeInput(session = session, inputId = "avgheatmapGeneSelec",
                                 choices = Gene,selected = CTKgenes, server = T)
            updateSelectizeInput(session = session, inputId = "scatterG1",
                                 choices = Gene,selected = Gene[1], server = T)
            updateSelectizeInput(session = session, inputId = "scatterG2",
                                 choices = Gene,selected = Gene[2], server = T)
            updateSelectizeInput(session = session, inputId = "userGeneSelec",
                                 choices = sort(as.vector(geneList[,1])),selected = NULL, server = T)
            updateSelectizeInput(session = session, inputId = "scatterGeneSelec",
                                 choices = rownames(expr),selected = NULL, server = T)
          } else if (input$ImmDeconvAnalysis == "Gene Expression") {
            req(ImmDeconv_react())
            expr <- expr_input()
            Gene <- rownames(expr)
            geneList <- as.data.frame(Gene)
            geneList_raw(geneList);
            Gene_raw(Gene)
            CTKgenes <- CTKgenes[which(CTKgenes %in% Gene)]
            updateSelectizeInput(session = session, inputId = "heatmapGeneSelec",
                                 choices = Gene,selected = CTKgenes, server = T)
            updateSelectizeInput(session = session, inputId = "avgheatmapGeneSelec",
                                 choices = Gene,selected = CTKgenes, server = T)
            updateSelectizeInput(session = session, inputId = "scatterG1",
                                 choices = Gene,selected = Gene[1], server = T)
            updateSelectizeInput(session = session, inputId = "scatterG2",
                                 choices = Gene,selected = Gene[2], server = T)
            updateSelectizeInput(session = session, inputId = "userGeneSelec",
                                 choices = sort(as.vector(geneList[,1])),selected = NULL, server = T)
            updateSelectizeInput(session = session, inputId = "scatterGeneSelec",
                                 choices = rownames(expr),selected = NULL, server = T)
          }
        }
        
        
      })
      
      output$FileCheckAlerts <- renderPrint({
        
        req(FileCheckAlerts_react())
        text <- paste(FileCheckAlerts_react(), collapse = "\n")
        cat(text)
        
      })
      
      output$ExprFile_Preview <- DT::renderDataTable({
        req(expr_input())
        req(input$ExprHead)
        
        if (immudecon_check & input$ImmDeconvAnalysis == "Immune Deconvolution") {
          req(ImmDeconv_react())
          expr <- ImmDeconv_react()
        } else {
          expr <- expr_input()
        }
        if (input$ExprHead == "View table head") {
          expr <- head(expr,c(100,100))
        }
        DT::datatable(expr,
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 10,
                                     scrollX = T),
                      rownames = T) %>%
          formatRound(columns = colnames(expr), digits = 4)
      })
      
      output$ClinFile_Preview <- DT::renderDataTable({
        req(meta_react())
        req(input$ClinHead)
        clin <- meta_react()
        if (input$ClinHead == "View table head") {
          clin <- head(clin,c(100,100))
        }
        DT::datatable(clin,
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 10,
                                     scrollX = T),
                      rownames = F)
      })
      
      output$renddownload_expr <- renderUI({
        if (isTruthy(expr_input()) | isTruthy(meta_input())) {
          downloadButton("download_expr","Download Expression")
        }
      })
      output$renddownload_meta <- renderUI({
        if (isTruthy(expr_input()) | isTruthy(meta_input())) {
          downloadButton("download_meta","Download Meta")
        }
      })
      output$renddownload_notes <- renderUI({
        if (length(FileCheckAlerts_react()) > 0) {
          downloadButton("download_notes","Download data processing notes")
        }
      })
      
      # Download handler for the button
      output$download_expr <- downloadHandler(
        filename = function() {
          paste(ProjectName_react(),"_expr_matrix_", ExpressionMatrix_file_react(), "_", Sys.Date(), ".txt", sep = "")
        },
        content = function(file) {
          expr <- as.data.frame(expr_react())
          expr <- cbind(rownames(expr),expr)
          colnames(expr)[1] <- "Gene"
          write.table(expr, file, sep = '\t', row.names = F)
        }
      )
      
      output$download_meta <- downloadHandler(
        filename = function() {
          paste(ProjectName_react(),"_meta_data_", MetaData_file_react(), "_", Sys.Date(), ".txt", sep = "")
        },
        content = function(file) {
          write.table(meta_react(), file, sep = '\t', row.names = F)
        }
      )
      
      output$download_notes <- downloadHandler(
        filename = function() {
          paste(ProjectName_react(), "_DataProcessingNotes_", Sys.Date(), ".txt", sep = "")
        },
        content = function(file) {
          df <- data.frame(Notes = FileCheckAlerts_react())
          write_tsv(df, file, col_names = F)
        }
      )
      
      
      # Gene Set Data -------------------------------------------------------------
      #gmt <- reactiveVal()
      #msigdb.gsea2 <- reactiveVal()
      #tab2 <- reactiveVal()
      #GeneSet2 <- reactiveVal()
      #gs <- reactiveVal()
      #gs2 <- reactiveVal()
      
      geneset_list <- reactiveVal()
      geneset_cat_tab <- reactiveVal()
      
      observe({
        
        if (isTruthy(MM_react())) {
          if (as.logical(MM_react())) {
            print("Reading mouse gene sets")
            geneset_mm <- loadRData(GeneSetMM_File)
            gs_cat_mm <- as.data.frame(fread(GeneSetTableMM_File))
            updateSelectInput(session,"GeneSetCat_Select",choices = unique(gs_cat_mm[,1]))
            updateSelectInput(session,"GeneSetCat_Select_heat",choices = c(unique(gs_cat_mm[,1]),"User Input"),
                              selected = c(unique(gs_cat_mm[,1]),"User Input")[1])
            geneset_list(geneset_mm)
            geneset_cat_tab(gs_cat_mm)
            #write in the name of your gene set list for shiny UI
            #userGSlist_name <- 'Cell Marker'
            #path to your gene set file .gmt or .txt/.tsv
            #userGS_file <- 'Genesets/CellMarker_gsNsym_MM.tsv'
            #path to your R data list object for ssGSEA
            #userRData_file <- 'Genesets/CellMarker_GS_MM.RData'
            #MSigDB gene set
            #msigdb <- 'Genesets/msigdb_gsNsym_MM.zip'
            #MSigDB gene set FOR UI
            #msigdb2 <- 'Genesets/msigdb_gsNcat_MM.tsv'
            #gene set list for ssGSEA
            #gs_file <- 'Genesets/msigdb_gs_MM.RData'
            #Cytokine genes for mouse
            CTKgenes <- c("Il2","Il12a","Il12b","Il17a","Ifna13","Ifnb1","Ifng","Ifngr1","Cd11b","Itgam",
                          "Cd33","Entpd1","Icosl","Icos","Tnfsf9","Tnfrsf9","Cd40","Cd40lg","Cd70","Cd27",
                          "Tnfsf18","Tnfrsf18","Tnfsf14","Tnfrsf14","Tnfsf4","Tnfrsf4","H2-K1","CD3G",
                          "Ceacam1","Cd80","Cd86","Ctla4","Cd276","Vtcn1","Pvr","Cd226","Tigit","Cd96","Lgals3",
                          "Lgals3bp","Lgals9","Lgals9c","Havcr2","Hhla2","Cd274","Pdcd1lg2","Pdcd1","Vsir")
          } else {
            print("Reading human gene sets")
            geneset_hs <- loadRData(GeneSetHS_File)
            gs_cat_hs <- as.data.frame(fread(GeneSetTableHS_File))
            updateSelectInput(session,"GeneSetCat_Select",choices = unique(gs_cat_hs[,1]))
            updateSelectInput(session,"GeneSetCat_Select_heat",choices = c(unique(gs_cat_hs[,1]),"User Input"),
                              selected = c(unique(gs_cat_hs[,1]),"User Input")[1])
            geneset_list(geneset_hs)
            geneset_cat_tab(gs_cat_hs)
            #write in the name of your gene set list for shiny UI
            #userGSlist_name <- 'LINCS L1000'
            #path to your gene set file .gmt or .txt/.tsv
            #userGS_file <- 'Genesets/CellMarker_gsNsym_HS_v2.txt'
            #path to your R data list object for ssGSEA
            #userRData_file <- 'Genesets/CellMarker_GS_HS_v2.RData'
            #MSigDB gene set
            #msigdb <- 'Genesets/msigdb_gsNsym_HS_v2.zip'
            #MSigDB gene set FOR UI
            #msigdb2 <- 'Genesets/msigdb_gsNcat_HS_v2.txt'
            #gene set list for ssGSEA
            gs_file <- 'Genesets/msigdb_gs_HS_v2.RData'
            CTKgenes <- c("IL2","IL12A","IL12B","IL17A","IFNA1","IFNB1","IFNG","IFNGR","CD11b",
                          "ITGAM","CD33","ENTPD1","ICOSLG","CD275","CD278","TNFSF9","TNFRSF9",
                          "CD40","CD40LG","CD70","CD27","TNFSF18","TNFRSF18","TNFSF14","TNFRSF14",
                          "TNFSF4","TNFRSF4","HLA-A","CD3","CEACAM1","CD80","CD86","CTLA4","CD276",
                          "VTCN1","PVR","CD226","TIGIT","CD96","LGALS3","LGALS3BP","LGALS9","LGALS9C",
                          "HAVCR2","HHLA2","TMIGD2","CD274","PDCD1LG2","PDCD1","VSIR")
          }
        } else {
          print("Reading human gene sets")
          geneset_hs <- loadRData(GeneSetHS_File)
          gs_cat_hs <- as.data.frame(fread(GeneSetTableHS_File))
          updateSelectInput(session,"GeneSetCat_Select",choices = unique(gs_cat_hs[,1]))
          updateSelectInput(session,"GeneSetCat_Select_heat",choices = c(unique(gs_cat_hs[,1]),"User Input"),
                            selected = c(unique(gs_cat_hs[,1]),"User Input")[1])
          geneset_list(geneset_hs)
          geneset_cat_tab(gs_cat_hs)
          #write in the name of your gene set list for shiny UI
          #userGSlist_name <- 'LINCS L1000'
          #path to your gene set file .gmt or .txt/.tsv
          #userGS_file <- 'Genesets/CellMarker_gsNsym_HS_v2.txt'
          #path to your R data list object for ssGSEA
          #userRData_file <- 'Genesets/CellMarker_GS_HS_v2.RData'
          #MSigDB gene set
          #msigdb <- 'Genesets/msigdb_gsNsym_HS_v2.zip'
          #MSigDB gene set FOR UI
          #msigdb2 <- 'Genesets/msigdb_gsNcat_HS_v2.txt'
          #gene set list for ssGSEA
          #gs_file <- 'Genesets/msigdb_gs_HS_v2.RData'
          CTKgenes <- c("IL2","IL12A","IL12B","IL17A","IFNA1","IFNB1","IFNG","IFNGR","CD11b",
                        "ITGAM","CD33","ENTPD1","ICOSLG","CD275","CD278","TNFSF9","TNFRSF9",
                        "CD40","CD40LG","CD70","CD27","TNFSF18","TNFRSF18","TNFSF14","TNFRSF14",
                        "TNFSF4","TNFRSF4","HLA-A","CD3","CEACAM1","CD80","CD86","CTLA4","CD276",
                        "VTCN1","PVR","CD226","TIGIT","CD96","LGALS3","LGALS3BP","LGALS9","LGALS9C",
                        "HAVCR2","HHLA2","TMIGD2","CD274","PDCD1LG2","PDCD1","VSIR")
        }
        
        
        
        #MSigDB gene sets
        #msigdb.gsea <- as.data.frame(read_delim(msigdb, delim = '\t'))
        #msigdb.gsea <- as.data.frame(fread(msigdb, sep = '\t'))
        #gmt <- msigdb.gsea
        #MSigDB gene sets FOR UI
        #msigdb.gsea2 <- as.data.frame(read_delim(msigdb2, delim = '\t'))
        #msigdb.gsea2 <- as.data.frame(fread(msigdb2, sep = '\t'))
        #tab2 User gene set
        #tab2 <- as.data.frame(read_delim(userGS_file, col_names = T, delim = '\t'))
        #tab2 <- as.data.frame(fread(userGS_file, header = T, sep = '\t'))
        #tab2 back end
        #GeneSet2 <- as.data.frame(unique(tab2[,1]))
        #rownames(GeneSet2) <- 1:nrow(GeneSet2)
        #colnames(GeneSet2)[1] <- "Gene_Set"
        #tab2 R Data list
        #gs <- loadRData(gs_file)
        #gs2 <- loadRData(userRData_file)
        
        #gmt(gmt)
        #msigdb.gsea2(msigdb.gsea2)
        #tab2(tab2)
        #GeneSet2(GeneSet2)
        #gs(gs)
        #gs2(gs2)
        
      })
      
      GeneSetTable_React <- reactive({
        GS_database <- input$GeneSetCat_Select
        GeneSetTable <- geneset_cat_tab()
        sub_tab <- GeneSetTable[which(GeneSetTable[,1] == GS_database),]
        new_tab <- sub_tab[,-1]
        new_tab
      })
      
      observe({
        
        req(input$GeneSetCat_Select_heat)
        ssgsea_heat_gs_cat <- input$GeneSetCat_Select_heat
        GeneSetTable <- geneset_cat_tab()
        gs_choices <- c()
        if (any(ssgsea_heat_gs_cat %in% GeneSetTable[,1])) {
          gs_choices <- GeneSetTable[which(GeneSetTable[,1] %in% ssgsea_heat_gs_cat),ncol(GeneSetTable)]
        }
        if ("User Input" %in% ssgsea_heat_gs_cat) {
          gs_choices <- c(userGeneSet_table()[,1],gs_choices)
        }
        updateSelectizeInput(session,"ssgseaHeatGS", choices = gs_choices, selected = gs_choices[1:3], server = T)
        
      })
      
      ## Render Gene Set Selection Table
      output$GeneSetTable <- DT::renderDataTable({
        if (input$GeneSetMultiOpt == "Select single gene set") {
          DT::datatable(GeneSetTable_React(),
                        options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                       pageLength = 10,
                                       scrollX = T),
                        selection = list(mode = 'single', selected = 1),
                        rownames = F)
        } else {
          DT::datatable(GeneSetTable_React(),
                        options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                       pageLength = 10,
                                       scrollX = T),
                        selection = list(selected = 1),
                        rownames = F)
        }
      })
      
      ## Render User Gene Set Table - Backend
      userGeneSet_loaded <- reactive({
        gs.u <- input$userGeneSet
        ext <- tools::file_ext(gs.u$datapath)
        req(gs.u)
        #validate(need(ext == c("gmt","tsv","txt", "RData"), "Please upload .gmt, .tsv, .txt, or .RData file"))
        # If user provides GMT file
        if (ext == "gmt") {
          gmt <- clusterProfiler::read.gmt(gs.u$datapath)
        } else if (ext == "RData") {
          gmt <- loadRData(gs.u$datapath)
        } else {
          #gmt <- as.data.frame(read_delim(gs.u$datapath, delim = '\t'))
          gmt <- as.data.frame(fread(gs.u$datapath))
        }
        gmt
      })
      
      userGeneSet_table <- reactive({
        gs.u <- input$userGeneSet
        ext <- tools::file_ext(gs.u$datapath)
        gmt <- userGeneSet_loaded()
        # If user provides GMT file
        if (ext == "gmt") {
          uGS_table <- as.data.frame(unique(gmt[,1]))
          colnames(uGS_table)[1] <- "GeneSet"
        } else if (ext == "RData") {
          uGS_table <- as.data.frame(names(gmt))
          colnames(uGS_table)[1] <- "GeneSet"
        } else {
          uGS_table <- as.data.frame(unique(gmt[,1]))
          colnames(uGS_table)[1] <- "GeneSet"
        }
        uGS_table
      })
      
      #observe({
      #  
      #  #req(userGeneSet_table())
      #  req(geneset_list())
      #  gs_list <- geneset_list()
      #  if (isTruthy(userGeneSet_table())) {
      #    gs_df <- userGeneSet_table()
      #    gs_list <- c(gs_df[,1], names(gs_list))
      #  }
      #  updateSelectizeInput(session,"ssgseaHeatGS", choices = gs_list, selected = gs_list[1:3], server = T)
      #  
      #  
      #})
      
      output$userGeneSetTable <- DT::renderDataTable({
        req(input$userGeneSet)
        uGS_table <- userGeneSet_table()
        DT::datatable(uGS_table,
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 10,
                                     scrollX = T),
                      selection = list(mode = 'single', selected = 1),
                      rownames = F)
      })
      
      
      
      gs_react_heat <- reactive({
        
        req(geneset_list())
        req(input$ssgseaHeatGS)
        heat_gs <- input$ssgseaHeatGS
        gs <- geneset_list()
        gs2 <- gs[heat_gs]
        
        if ("User Input" %in% input$GeneSetCat_Select_heat) {
          gmt <- userGeneSet_loaded()
          gs.u <- input$userGeneSet
          ext <- tools::file_ext(gs.u$datapath)
          if (ext == "gmt") {
            if (any(heat_gs %in% gmt[,1])) {
              gmt[which(gmt[,1] %in% heat_gs),]
              geneset <- list()
              for (i in unique(gmt[,1])){
                geneset[[i]] <- gmt[gmt[,1] == i,][,1]
              }
              gs2 <- c(gs2,geneset)
            }
          } else if (ext == "RData") {
            if (any(heat_gs %in% names(gmt))) {
              geneset <- gmt[heat_gs]
              gs2 <- c(gs2,geneset)
            }
          } else {
            if (any(heat_gs %in% gmt[,1])) {
              gmt[which(gmt[,1] %in% heat_gs),]
              geneset <- list()
              for (i in unique(gmt[,1])){
                geneset[[i]] <- gmt[gmt[,1] == i,][,1]
              }
              gs2 <- c(gs2,geneset)
            }
          }
        }
        
        gs2
      })
      
      gs_react <- reactive({
        req(input$tables)
        gs <- geneset_list()
        if (input$tables == 1) {
          req(GeneSetTable_React())
          geneset_name <- GeneSetTable_React()[input$GeneSetTable_rows_selected,ncol(GeneSetTable_React())]
          geneset <- gs[geneset_name]
          geneset
        } else if (input$tables == 5) {
          req(input$userGeneSet)
          gs.u <- input$userGeneSet
          ext <- tools::file_ext(gs.u$datapath)
          gmt <- userGeneSet_loaded()
          uGS_table <- userGeneSet_table()
          geneset_name <- uGS_table[input$userGeneSetTable_rows_selected,1]
          if (isTruthy(geneset_name)) {
            if (ext == "gmt") {
              geneset <- list(geneset_name = gmt[which(gmt[,1] == geneset_name),2])
              names(geneset) <- geneset_name
            } else if (ext == "RData") {
              geneset <- gmt[geneset_name]
            } else {
              geneset <- list(geneset_name = gmt[which(gmt[,1] == geneset_name),2])
              names(geneset) <- geneset_name
            }
            geneset
          }
        }
      })
      
      geneset_gmt <- reactive({
        
        req(input$tables)
        gs <- geneset_list()
        GS_database <- input$GeneSetCat_Select
        GeneSetTable <- geneset_cat_tab()
        sub_tab <- GeneSetTable[which(GeneSetTable[,1] == GS_database),]
        gs <- gs[sub_tab[,ncol(sub_tab)]]
        withProgress(message = "Processing", value = 0, {
          incProgress(0.25, detail = "Generating GMT table of geneset")
          if (input$tables == 1) {
            geneset <- gs
            term <- rep(names(geneset), times = sapply(geneset, length))
            gene <- unlist(geneset)
            geneset <- data.frame(term = term, gene = gene)
            geneset
          } else if (input$tables == 5) {
            req(input$userGeneSet)
            gs.u <- input$userGeneSet
            ext <- tools::file_ext(gs.u$datapath)
            geneset <- userGeneSet_loaded()
            if (isTruthy(geneset)) {
              if (ext == "RData") {
                term <- rep(names(geneset), times = sapply(geneset, length))
                gene <- unlist(geneset)
                geneset <- data.frame(term = term, gene = gene)
              } else {
                geneset
              }
              geneset
            }
          }
          incProgress(0.5, detail = "Complete!")
        })
        geneset
        
      })
      
      
      gmt_react <- reactive({
        req(input$tables)
        gs <- geneset_list()
        if (input$tables == 1) {
          req(GeneSetTable_React())
          geneset_name <- GeneSetTable_React()[input$GeneSetTable_rows_selected,ncol(GeneSetTable_React())]
          geneset <- gs[geneset_name]
          term <- rep(names(geneset), times = sapply(geneset, length))
          gene <- unlist(geneset)
          geneset <- data.frame(term = term, gene = gene)
          #geneset <- data.frame(term = names(geneset),
          #                      gene = geneset[[1]])
          geneset
        } else if (input$tables == 5) {
          req(input$userGeneSet)
          gs.u <- input$userGeneSet
          ext <- tools::file_ext(gs.u$datapath)
          geneset <- userGeneSet_loaded()
          uGS_table <- userGeneSet_table()
          geneset_name <- uGS_table[input$userGeneSetTable_rows_selected,1]
          if (isTruthy(geneset_name)) {
            if (ext == "gmt") {
              geneset <- geneset[which(geneset[,1] %in% geneset_name),]
            } else if (ext == "RData") {
              geneset <- geneset[geneset_name]
              term <- rep(names(geneset), times = sapply(geneset, length))
              gene <- unlist(geneset)
              geneset <- data.frame(term = term, gene = gene)
              #geneset <- data.frame(term = names(geneset),
              #                      gene = geneset[[1]])
            } else {
              geneset <- geneset[which(geneset[,1] %in% geneset_name),]
            }
            geneset
          }
        }
      })
      
      
      
      
      
      
      # Render UI ------------------------------------------------------------------
      
      ## Data Exploration ----------------------------------------------------------
      observe({
        
        req(meta_input())
        meta <- meta_input()
        if (ncol(meta) > 2) {
          #CharCols <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
          CharCols <- c("Select All Samples",get_feat_cols(meta, keep_num = FALSE))
          if (isTruthy(SubCol_reactVal())) {
            preSelect <- SubCol_reactVal()
          } else {
            preSelect <- "Select All Samples"
          }
          updateSelectizeInput(session,"SubsetCol","Subset Samples By:",
                               choices = CharCols,
                               selected = preSelect, server = T)
        }
        
      })
      #output$rendSubsetCol <- renderUI({
      #  
      #  req(meta_input())
      #  meta <- meta_input()
      #  if (ncol(meta) > 2) {
      #    #CharCols <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
      #    CharCols <- c("Select All Samples",get_feat_cols(meta))
      #    selectInput("SubsetCol","Subset Samples By:",
      #                choices = CharCols,
      #                selected = SubCol_reactVal())
      #  }
      #  
      #})
      #observeEvent(input$DEGSubsetCol, {
      #  updateSelectInput(session, "SubsetCol",selected = isolate(input$DEGSubsetCol))
      #})
      #observeEvent(input$GSEASubsetCol, {
      #  updateSelectInput(session, "SubsetCol",selected = isolate(input$GSEASubsetCol))
      #})
      #observeEvent(input$SubsetCol,{SubCol_reactVal(input$SubsetCol)})
      
      
      output$rendSubsetCrit <- renderUI({
        
        req(input$SubsetCol)
        meta <- meta_input()
        if (input$SubsetCol != "Select All Samples") {
          SubCrit <- unique(meta[,input$SubsetCol])
          #selectInput("SubsetCrit","Sample Criteria:",choices = SubCrit, selected = SubCrit_reactVal())
          selectInput("SubsetCrit","Sample Criteria:",choices = SubCrit)
        }
        
      })
      #observeEvent(input$DEGSubsetCrit, {
      #  updateSelectInput(session, "SubsetCrit",selected = isolate(input$DEGSubsetCrit))
      #})
      #observeEvent(input$GSEASubsetCrit, {
      #  updateSelectInput(session, "SubsetCrit",selected = isolate(input$GSEASubsetCrit))
      #})
      observeEvent(input$SubsetCrit,{SubCrit_reactVal(input$SubsetCrit)})
      
      output$rendGroupCol <- renderUI({
        
        req(input$SubsetCol)
        meta <- meta_react()
        #FeatureChoices <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
        #FeatureChoices <- suppressWarnings(column_type_check(meta,"character"))[-1]
        FeatureChoices <- get_feat_cols(meta)
        if (input$SubsetCol != "Select All Samples") {
          FeatureChoices <- FeatureChoices[which(FeatureChoices!=input$SubsetCol)]
        }
        selectInput("GroupCol","Select Color Grouping Column:",
                    choices = FeatureChoices,
                    #selected = Feature_Selected,
                    selected = metacol_reactVal())
        
      })
      
      observeEvent(input$DEGGroupCol, {
        updateSelectInput(session, "GroupCol",selected = isolate(input$DEGGroupCol))
      })
      observeEvent(input$GSEAGroupCol, {
        updateSelectInput(session, "GroupCol",selected = isolate(input$GSEAGroupCol))
      })
      observeEvent(input$GroupColAvgHeat, {
        updateSelectInput(session, "GroupCol",selected = isolate(input$GroupColAvgHeat))
      })
      observeEvent(input$ssGSEAGroupColAvg, {
        updateSelectInput(session, "GroupCol",selected = isolate(input$ssGSEAGroupColAvg))
      })
      observeEvent(input$DEGcolHeat, {
        updateSelectInput(session, "GroupCol",selected = isolate(input$DEGcolHeat))
      })
      observeEvent(input$GroupCol,{metacol_reactVal(input$GroupCol)})
      output$rendGroupColMulti <- renderUI({
        
        req(input$SubsetCol)
        req(meta_react())
        meta <- meta_react()
        #FeatureChoices <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
        #FeatureChoices <- suppressWarnings(column_type_check(meta,"character"))[-1]
        FeatureChoices <- get_feat_cols(meta)
        if (input$SubsetCol != "Select All Samples") {
          FeatureChoices <- FeatureChoices[which(FeatureChoices!=input$SubsetCol)]
        }
        if (!isTruthy(Feature_Selected_react())) {
          if (!isTruthy(metacol_reactVal())) {
            initSelect <- FeatureChoices[1]
          } else {
            if (!metacol_reactVal() %in% FeatureChoices) {
              initSelect <- FeatureChoices[1]
            } else {
              initSelect <- metacol_reactVal()
            }
          }
        } else {
          initSelect <- metacol_reactVal()
        }
        selectInput("GroupColMulti","Select Annotation Data:",
                    choices = FeatureChoices,
                    selected = initSelect,
                    multiple = T)
        
      })
      output$rendGroupColAvgHeat <- renderUI({
        
        req(input$SubsetCol)
        meta <- meta_react()
        #FeatureChoices <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
        #FeatureChoices <- suppressWarnings(column_type_check(meta,"character"))[-1]
        FeatureChoices <- get_feat_cols(meta)
        if (input$SubsetCol != "Select All Samples") {
          FeatureChoices <- FeatureChoices[which(FeatureChoices!=input$SubsetCol)]
        }
        selectInput("GroupColAvgHeat","Select Annotation Data:",
                    choices = FeatureChoices,
                    #selected = Feature_Selected,
                    selected = metacol_reactVal())
        
      })
      observeEvent(input$GroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "GroupColAvgHeat",selected = isolate(input$GroupCol))
      })
      observeEvent(input$DEGGroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "GroupColAvgHeat",selected = isolate(input$DEGGroupCol))
      })
      observeEvent(input$GSEAGroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "GroupColAvgHeat",selected = isolate(input$GSEAGroupCol))
      })
      observeEvent(input$ssGSEAGroupColAvg, {
        updateSelectInput(session, "GroupColAvgHeat",selected = isolate(input$ssGSEAGroupColAvg))
      })
      observeEvent(input$DEGcolHeat, {
        updateSelectInput(session, "GroupColAvgHeat",selected = isolate(input$DEGcolHeat))
      })
      observeEvent(input$GroupColAvgHeat,{metacol_reactVal(input$GroupColAvgHeat)})
      
      output$rendDEGcolHeat <- renderUI({
        meta <- meta_react()
        selectInput("DEGcolHeat", "Differential Expression Feature:",
                    choices = colnames(meta)[-1], selected = metacol_reactVal())
      })
      observeEvent(input$GroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "DEGcolHeat",selected = isolate(input$GroupCol))
      })
      observeEvent(input$DEGGroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "DEGcolHeat",selected = isolate(input$DEGGroupCol))
      })
      observeEvent(input$GSEAGroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "DEGcolHeat",selected = isolate(input$GSEAGroupCol))
      })
      observeEvent(input$ssGSEAGroupColAvg, {
        updateSelectInput(session, "DEGcolHeat",selected = isolate(input$ssGSEAGroupColAvg))
      })
      observeEvent(input$GroupColAvgHeat, {
        updateSelectInput(session, "DEGcolHeat",selected = isolate(input$GroupColAvgHeat))
      })
      observeEvent(input$DEGcolHeat,{metacol_reactVal(input$DEGcolHeat)})
      
      
      ## DEG -----------------------------------------------------------------------
      #output$rendDEGSubsetCol <- renderUI({
      #  
      #  meta <- meta_react()
      #  if (ncol(meta) > 2) {
      #    #CharCols <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
      #    CharCols <- c("Select All Samples",get_feat_cols(meta))
      #    selectInput("DEGSubsetCol","Subset Samples By:",
      #                choices = CharCols,
      #                selected = SubCol_reactVal())
      #  }
      #  
      #})
      #observeEvent(input$SubsetCol, ignoreInit = TRUE, {
      #  updateSelectInput(session, "DEGSubsetCol",selected = isolate(input$SubsetCol))
      #})
      #observeEvent(input$GSEASubsetCol, ignoreInit = TRUE, {
      #  updateSelectInput(session, "DEGSubsetCol",selected = isolate(input$GSEASubsetCol))
      #})
      #observeEvent(input$DEGSubsetCol,{SubCol_reactVal(input$DEGSubsetCol)})
      
      #output$rendDEGSubsetCrit <- renderUI({
      #  
      #  req(input$DEGSubsetCol)
      #  if (input$DEGSubsetCol != "Select All Samples") {
      #    meta <- meta_input()
      #    SubCrit <- unique(meta[,input$DEGSubsetCol])
      #    selectInput("DEGSubsetCrit","Sample Criteria:",choices = SubCrit, selected = SubCrit_reactVal())
      #  }
      #  
      #})
      #observeEvent(input$SubsetCrit, ignoreInit = TRUE, {
      #  updateSelectInput(session, "DEGSubsetCrit",selected = isolate(input$SubsetCrit))
      #})
      #observeEvent(input$GSEASubsetCrit, ignoreInit = TRUE, {
      #  updateSelectInput(session, "DEGSubsetCrit",selected = isolate(input$GSEASubsetCrit))
      #})
      #observeEvent(input$DEGSubsetCrit,{SubCrit_reactVal(input$DEGSubsetCrit)})
      
      output$rendDEGGroupCol <- renderUI({
        
        req(input$SubsetCol)
        req(meta_react())
        meta <- meta_react()
        #FeatureChoices <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
        #FeatureChoices <- suppressWarnings(column_type_check(meta,"character"))[-1]
        FeatureChoices <- get_feat_cols(meta, keep_num = FALSE)
        if (input$SubsetCol != "Select All Samples") {
          FeatureChoices <- FeatureChoices[which(FeatureChoices!=input$SubsetCol)]
        }
        
        selectInput("DEGGroupCol","Select Differential Expression Feature:",
                    choices = FeatureChoices,
                    selected = metacol_reactVal())
        
      })
      observeEvent(input$GroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "DEGGroupCol",selected = isolate(input$GroupCol))
      })
      observeEvent(input$GSEAGroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "DEGGroupCol",selected = isolate(input$GSEAGroupCol))
      })
      observeEvent(input$GroupColAvgHeat, {
        updateSelectInput(session, "DEGGroupCol",selected = isolate(input$GroupColAvgHeat))
      })
      observeEvent(input$ssGSEAGroupColAvg, {
        updateSelectInput(session, "DEGGroupCol",selected = isolate(input$ssGSEAGroupColAvg))
      })
      observeEvent(input$DEGcolHeat, {
        updateSelectInput(session, "DEGGroupCol",selected = isolate(input$DEGcolHeat))
      })
      observeEvent(input$AvgGroupCol, {
        updateSelectInput(session, "DEGGroupCol",selected = isolate(input$AvgGroupCol))
      })
      observeEvent(input$DEGtabGroupCol, {
        updateSelectInput(session, "DEGGroupCol",selected = isolate(input$DEGtabGroupCol))
      })
      observeEvent(input$DEGGroupCol,{metacol_reactVal(input$DEGGroupCol)})
      
      observe({
        req(input$DEGGroupCol)
        req(meta_react())
        meta <- meta_react()
        #FeatureChoices <- suppressWarnings(column_type_check(meta,"character"))[-1]
        FeatureChoices <- get_feat_cols(meta, keep_num = FALSE)
        if (input$SubsetCol != "Select All Samples") {
          FeatureChoices <- FeatureChoices[which(FeatureChoices!=input$SubsetCol)]
        }
        FeatureChoices <- FeatureChoices[which(FeatureChoices != input$DEGGroupCol)]
        updateSelectizeInput(session,"DEGCovarSelect", choices = FeatureChoices, server = T)
      })
      
      output$rendAvgGroupCol <- renderUI({
        
        req(input$SubsetCol)
        meta <- meta_react()
        #FeatureChoices <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
        #FeatureChoices <- suppressWarnings(column_type_check(meta,"character"))[-1]
        FeatureChoices <- get_feat_cols(meta, keep_num = FALSE)
        if (input$SubsetCol != "Select All Samples") {
          FeatureChoices <- FeatureChoices[which(FeatureChoices!=input$SubsetCol)]
        }
        
        selectInput("AvgGroupCol","Select Comparison Feature:",
                    choices = FeatureChoices,
                    selected = metacol_reactVal())
        
      })
      observeEvent(input$GroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "AvgGroupCol",selected = isolate(input$GroupCol))
      })
      observeEvent(input$GSEAGroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "AvgGroupCol",selected = isolate(input$GSEAGroupCol))
      })
      observeEvent(input$GroupColAvgHeat, {
        updateSelectInput(session, "AvgGroupCol",selected = isolate(input$GroupColAvgHeat))
      })
      observeEvent(input$ssGSEAGroupColAvg, {
        updateSelectInput(session, "AvgGroupCol",selected = isolate(input$ssGSEAGroupColAvg))
      })
      observeEvent(input$DEGcolHeat, {
        updateSelectInput(session, "AvgGroupCol",selected = isolate(input$DEGcolHeat))
      })
      observeEvent(input$DEGGroupCol, {
        updateSelectInput(session, "AvgGroupCol",selected = isolate(input$DEGGroupCol))
      })
      observeEvent(input$DEGtabGroupCol, {
        updateSelectInput(session, "AvgGroupCol",selected = isolate(input$DEGtabGroupCol))
      })
      observeEvent(input$AvgGroupCol,{metacol_reactVal(input$AvgGroupCol)})
      
      output$rendDEGtabGroupCol <- renderUI({
        
        req(input$SubsetCol)
        meta <- meta_react()
        #FeatureChoices <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
        #FeatureChoices <- suppressWarnings(column_type_check(meta,"character"))[-1]
        FeatureChoices <- get_feat_cols(meta, keep_num = FALSE)
        if (input$SubsetCol != "Select All Samples") {
          FeatureChoices <- FeatureChoices[which(FeatureChoices!=input$SubsetCol)]
        }
        
        selectInput("DEGtabGroupCol","Select Differential Expression Feature:",
                    choices = FeatureChoices,
                    selected = metacol_reactVal())
        
      })
      observeEvent(input$GroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "DEGtabGroupCol",selected = isolate(input$GroupCol))
      })
      observeEvent(input$GSEAGroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "DEGtabGroupCol",selected = isolate(input$GSEAGroupCol))
      })
      observeEvent(input$GroupColAvgHeat, {
        updateSelectInput(session, "DEGtabGroupCol",selected = isolate(input$GroupColAvgHeat))
      })
      observeEvent(input$ssGSEAGroupColAvg, {
        updateSelectInput(session, "DEGtabGroupCol",selected = isolate(input$ssGSEAGroupColAvg))
      })
      observeEvent(input$DEGcolHeat, {
        updateSelectInput(session, "DEGtabGroupCol",selected = isolate(input$DEGcolHeat))
      })
      observeEvent(input$AvgGroupCol, {
        updateSelectInput(session, "DEGtabGroupCol",selected = isolate(input$AvgGroupCol))
      })
      observeEvent(input$DEGtabGroupCol,{metacol_reactVal(input$DEGtabGroupCol)})
      
      
      
      observe({
        req(input$DEGcolHeat)
        req(meta_react())
        meta <- meta_react()
        #FeatureChoices <- suppressWarnings(column_type_check(meta,"character"))[-1]
        FeatureChoices <- get_feat_cols(meta, keep_num = FALSE)
        if (isTruthy(input$SubsetCol)) {
          if (input$SubsetCol != "Select All Samples") {
            FeatureChoices <- FeatureChoices[which(FeatureChoices!=input$SubsetCol)]
          }
        }
        
        FeatureChoices <- FeatureChoices[which(FeatureChoices != input$DEGcolHeat)]
        updateSelectizeInput(session,"DEGCovarSelectHeat", choices = FeatureChoices, server = T)
      })
      
      observe({
        req(input$DEGGroupCol)
        req(meta_react())
        meta <- meta_react()
        #FeatureChoices <- suppressWarnings(column_type_check(meta,"character"))[-1]
        FeatureChoices <- get_feat_cols(meta, keep_num = FALSE)
        if (input$DEGGroupCol != "Select All Samples") {
          FeatureChoices <- FeatureChoices[which(FeatureChoices!=input$DEGGroupCol)]
        }
        FeatureChoices <- FeatureChoices[which(FeatureChoices != input$DEGGroupCol)]
        updateSelectizeInput(session,"DEGCovarSelectTab", choices = FeatureChoices, server = T)
      })
      
      
      
      ## GSEA ----------------------------------------------------------------------
      
      output$GroupSelectionError <- renderText({
        
        req(input$comparisonA)
        req(input$comparisonB)
        if (input$comparisonA == input$comparisonB) {
          paste("Comparison groups must be unique.")
        } else {
          req(meta_react())
          req(input$GSEAGroupCol)
          meta <- meta_react()
          groupCol <- input$GSEAGroupCol
          if (length(unique(meta[,groupCol])) == length(meta[,groupCol])) {
            paste("Must select a catagorical comparison group")
          } else {
            if (length(which(meta[,groupCol] == input$comparisonA)) <= 1 | length(which(meta[,groupCol] == input$comparisonB)) <= 1) {
              paste("Must select a comparison group with more than 1 observation")
            }
          }
          
        }
        
        
      })
      
      observe({
        req(meta_input())
        meta <- meta_input()
        if (ncol(meta) == 2) {
          metacol_reactVal(colnames(meta)[2])
        }
      })
      
      #output$rendGSEASubsetCol <- renderUI({
      #  
      #  meta <- meta_react()
      #  if (ncol(meta) > 2) {
      #    #if ()
      #    #CharCols <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
      #    CharCols <- c("Select All Samples",get_feat_cols(meta))
      #    selectInput("GSEASubsetCol","Subset Samples By:",
      #                choices = CharCols,
      #                selected = SubCol_reactVal())
      #  }
      #  
      #})
      #observeEvent(input$SubsetCol, ignoreInit = TRUE, {
      #  updateSelectInput(session, "GSEASubsetCol",selected = isolate(input$SubsetCol))
      #})
      #observeEvent(input$DEGSubsetCol, ignoreInit = TRUE, {
      #  updateSelectInput(session, "GSEASubsetCol",selected = isolate(input$DEGSubsetCol))
      #})
      #observeEvent(input$GSEASubsetCol,{SubCol_reactVal(input$GSEASubsetCol)})
      
      #output$rendGSEASubsetCrit <- renderUI({
      #  
      #  req(input$GSEASubsetCol)
      #  if (input$GSEASubsetCol != "Select All Samples") {
      #    meta <- meta_input()
      #    SubCrit <- unique(meta[,input$GSEASubsetCol])
      #    selectInput("GSEASubsetCrit","Sample Criteria:",choices = SubCrit, selected = SubCrit_reactVal())
      #  }
      #  
      #})
      #observeEvent(input$SubsetCrit, ignoreInit = TRUE, {
      #  updateSelectInput(session, "GSEASubsetCrit",selected = isolate(input$SubsetCrit))
      #})
      #observeEvent(input$DEGSubsetCrit, ignoreInit = TRUE, {
      #  updateSelectInput(session, "GSEASubsetCrit",selected = isolate(input$DEGSubsetCrit))
      #})
      #observeEvent(input$GSEASubsetCrit,{SubCrit_reactVal(input$GSEASubsetCrit)})
      
      output$rendGSEAGroupCol <- renderUI({
        
        req(input$SubsetCol)
        meta <- meta_react()
        #FeatureChoices <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
        #FeatureChoices <- suppressWarnings(column_type_check(meta,"character"))[-1]
        FeatureChoices <- get_feat_cols(meta, keep_num = FALSE)
        if (input$SubsetCol != "Select All Samples") {
          FeatureChoices <- FeatureChoices[which(FeatureChoices!=input$SubsetCol)]
        }
        selectInput("GSEAGroupCol","Select GSEA Grouping Column:",
                    choices = FeatureChoices,
                    selected = metacol_reactVal())
        
      })
      observeEvent(input$GroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "GSEAGroupCol",selected = isolate(input$GroupCol))
      })
      observeEvent(input$DEGGroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "GSEAGroupCol",selected = isolate(input$DEGGroupCol))
      })
      observeEvent(input$GroupColAvgHeat, {
        updateSelectInput(session, "GSEAGroupCol",selected = isolate(input$GroupColAvgHeat))
      })
      observeEvent(input$ssGSEAGroupColAvg, {
        updateSelectInput(session, "GSEAGroupCol",selected = isolate(input$ssGSEAGroupColAvg))
      })
      observeEvent(input$DEGcolHeat, {
        updateSelectInput(session, "GSEAGroupCol",selected = isolate(input$DEGcolHeat))
      })
      observeEvent(input$GSEAGroupCol,{metacol_reactVal(input$GSEAGroupCol)})
      output$rendssGSEAGroupColAvg <- renderUI({
        
        req(input$SubsetCol)
        meta <- meta_react()
        #FeatureChoices <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
        #FeatureChoices <- suppressWarnings(column_type_check(meta,"character"))[-1]
        FeatureChoices <- get_feat_cols(meta, keep_num = FALSE)
        if (input$SubsetCol != "Select All Samples") {
          FeatureChoices <- FeatureChoices[which(FeatureChoices!=input$SubsetCol)]
        }
        selectInput("ssGSEAGroupColAvg","Select Grouping Data:",
                    choices = FeatureChoices,
                    selected = metacol_reactVal())
        
      })
      observeEvent(input$GroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "ssGSEAGroupColAvg",selected = isolate(input$GroupCol))
      })
      observeEvent(input$DEGGroupCol, ignoreInit = TRUE, {
        updateSelectInput(session, "ssGSEAGroupColAvg",selected = isolate(input$DEGGroupCol))
      })
      observeEvent(input$GroupColAvgHeat, {
        updateSelectInput(session, "ssGSEAGroupColAvg",selected = isolate(input$GroupColAvgHeat))
      })
      observeEvent(input$ssGSEAGroupColAvg, {
        updateSelectInput(session, "ssGSEAGroupColAvg",selected = isolate(input$ssGSEAGroupColAvg))
      })
      observeEvent(input$DEGcolHeat, {
        updateSelectInput(session, "ssGSEAGroupColAvg",selected = isolate(input$DEGcolHeat))
      })
      observeEvent(input$ssGSEAGroupColAvg,{metacol_reactVal(input$ssGSEAGroupColAvg)})
      output$rendssGSEAGroupColMulti <- renderUI({
        
        req(input$SubsetCol)
        meta <- meta_react()
        #FeatureChoices <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
        #FeatureChoices <- suppressWarnings(column_type_check(meta,"character"))[-1]
        FeatureChoices <- get_feat_cols(meta)
        if (input$SubsetCol != "Select All Samples") {
          FeatureChoices <- FeatureChoices[which(FeatureChoices!=input$SubsetCol)]
        }
        selectInput("ssGSEAGroupColMulti","Select Annotation Data:",
                    choices = FeatureChoices,
                    selected = metacol_reactVal(),
                    multiple = T)
        
      })
      output$rendGSEAGroupColMulti <- renderUI({
        
        req(input$SubsetCol)
        meta <- meta_react()
        #FeatureChoices <- c("Select All Samples",suppressWarnings(column_type_check(meta,"character"))[-1])
        #FeatureChoices <- suppressWarnings(column_type_check(meta,"character"))[-1]
        FeatureChoices <- get_feat_cols(meta)
        if (input$SubsetCol != "Select All Samples") {
          FeatureChoices <- FeatureChoices[which(FeatureChoices!=input$SubsetCol)]
        }
        selectInput("GSEAGroupColMulti","Select Annotation Data:",
                    choices = FeatureChoices,
                    selected = metacol_reactVal(),
                    multiple = T)
        
      })
      
      
      meta_react <- reactive({
        
        req(meta_input())
        meta <- meta_input()
        
        if (ncol(meta) > 2) {
          req(input$SubsetCol)
          SubCol <- input$SubsetCol
          SubCrit <- input$SubsetCrit
          if (SubCol != "Select All Samples") {
            if (SubCol %in% colnames(meta)){
              if (isTruthy(SubCrit)) {
                meta <- meta[which(meta[,SubCol] == SubCrit),]
                meta
              }
            }
          } #else {
            #meta
          #}
        } #else {
          #meta
        #}
        
        if (isTruthy(ImmDeconv_react())) {
          imm_deconv <- ImmDeconv_react()
          imm_deconv <- as.data.frame(t(imm_deconv))
          meta <- merge(meta,imm_deconv, by.x = colnames(meta)[1], by.y = 0, all.x = T)
        }
        
        meta
      })
      
      #observe({
      #  
      #  req(meta_react())
      #  meta <- meta_react()
      #  updateSelectizeInput(session = session, inputId = "userheatsamp2",
      #                       choices = meta[,1],selected = meta[,1], server = T)
      #  updateSelectizeInput(session = session, inputId = "userheatsampSS",
      #                       choices = meta[,1],selected = meta[,1], server = T)
      #  
      #})
      
      observe({
        
        req(meta_react())
        req(input$GroupColMulti)
        meta <- meta_react()
        metacol <- input$GroupColMulti[1]
        meta2 <- meta[,c(colnames(meta)[1],metacol)]
        meta2[,2] <- as.factor(meta2[,2])
        meta2 <- meta2[order(meta2[,2]),]
        
        updateSelectizeInput(session,"userheatsamp2",choices = meta2[,1], selected = meta2[,1])
      })
      
      
      expr_react <- reactive({
        req(meta_react())
        meta <- meta_react()
        if (immudecon_check & input$ImmDeconvAnalysis == "Immune Deconvolution") {
          req(ImmDeconv_react())
          expr <- ImmDeconv_react()
        } else {
          req(expr_input())
          expr <- expr_input()
        }
        expr2 <- expr[,meta[,1]]
        expr2
        
      })
      observe({
        if (immudecon_check & input$ImmDeconvAnalysis == "Immune Deconvolution") {
          req(ImmDeconv_react())
          req(meta_react())
          meta <- meta_react()
          req(expr_input())
          expr <- expr_input()
          expr <- expr[,meta[,1]]
        } else {
          req(expr_react())
          expr <- expr_react()
        }
        A <- as.matrix(expr)
        A_raw(A)
      })
      
      #expr_mat_react <- reactive({
      #
      #  meta <- meta_react()
      #  A <- A_raw()
      #  A2 <- A[,meta[,1]]
      #  A2
      #
      #})
      
      #observe({
      #  Gene <- Gene_raw()
      #  updateSelectizeInput(session = session, inputId = "heatmapGeneSelec",
      #                       choices = Gene,selected = CTKgenes, server = T)
      #})
      #observeEvent(input$heatmapGeneSelec,{
      #  heatmapGeneSelec_react <- reactive({
      #    input$heatmapGeneSelec
      #  })
      #})
      #observe({
      #  meta <- meta_react()
      #  expr <- expr_react()
      #  sampsames <- intersect(colnames(expr),meta[,1])
      #  updateSelectizeInput(session = session, inputId = "userheatsamp2",
      #                       choices = sampsames,selected = sampsames, server = T)
      #})
      #observe({
      #  Gene <- Gene_raw()
      #  updateSelectizeInput(session = session, inputId = "avgheatmapGeneSelec",
      #                       choices = Gene,selected = CTKgenes, server = T)
      #})
      #observe({
      #  Gene <- Gene_raw()
      #  updateSelectizeInput(session = session, inputId = "scatterG1",
      #                       choices = Gene,selected = Gene[1], server = T)
      #})
      #observe({
      #  Gene <- Gene_raw()
      #  updateSelectizeInput(session = session, inputId = "scatterG2",
      #                       choices = Gene,selected = Gene[2], server = T)
      #})
      #observe({
      #  geneList <- geneList_raw()
      #  updateSelectizeInput(session = session, inputId = "userGeneSelec",
      #                       choices = sort(as.vector(geneList[,1])),selected = NULL, server = T)
      #})
      #observe({
      #  updateSelectizeInput(session = session, inputId = "scatterGeneSelec",
      #                       choices = rownames(expr),selected = NULL, server = T)
      #})
      #observe({
      #  meta <- meta_react()
      #  expr <- expr_react()
      #  sampsames <- intersect(colnames(expr),meta[,1])
      #  updateSelectizeInput(session = session, inputId = "userheatsampSS",
      #                       choices = sampsames,selected = sampsames, server = T)
      #})
      
      #output$rendssgseaHeatGS <- renderUI({
      #  
      #  GS <- geneset_list()
      #  gs_names <- names(GS)
      #  
      #  #if (input$tables == 1) {
      #  #  gs_names <- c(grep("^HALLMARK",names(gs()),value = T),grep("^HALLMARK",names(gs()),value = T, invert = T))
      #  #}
      #  #if (input$tables == 3) {
      #  #  gs_names <- names(gs2())
      #  #}
      #  #if (input$tables == 5) {
      #  #  gmt <- GStable.ubg()
      #  #  colnames(gmt) <- c("term","gene")
      #  #  #gs_name <- user_gs_mirror()[input$GStable.u_rows_selected,1]
      #  #  #gmt_sub <- gmt[which(gmt$term == gs_name),]
      #  #  GS <- list()
      #  #  for (i in unique(gmt[,1])){
      #  #    GS[[i]] <- gmt[gmt[,1] == i,]$gene
      #  #  }
      #  #  #GS <- RDataListGen()[geneset_names]()[(user_gs_mirror()[input$GStable.u_rows_selected,1])]
      #  #  gs_names <- names(GS)
      #  #}
      #  # take NA out if number of gs is less than 50
      #  gs_selected <- gs_names[1:3]
      #  gs_selected <- gs_selected[!is.na(gs_selected)]
      #  selectInput("ssgseaHeatGS", "Select Gene Sets:",
      #              choices = gs_names, selected = gs_selected, multiple = T)
      #  
      #})
      
      # Render cluster method selection for average heatmap
      output$ClusterMethodMVG <- renderUI({
        
        if (input$clustrowsMVG == TRUE | input$clustcolsMVG == TRUE) {
          
          fluidRow(
            column(6,
                   selectInput("ClusteringMethod",
                               "Select Clustering Method",
                               choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
            ),
            column(6,
                   numericInput("NumClusters", step = 1, label = "Num of Row Clusters (k)", value = NA)
            )
          )
          
          
          
        }
        
      })
      
      #output$rendMVGheatAnnoCol <- renderUI({
      #
      #  if (ncol(meta) > 2) {
      #    selectInput("MVGheatAnnoCol","Select Annotation Column:",
      #                choices = colnames(meta[,2:ncol(meta)]),
      #                selected = Feature_Selected)
      #  }
      #
      #})
      #
      #output$rendCustomheatAnnoCol <- renderUI({
      #
      #  if (ncol(meta) > 2) {
      #    selectInput("CustomheatAnnoCol","Select Annotation Column:",
      #                choices = colnames(meta[,2:ncol(meta)]),
      #                selected = Feature_Selected)
      #  }
      #
      #})
      
      #output$rendMetaColDEGHeat <- renderUI({
      #
      #  if (ncol(meta) > 2){
      #    selectInput("MetaColDEGHeat","Meta Column:",
      #                choices = colnames(meta)[2:ncol(meta)],
      #                selected = Feature_Selected)
      #  }
      #
      #})
      #
      #output$rendBoxPlot1MetaCol <- renderUI({
      #
      #  if (ncol(meta) > 2) {
      #    selectInput("BoxPlot1MetaCol","Select Annotation Column:",
      #                choices = colnames(meta[,2:ncol(meta)]),
      #                selected = Feature_Selected)
      #  }
      #
      #})
      #
      #output$rendBarPlot1MetaCol <- renderUI({
      #
      #  if (ncol(meta) > 2) {
      #    selectInput("BarPlot1MetaCol","Select Annotation Column:",
      #                choices = colnames(meta[,2:ncol(meta)]),
      #                selected = Feature_Selected)
      #  }
      #
      #})
      #
      #output$rendAvgHeatMetaCol <- renderUI({
      #
      #  if (ncol(meta) > 2) {
      #    selectInput("AvgHeatMetaCol","Select Meta Column:",
      #                choices = colnames(meta[,2:ncol(meta)]),
      #                selected = Feature_Selected)
      #  }
      #
      #})
      
      output$rendSampCondSelectionCust <- renderUI({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        metagroups <- unique(meta[,metacol])
        selectInput("SampCondSelectionCust", "Sample Conditon Selection",
                    choices = metagroups, selected = metagroups, multiple = T)
        #if (ncol(meta) > 2) {
        #  metagroups_new <- as.vector(levels(factor(meta[,input$AvgHeatMetaCol])))
        #  selectInput("SampCondSelectionCust","Sample Coniditon Selection",
        #              choices = metagroups_new, selected = metagroups_new[c(1:length(metagroups_new))],
        #              multiple = T)
        #}
        #else if (ncol(meta) == 2) {
        #  selectInput("SampCondSelectionCust","Sample Coniditon Selection",
        #              choices = metagroups, selected = metagroups[c(1:length(metagroups))],
        #              multiple = T)
        #}
        
      })
      
      output$rendDEGcolHeat <- renderUI({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        selectInput("DEGcolHeat", "Differential Expression Feature:",
                    choices = colnames(meta)[-1], selected = metacol[1])
      })
      
      output$rendcomparisonA2_h <- renderUI({
        
        req(input$DEGcolHeat)
        meta <- meta_react()
        #metacol <- metacol_reactVal()
        metacol <- input$DEGcolHeat
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonA2_h", "Comparison: GroupA",
                    choices = metagroups, selected = metagroups[1])
      })
      
      output$rendcomparisonA2_h_one <- renderUI({
        
        req(input$DEGcolHeat)
        meta <- meta_react()
        #metacol <- metacol_reactVal()
        metacol <- input$DEGcolHeat
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonA2_h_one", "Comparison: GroupA",
                    choices = metagroups, selected = metagroups[1])
      })
      
      #output$rendGSEAmetaCol <- renderUI({
      #
      #  if (ncol(meta) > 2) {
      #    selectInput("GSEAmetaCol","Select Meta Column:",
      #                choices = colnames(meta[,2:ncol(meta)]),
      #                selected = Feature_Selected)
      #  }
      #
      #})
      
      output$rendcomparisonA <- renderUI({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonA", "Comparison: GroupA",
                    choices = metagroups, selected = metagroups[1])
      })
      
      output$rendcomparisonA_one <- renderUI({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonA_one", "Comparison Group",
                    choices = metagroups, selected = metagroups[1])
      })
      
      output$rendcomparisonB <- renderUI({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonB", "Comparison: GroupB",
                    choices = metagroups, selected = metagroups[2])
        #if (ncol(meta) == 2){
        #  selectInput("comparisonB", "Comparison: GroupB",
        #              choices = metagroups, selected = metagroups[2])
        #}
        #else if (ncol(meta) > 2){
        #  metagroups_new <- as.vector(levels(factor(meta[,input$GSEAmetaCol])))
        #  selectInput("comparisonB", "Comparison: GroupB",
        #              choices = metagroups_new, selected = metagroups_new[2])
        #}
        
      })
      
      output$rendcomparisonB2_h <- renderUI({
        
        req(input$DEGcolHeat)
        meta <- meta_react()
        #metacol <- metacol_reactVal()
        metacol <- input$DEGcolHeat
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonB2_h", "Comparison: GroupB",
                    choices = metagroups, selected = metagroups[2])
        #if (ncol(meta) == 2){
        #  selectInput("comparisonB2_h", "Comparison: GroupB",
        #              choices = metagroups, selected = metagroups[2])
        #}
        #else if (ncol(meta) > 2){
        #  metagroups_new <- as.vector(levels(factor(meta[,input$MetaColDEGHeat])))
        #  selectInput("comparisonB2_h", "Comparison: GroupB",
        #              choices = metagroups_new, selected = metagroups_new[2])
        #}
        
      })
      
      #output$rendDEGMetaCol <- renderUI({
      #
      #  if (ncol(meta) > 2) {
      #    selectInput("DEGMetaCol","Select Meta Column:",
      #                choices = colnames(meta[,2:ncol(meta)]),
      #                selected = Feature_Selected)
      #  }
      #
      #})
      
      output$rendcomparisonA2 <- renderUI({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonA2", "Comparison: GroupA",
                    choices = metagroups, selected = metagroups[1])
        
      })
      
      output$rendcomparisonA2_one <- renderUI({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonA2_one", "Comparison Group",
                    choices = metagroups, selected = metagroups[1])
        
      })
      
      output$rendcomparisonB2 <- renderUI({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonB2", "Comparison: GroupB",
                    choices = metagroups, selected = metagroups[2])
        #if (ncol(meta) == 2){
        #  selectInput("comparisonB2", "Comparison: GroupB",
        #              choices = metagroups, selected = metagroups[2])
        #}
        #else if (ncol(meta) > 2){
        #  metagroups_new <- as.vector(levels(factor(meta[,input$DEGMetaCol])))
        #  selectInput("comparisonB2", "Comparison: GroupB",
        #              choices = metagroups_new, selected = metagroups_new[2])
        #}
        
      })
      
      #output$rendAvgExprScatterMetaCol <- renderUI({
      #
      #  if (ncol(meta) > 2) {
      #    selectInput("AvgExprScatterMetaCol","Select Meta Column:",
      #                choices = colnames(meta[,2:ncol(meta)]),
      #                selected = Feature_Selected)
      #  }
      #
      #})
      
      output$rendcomparisonA2.avg <- renderUI({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonA2.avg", "Comparison: GroupA",
                    choices = metagroups, selected = metagroups[1])
        #if (ncol(meta) == 2){
        #  selectInput("comparisonA2.avg", "Comparison: GroupA",
        #              choices = metagroups, selected = metagroups[1])
        #}
        #else if (ncol(meta) > 2){
        #  metagroups_new <- as.vector(levels(factor(meta[,input$AvgExprScatterMetaCol])))
        #  selectInput("comparisonA2.avg", "Comparison: GroupA",
        #              choices = metagroups_new, selected = metagroups_new[1])
        #}
        
      })
      
      output$rendcomparisonB2.avg <- renderUI({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonB2.avg", "Comparison: GroupB",
                    choices = metagroups, selected = metagroups[2])
        #if (ncol(meta) == 2){
        #  selectInput("comparisonB2.avg", "Comparison: GroupB",
        #              choices = metagroups, selected = metagroups[2])
        #}
        #else if (ncol(meta) > 2){
        #  metagroups_new <- as.vector(levels(factor(meta[,input$AvgExprScatterMetaCol])))
        #  selectInput("comparisonB2.avg", "Comparison: GroupB",
        #              choices = metagroups_new, selected = metagroups_new[2])
        #}
        
      })
      
      #output$rendScatterPlotMetaCol <- renderUI({
      #
      #  if (ncol(meta) > 2){
      #    selectInput("ScatterPlotMetaCol","Select Sample Condition",
      #                choices = colnames(meta)[2:ncol(meta)],
      #                selected = Feature_Selected)
      #  }
      #
      #})
      #
      #output$rendEnrichRPathMetaCol <- renderUI({
      #
      #  if (ncol(meta) > 2){
      #    selectInput("EnrichRPathMetaCol","Select Sample Condition",
      #                choices = colnames(meta)[2:ncol(meta)],
      #                selected = Feature_Selected)
      #  }
      #
      #})
      
      #output$rendcomparisonA2.path <- renderUI({
      #
      #  if (ncol(meta) == 2){
      #    selectInput("comparisonA2.path", "Comparison: GroupA",
      #                choices = metagroups, selected = metagroups[1])
      #  }
      #  else if (ncol(meta) > 2){
      #    metagroups_new <- as.vector(levels(factor(meta[,input$EnrichRPathMetaCol])))
      #    selectInput("comparisonA2.path", "Comparison: GroupA",
      #                choices = metagroups_new, selected = metagroups_new[1])
      #  }
      #
      #})
      #
      #output$rendcomparisonB2.path <- renderUI({
      #
      #  if (ncol(meta) == 2){
      #    selectInput("comparisonB2.path", "Comparison: GroupB",
      #                choices = metagroups, selected = metagroups[2])
      #  }
      #  else if (ncol(meta) > 2){
      #    metagroups_new <- as.vector(levels(factor(meta[,input$EnrichRPathMetaCol])))
      #    selectInput("comparisonB2.path", "Comparison: GroupB",
      #                choices = metagroups_new, selected = metagroups_new[2])
      #  }
      #
      #})
      
      #output$rendDEGtableMetaCol <- renderUI({
      #
      #  if (ncol(meta) > 2){
      #    selectInput("DEGtableMetaCol","Select Sample Condition",
      #                choices = colnames(meta)[2:ncol(meta)],
      #                selected = Feature_Selected)
      #  }
      #
      #})
      
      output$rendcomparisonA2.DEG <- renderUI({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonA2.DEG", "Comparison: GroupA",
                    choices = metagroups, selected = metagroups[1])
      })
      
      output$rendcomparisonA2.DEG_one <- renderUI({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonA2.DEG_one", "Comparison Group",
                    choices = metagroups, selected = metagroups[1])
      })
      
      output$rendcomparisonB2.DEG <- renderUI({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        metagroups <- unique(meta[,metacol])
        metagroups <- sort(metagroups,na.last = TRUE)
        selectInput("comparisonB2.DEG", "Comparison: GroupB",
                    choices = metagroups, selected = metagroups[2])
        #if (ncol(meta) == 2){
        #  selectInput("comparisonB2.DEG", "Comparison: GroupB",
        #              choices = metagroups, selected = metagroups[2])
        #}
        #else if (ncol(meta) > 2){
        #  metagroups_new <- as.vector(levels(factor(meta[,input$DEGtableMetaCol])))
        #  selectInput("comparisonB2.DEG", "Comparison: GroupB",
        #              choices = metagroups_new, selected = metagroups_new[2])
        #}
        
      })
      
      #output$rendBoxPlotMetaColSelec <- renderUI({
      #
      #  if (ncol(meta) > 2){
      #    selectInput("BoxPlotMetaColSelec","Select Sample Condition",
      #                choices = colnames(meta)[2:ncol(meta)],
      #                selected = Feature_Selected)
      #  }
      #
      #})
      
      # Render cluster method selection for average heatmap
      output$rendClustMethodsCust <- renderUI({
        
        if (input$ClusRowOptCust == TRUE | input$ClusColOptCust == TRUE) {
          
          selectInput("ClusteringMethodCust",
                      "Select Clustering Method",
                      choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
          
        }
        
      })
      
      # Render cluster method selection for average heatmap
      output$rendClustMethodsDEG <- renderUI({
        
        if (input$ClusRowOptDEG == TRUE | input$ClusColOptDEG == TRUE) {
          
          selectInput("ClusteringMethod_degh",
                      "Select Clustering Method",
                      choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
          
        }
        
      })
      
      output$rendAvgClustMethodsCust <- renderUI({
        
        if (input$AvgClusRowOptCust == TRUE | input$AvgClusColOptCust == TRUE) {
          
          selectInput("AvgClusteringMethodCust",
                      "Select Clustering Method",
                      choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
          
        }
        
      })
      
      # Render cluster method selection for average heatmap
      output$rendClusterMethodSSheat <- renderUI({
        
        if (input$clustcolsSSheat == TRUE | input$clustrowsSSheat == TRUE) {
          
          selectInput("ClusterMethodSSheat",
                      "Select Clustering Method",
                      choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
          
        }
        
      })
      
      
      
      #render enrich signature table generation button for msigdb data
      #output$genESTbutton <- renderUI({
      #  if (input$tables == 1) {
      #    
      #  }
      #})
      
      #render user gmt data upload if indicated
      output$user.gmt <- renderUI({
        fileInput("user.gmt.file", "Gene Set File (.gmt, .tsv, or .txt)", accept = c(".gmt",".tsv",".txt"))
      })
      
      #render gene set table based off gmt file given
      output$user.GStable <- renderUI({
        req(input$user.gmt.file)
        div(DT::dataTableOutput("GStable.u"), style = "font-size:10px; height:500px; overflow-X: scroll")
      })
      
      #Warning message for RData list generation
      output$RDataMessage <- renderUI({
        req(input$user.gmt.file)
        p("Generating RData list may take up to several minutes depending on size of GMT file.")
      })
      
      #render action button to create RData list for ssGSEA boxplots
      #output$user.RDataButton <- renderUI({
      #  req(input$user.gmt.file)
      #  actionButton("user.RData.Gen", "Generate RData list for ssGSEA Boxplot and Heatmap")
      #})
      
      observe({
        req(input$PorAdjPval)
        if (input$PorAdjPval == "-log10(P.value)") {
          updateNumericInput(session,"p_cutoff", label = "P.Value Threshold:",
                             min = 0, max = 10, step = 0.1, value = 0.05)
        } else if (input$PorAdjPval == "-log10(Adjusted P.value)") {
          updateNumericInput(session,"p_cutoff", label = "Adj.P.Value Threshold:",
                             min = 0, max = 10, step = 0.1, value = 1)
        }
      })
      
      #render UI for hover text in volcano plot
      output$hover_info <- renderUI({
        req(topgenereact())
        top2 <- topgenereact()
        df <- top2 %>%
          select(GeneName,logFC,P.Value,adj.P.Val)
        colnames(df)[3] <- "-log10(P.Value)"
        df$P.Value <- df$'-log10(P.Value)'
        df$`-log10(P.Value)` <- -log10(df$`-log10(P.Value)`)
        
        colnames(df)[4] <- "-log10(adj.P.Val)"
        df$adj.P.Val <- df$'-log10(adj.P.Val)'
        df$`-log10(adj.P.Val)` <- -log10(df$`-log10(adj.P.Val)`)
        hover <- input$plot_hover
        point <- nearPoints(df, hover, threshold = 10, maxpoints = 1, addDist = FALSE)
        if (nrow(point) == 0) return(NULL)
        wellPanel(
          p(HTML(paste0("<b> Name: </b>", point[1], "<br/>",
                        "<b> Fold change: </b>", round(point[2], digits = 4), "<br/>",
                        "<b> P Value: </b>", point[5], "<br/>",
                        "<b> adj P Value: </b>", point[4], "<br/>",
                        NULL
          ))))
        
        
      })
      
      #render UI for hover text in MA plot
      output$hover_info2 <- renderUI({
        req(topgenereact())
        top2 <- topgenereact()
        if (input$volcanoCompChoice == "DESeq2") {
          top2$AveExpr <- log2(top2$AveExpr+1)
        }
        df <- top2 %>%
          select(GeneName,AveExpr,logFC,P.Value,adj.P.Val)
        colnames(df)[4] <- "-log10(P.Value)"
        df$P.Value <- df$'-log10(P.Value)'
        df$`-log10(P.Value)` <- -log10(df$`-log10(P.Value)`)
        hover <- input$plot_hover2
        point <- nearPoints(df, hover, threshold = 10, maxpoints = 1, addDist = FALSE)
        if (nrow(point) == 0) return(NULL)
        wellPanel(
          p(HTML(paste0("<b> Name: </b>", point[1], "<br/>",
                        "<b> Fold change: </b>", round(point[3], digits = 4), "<br/>",
                        "<b> P Value: </b>", point[6], "<br/>",
                        "<b> adj P Value: </b>", point[5], "<br/>",
                        "<b> Avg. Expression: </b>", round(point[2], digits = 4), "<br/>",
                        NULL
          ))))
      })
      
      AvgExprReact <- reactive({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        expr <- expr_react()
        A_choice <- input$comparisonA2.avg
        B_choice <- input$comparisonB2.avg
        log_choice <- input$AvgExpLogFC
        #if (ncol(meta) > 2) {
        #  metacol <- input$AvgExprScatterMetaCol
        #}
        #else if (ncol(meta) == 2) {
        #  metacol <- colnames(meta)[2]
        #}
        # Get A and B sample names
        A <- meta[which(meta[,metacol] == A_choice),1]
        B <- meta[which(meta[,metacol] == B_choice),1]
        # Make A and B expression Matrices
        mat_A <- expr[,A]
        mat_B <- expr[,B]
        logFC_text <- ""
        # Log if user designates
        if (log_choice == TRUE) {
          mat_A <- log2(mat_A + 1)
          mat_B <- log2(mat_B + 1)
          logFC_text <- " (log2 +1)"
        }
        mat_A <- as.data.frame(mat_A)
        mat_B <- as.data.frame(mat_B)
        # Get avg expression of each gene
        mat_A$AvgExpression_GroupA <- rowMeans(mat_A)
        mat_B$AvgExpression_GroupB <- rowMeans(mat_B)
        # Add Gene Names as column
        mat_A$GeneSymbol <- rownames(mat_A)
        mat_B$GeneSymbol <- rownames(mat_B)
        # Merge average columns
        AvgExpr_Table <- merge(mat_A[,c("AvgExpression_GroupA","GeneSymbol")],
                               mat_B[,c("AvgExpression_GroupB","GeneSymbol")],
                               by = "GeneSymbol",
                               all = T)
        AvgExpr_Table
        
      })
      
      #render UI for hover text in volcano plot
      output$hover_info3 <- renderUI({
        req(AvgExprReact())
        
        #top2 <- topgenereact()
        #df2 <- top2 %>%
        #  select(GeneName,AveExpr,logFC,P.Value,adj.P.Val)
        #colnames(df2)[4] <- "-log10(P.Value)"
        #df2$P.Value <- df2$'-log10(P.Value)'
        #df2$`-log10(P.Value)` <- -log10(df2$`-log10(P.Value)`)
        
        A_choice <- input$comparisonA2.avg
        B_choice <- input$comparisonB2.avg
        df <- AvgExprReact()
        hover <- input$plot_hover3
        point <- nearPoints(df, hover, threshold = 10, maxpoints = 1, addDist = FALSE)
        if (nrow(point) == 0) return(NULL)
        wellPanel(
          p(HTML(paste0("<b> Gene: </b>", point[1], "<br/>",
                        "<b> ",A_choice," Average Expression: </b>", round(point[2], digits = 4), "<br/>",
                        "<b> ",B_choice," Average Expression: </b>", round(point[3], digits = 4), "<br/>",
                        NULL
          ))))
      })
      
      ####----Reactives----####
      
      
      ssGSEA_Heat_GS <- reactive({
        
        gs <- geneset_list()
        gs_names_start <- names(gs)
        gs_names_start <- gs_names_start[1:3]
        #if (input$tables == 1) {
        #  ## starting ssgsea heatmap genes
        #  gs_names_start <- c(grep("^HALLMARK",names(gs()),value = T),grep("^HALLMARK",names(gs()),value = T, invert = T))
        #  gs_names_start <- gs_names_start[1:3]
        #}
        #else if (input$tables == 3) {
        #  gs_names_start <- names(gs2())[1:3]
        #}
        #else if (input$tables == 5) {
        #  gmt <- GStable.ubg()
        #  colnames(gmt) <- c("term","gene")
        #  #gs_name <- user_gs_mirror()[input$GStable.u_rows_selected,1]
        #  #gmt_sub <- gmt[which(gmt$term == gs_name),]
        #  GS <- list()
        #  for (i in unique(gmt[,1])){
        #    GS[[i]] <- gmt[gmt[,1] == i,]$gene
        #  }
        #  #GS <- RDataListGen()[geneset_names]()[(user_gs_mirror()[input$GStable.u_rows_selected,1])]
        #  gs_names_start <- names(GS)[1:3]
        #}
        gs_names_start
        
      })
      
      GeneratedMSigDBEST <- eventReactive(input$GenerateEST,{
        #if (input$GenerateEST == TRUE) {
        #req(input$comparisonA)
        #req(input$comparisonB)
        req(geneset_gmt())
        req(meta_react())
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        
        if (input$volcanoCompChoice3 == "Two groups") {
          Acomp <- input$comparisonA
          Bcomp <- input$comparisonB
          groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
          groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
        } else if (input$volcanoCompChoice3 == "One group") {
          Acomp <- input$comparisonA_one
          Bcomp <- paste0("Not_",Acomp)
          meta[,paste0(metacol,"_Dichot")] <- NA
          meta[which(meta[,metacol] == Acomp),paste0(metacol,"_Dichot")] <- Acomp
          meta[which(meta[,metacol] != Acomp),paste0(metacol,"_Dichot")] <- paste0("Not_",Acomp)
          metacol <- paste0(metacol,"_Dichot")
          groupA <- meta[which(meta[,metacol] == Acomp),1]
          groupB <- meta[which(meta[,metacol] == Bcomp),1]
        }
        
        
        
        
        #meta <- meta_react()
        #metacol <- metacol_reactVal()
        #if (length(colnames(meta)) == 2) {
        #  metacol <- colnames(meta)[2]
        #}
        A <- A_raw()
        gmt <- geneset_gmt()
        #withProgress(message = "Processing", value = 0, {
        #  incProgress(0.25, detail = "Calculating Signal-to-Noise")
        ##----Signal-to-Noise Calculation----##
        A <- A + 0.00000001
        P = as.matrix(as.numeric(colnames(A) %in% groupA))
        n1 <- sum(P[,1])
        M1 <- A %*% P
        M1 <- M1/n1
        A2 <- A*A
        S1 <- A2 %*% P
        S1 <- S1/n1 - M1*M1
        S1 <- sqrt(abs((n1/(n1-1)) * S1))
        P = as.matrix(as.numeric(colnames(A) %in% groupB))
        n2 <- sum(P[,1])
        M2 <- A %*% P
        M2 <- M2/n2
        A2 <- A*A
        S2 <- A2 %*% P
        S2 <- S2/n2 - M2*M2
        S2 <- sqrt(abs((n2/(n2-1)) * S2))
        rm(A2)
        # small sigma "fix" as used in GeneCluster
        S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
        S2 <- ifelse(S2 == 0, 0.2, S2)
        S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
        S1 <- ifelse(S1 == 0, 0.2, S1)
        M1 <- M1 - M2
        rm(M2)
        S1 <- S1 + S2
        rm(S2)
        s2n.matrix <- M1/S1
        ##----Reformatting----##
        #incProgress(0.25, detail = "Formatting Data")
        s2n.df <- as.data.frame(s2n.matrix)
        s2n.df$GeneID <- rownames(s2n.df)
        rownames(s2n.df) <- NULL
        data <- dplyr::select(s2n.df, GeneID, V1)
        data.gsea <- data$V1
        names(data.gsea) <- as.character(data$GeneID)
        s2n.matrix.s <- sort(data.gsea, decreasing = T)
        #incProgress(0.25, detail = "Performing GSEA")
        ####----GSEA----####
        #perform GSEA
        gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, pvalueCutoff = input$userPval, eps = NA)
        #gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt(), verbose = F, pvalueCutoff = 1)
        #extract results and convert to tibble
        gsea.df <- gsea.res@result
        gsea.df
        #  incProgress(0.25, detail = "Complete!")
        #})
        gsea.df
        #}
      })
      
      #reactive to generate RData list
      #RDataListGen <- eventReactive(input$user.RData.Gen, {
      #  gmt <- GStable.ubg()
      #  colnames(gmt) <- c("term","gene")
      #  RData.u <- list()
      #  for (i in unique(gmt[,1])){
      #    RData.u[[i]] <- gmt[gmt[,1] == i,]$gene
      #  }
      #  RData.u
      #})
      
      #reactive for ssGSEA function
      ssGSEAfunc <- reactive({
        A <- A_raw()
        GS <- gs_react()
        scoreMethod <- input$ssGSEAtype
        #if (input$tables == 1) {
        #  GS <- gs()[(msigdb.gsea2()[input$msigdbTable_rows_selected,3])]
        #}
        #if (input$tables == 3) {
        #  GS <- gs2()[(GeneSet2()[input$tab2table_rows_selected,1])]
        #}
        #if (input$tables == 5) {
        #  gmt <- GStable.ubg()
        #  colnames(gmt) <- c("term","gene")
        #  gs_name <- user_gs_mirror()[input$GStable.u_rows_selected,1]
        #  gmt_sub <- gmt[which(gmt$term == gs_name),]
        #  GS <- list()
        #  for (i in unique(gmt_sub[,1])){
        #    GS[[i]] <- gmt_sub[gmt_sub[,1] == i,]$gene
        #  }
        #  #GS <- RDataListGen()[geneset_names]()[(user_gs_mirror()[input$GStable.u_rows_selected,1])]
        #}
        if (as.numeric(tools::file_path_sans_ext(packageVersion("GSVA"))) >= 1.5) {
          if (scoreMethod == "ssgsea") {
            ssGSEA_param <- GSVA::ssgseaParam(A,GS)
          } else if (scoreMethod == "gsva") {
            ssGSEA_param <- GSVA::gsvaParam(A,GS)
          } else if (scoreMethod == "plage") {
            ssGSEA_param <- GSVA::plageParam(A,GS)
          } else if (scoreMethod == "zscore") {
            ssGSEA_param <- GSVA::zscoreParam(A,GS)
          }
          ssGSEA <- GSVA::gsva(ssGSEA_param)
          #ssGSEA <- as.data.frame(t(ssGSEA))
          #ssGSEA[,SampleNameCol] <- rownames(ssGSEA)
        } else {
          ssGSEA <- gsva(A,GS,method = scoreMethod, verbose = FALSE)
          #ssGSEA <- as.data.frame(t(ssGSEA))
          #ssGSEA[,SampleNameCol] <- rownames(ssGSEA)
        }
        ssGSEA
        #GSVA::gsva(A, GS, method = input$ssGSEAtype, verbose = F)
      })
      
      #create background GMT from user input gene set table
      GStable.ubg <- reactive({
        gmt.u <- input$user.gmt.file
        ext <- tools::file_ext(gmt.u$datapath)
        req(gmt.u)
        headerCheck <- input$UserGSheaderCheck
        validate(need(ext == c("gmt","tsv","txt"), "Please upload .gmt, .tsv, or .txt file"))
        if (ext == "gmt") {
          read.gmt(gmt.u$datapath)
        }
        else {
          #as.data.frame(read_delim(gmt.u$datapath, delim = '\t', col_names = headerCheck))
          as.data.frame(fread(gmt.u$datapath))
        }
      })
      
      #gs mirror from user input for selection help on back end
      user_gs_mirror <- reactive({
        GeneSet <- as.data.frame(unique(GStable.ubg()[1]))
        rownames(GeneSet) <- 1:nrow(GeneSet)
        colnames(GeneSet)[1] <- "Gene_Set"
        GeneSet
      })
      
      #perform sig2noise calculation and create GSEA result from user chosen gene set
      datasetInput <- reactive({
        req(input$comparisonA)
        req(input$comparisonB)
        req(gmt_react())
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        
        if (input$volcanoCompChoice3 == "Two groups") {
          groupA <- input$comparisonA
          groupB <- input$comparisonB
        } else if (input$volcanoCompChoice3 == "One group") {
          groupA <- input$comparisonA_one
          groupB <- paste0("Not_",groupA)
          meta[,paste0(metacol,"_Dichot")] <- NA
          meta[which(meta[,metacol] == groupA),paste0(metacol,"_Dichot")] <- groupA
          meta[which(meta[,metacol] != groupA),paste0(metacol,"_Dichot")] <- paste0("Not_",groupA)
          metacol <- paste0(metacol,"_Dichot")
        }
        
        
        A <- A_raw()
        gmt <- gmt_react()
        if (isTruthy(groupA) & isTruthy(groupB)) {
          if (groupA != groupB) {
            if (length(which(meta[,metacol] == groupA)) > 1 | length(which(meta[,metacol] == groupB)) > 1) {
              groupA <- meta[which(meta[,metacol] == groupA),1]
              groupB <- meta[which(meta[,metacol] == groupB),1]
              ##----Signal-to-Noise Calculation----##
              A <- A + 0.00000001
              P = as.matrix(as.numeric(colnames(A) %in% groupA))
              n1 <- sum(P[,1])
              M1 <- A %*% P
              M1 <- M1/n1
              A2 <- A*A
              S1 <- A2 %*% P
              S1 <- S1/n1 - M1*M1 #
              S1 <- sqrt(abs((n1/(n1-1)) * S1))
              P = as.matrix(as.numeric(colnames(A) %in% groupB))
              n2 <- sum(P[,1])
              M2 <- A %*% P
              M2 <- M2/n2
              A2 <- A*A
              S2 <- A2 %*% P
              S2 <- S2/n2 - M2*M2
              S2 <- sqrt(abs((n2/(n2-1)) * S2))
              rm(A2)
              # small sigma "fix" as used in GeneCluster
              S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
              S2 <- ifelse(S2 == 0, 0.2, S2)
              S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
              S1 <- ifelse(S1 == 0, 0.2, S1)
              M1 <- M1 - M2
              rm(M2)
              S1 <- S1 + S2
              rm(S2)
              s2n.matrix <- M1/S1
              ##----Reformatting----##
              s2n.df <- as.data.frame(s2n.matrix)
              s2n.df$GeneID <- rownames(s2n.df)
              rownames(s2n.df) <- NULL
              data <- dplyr::select(s2n.df, GeneID, V1)
              data.gsea <- data$V1
              names(data.gsea) <- as.character(data$GeneID)
              s2n.matrix.s <- sort(data.gsea, decreasing = T)
              ##----GSEA----##
              GSEA(s2n.matrix.s, TERM2GENE = gmt,minGSSize = 5, maxGSSize  = 1000,
                   verbose = F, pvalueCutoff = 1, eps = NA)
              #if (input$tables == 1){
              #  GSEA(s2n.matrix.s, TERM2GENE = gmt()[which(gmt()[,1] == as.character(msigdb.gsea2()[input$msigdbTable_rows_selected,3])),],
              #       verbose = F, pvalueCutoff = input$userPval)
              #}
              #else if (input$tables == 3){
              #  GSEA(s2n.matrix.s, TERM2GENE = tab2()[which(tab2()[,1] == as.character(GeneSet2()[input$tab2table_rows_selected,1])),],
              #       verbose = F, pvalueCutoff = input$userPval)
              #}
              #else if (input$tables == 5){
              #  GSEA(s2n.matrix.s, TERM2GENE = GStable.ubg()[which(GStable.ubg()[,1] == as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])),],
              #       verbose = F, pvalueCutoff = input$userPval)
              #}
            }
            
          }
        }
        
        
      })
      
      
      
      #top genes data frame reactive
      topgenereact <- reactive({
        req(meta_react())
        req(metacol_reactVal())
        req(expr_react())
        req(input$volcanoCompChoice)
        #req((isTruthy(input$comparisonA2) && isTruthy(input$comparisonB2)) || isTruthy(input$comparisonA2_one))
        
        if (input$volcanoCompChoice == "DESeq2") {
          if (isTruthy(metacol_reactVal())) {
            desCol <- metacol_reactVal()
            expr <- expr_raw()
            meta <- meta_react()
            covars <- input$DEGCovarSelect
            meta <- meta[,c(colnames(meta)[1],desCol,covars)]
            meta <- meta[complete.cases(meta),]
            rownames(expr) <- expr[,1]
            expr <- expr[,-1]
            mat <- as.matrix(expr)
            sortedIndices <- match(colnames(mat), meta[,1])
            meta_deseq <- meta[sortedIndices,,drop=FALSE]
            rownames(meta_deseq) <- meta_deseq[,1]
            if (desCol != "Select column") {
              desRef <- input$DESeqDesignColRef_vol
              if (isTruthy(desRef)) {
                meta_deseq[,desCol] <- relevel(factor(meta_deseq[,desCol]), ref = desRef)
                colnames(meta_deseq) <- gsub(" ","_",colnames(meta_deseq))
                colnames(meta_deseq) <- gsub("[[:punct:]]","_",colnames(meta_deseq))
                desCol <- gsub(" ","_",desCol)
                desCol <- gsub("[[:punct:]]","_",desCol)
                covars <- gsub(" ","_",covars)
                covars <- gsub("[[:punct:]]","_",covars)
                
                if (isTruthy(covars)) {
                  meta_deseq[,covars] <- lapply(meta_deseq[,covars,drop = F], factor)
                  desCol <- paste(c(covars, desCol),collapse = "+")
                }
                
                
                dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat,
                                                      colData = meta_deseq,
                                                      design = as.formula(paste0("~ ",desCol)))
                
                
                dds <- DESeq(dds)
                res <- results(dds)
                print(as.formula(paste0("~ ",desCol)));
                
                resOrdered <- res[order(res$padj),]
                top1 <- as.data.frame(resOrdered)
                
                top2 <- top1
                colnames(top2)[c(1,2,5,6)] <- c("AveExpr","logFC","P.Value","adj.P.Val")
                top2["GeneName"] <- rownames(top2)
                top2['group'] <- "NotSignificant"
                p_cut <- input$p_cutoff
                f_cut <- input$fc_cutoff
                top2[which(top2$P.Value < p_cut & abs(top2$logFC) < abs(f_cut)), "group"] <- "Significant"
                top2[which(top2$P.Value > p_cut & abs(top2$logFC) > abs(f_cut)), "group"] <- "FoldChange"
                top2[which(top2$P.Value < p_cut & abs(top2$logFC) > abs(f_cut)), "group"] <- "Significant&FoldChange"
                top2['FCgroup'] <- "NotSignificant"
                top2[which(abs(top2$logFC) > abs(f_cut)), "group2"] <- "FoldChange"
                top2
              }
            }
          }
        } else {
          meta <- meta_react()
          metacol <- metacol_reactVal()
          expr <- expr_react()
          covars <- input$DEGCovarSelect
          if (length(colnames(meta)) == 2) {
            metacol <- colnames(meta)[2]
          }
          
          if (all(c(colnames(meta)[1],metacol,covars) %in% colnames(meta))) {
            metaSub <- meta[,c(colnames(meta)[1],metacol,covars)]
            metaSub <- metaSub[complete.cases(metaSub),]
            
            
            if (input$volcanoCompChoice == "Limma: Two groups") {
              groupA <- input$comparisonA2
              groupB <- input$comparisonB2
              A <- metaSub[which(metaSub[,metacol] == groupA),1]
              B <- metaSub[which(metaSub[,metacol] == groupB),1]
            } else if (input$volcanoCompChoice == "Limma: One group") {
              groupA <- input$comparisonA2_one
              groupB <- paste0("Not_",groupA)
              metaSub[,paste0(metacol,"_Dichot")] <- NA
              metaSub[which(metaSub[,metacol] == groupA),paste0(metacol,"_Dichot")] <- groupA
              metaSub[which(metaSub[,metacol] != groupA),paste0(metacol,"_Dichot")] <- paste0("Not_",groupA)
              A <- metaSub[which(metaSub[,paste0(metacol,"_Dichot")] == groupA),1]
              B <- metaSub[which(metaSub[,paste0(metacol,"_Dichot")] == paste0("Not_",groupA)),1]
              metacol <- paste0(metacol,"_Dichot")
            }
            
            
            
            metaSub <- metaSub[which(metaSub[,1] %in% c(A,B)),]
            colnames(metaSub) <- gsub(" ","_",colnames(metaSub))
            colnames(metaSub) <- gsub("[[:punct:]]","_",colnames(metaSub))
            metacol <- gsub(" ","_",metacol)
            metacol <- gsub("[[:punct:]]","_",metacol)
            covars <- gsub(" ","_",covars)
            covars <- gsub("[[:punct:]]","_",covars)
            metaSub[,metacol] <- factor(metaSub[,metacol], levels = c(groupA,groupB))
            
            mat <- expr[,metaSub[,1]]
            mat <- log2(mat + 1.0)
            if (isTruthy(covars)) {
              metaSub[,covars] <- lapply(metaSub[,covars,drop = F], factor)
              form <- as.formula(paste0("~0 +",paste(c(metacol,covars),collapse = "+")))
            } else {
              #metaSub[,metacol] <- factor(metaSub[,metacol])
              form <- as.formula(paste0("~0 +",metacol))
            }
            
            if (ncol(mat) > 0) {
              
              if(nlevels(metaSub[,metacol]) >= 2){
                designA <- eval(substitute(model.matrix(form, data = metaSub)))
                fit <- lmFit(mat, design = designA)
                colnames(designA) <- gsub(" ","_",colnames(designA))
                colnames(designA) <- gsub("[[:punct:]]","_",colnames(designA))
                contrast.matrix <- makeContrasts(contrasts = paste0(colnames(designA)[1],"-",colnames(designA)[2]), levels = designA)
                fit2 <- contrasts.fit(fit, contrast.matrix)
                fit2 <- eBayes(fit2)
                options(digits = 4)
                top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
                top2 <- top1
                
                top2["GeneName"] <- rownames(top2)
                top2['group'] <- "NotSignificant"
                p_cut <- input$p_cutoff
                f_cut <- input$fc_cutoff
                top2[which(top2$P.Value < p_cut & abs(top2$logFC) < abs(f_cut)), "group"] <- "Significant"
                top2[which(top2$P.Value > p_cut & abs(top2$logFC) > abs(f_cut)), "group"] <- "FoldChange"
                top2[which(top2$P.Value < p_cut & abs(top2$logFC) > abs(f_cut)), "group"] <- "Significant&FoldChange"
                top2['FCgroup'] <- "NotSignificant"
                top2[which(abs(top2$logFC) > abs(f_cut)), "group2"] <- "FoldChange"
                top2
              }
            }
            
          }
          
        }
        
        
        
        
        
      })
      
      topgenereact2 <- reactive({
        
        req(metacol_reactVal())
        withProgress(message = "Processing", value = 0, {
          incProgress(0.5, detail = "Calculating Differential Expression")
          
          if (input$volcanoCompChoice4 == "DESeq2") {
            if (isTruthy(metacol_reactVal())) {
              desCol <- metacol_reactVal()
              expr <- expr_raw()
              meta <- meta_react()
              covars <- input$DEGCovarSelectHeat
              meta <- meta[,c(colnames(meta)[1],desCol,covars)]
              meta <- meta[complete.cases(meta),]
              rownames(expr) <- expr[,1]
              expr <- expr[,-1]
              mat <- as.matrix(expr)
              sortedIndices <- match(colnames(mat), meta[,1])
              meta_deseq <- meta[sortedIndices,,drop=FALSE]
              rownames(meta_deseq) <- meta_deseq[,1]
              if (desCol != "Select column") {
                desRef <- input$DESeqDesignColRef_heat
                if (isTruthy(desRef)) {
                  meta_deseq[,desCol] <- relevel(factor(meta_deseq[,desCol]), ref = desRef)
                  colnames(meta_deseq) <- gsub(" ","_",colnames(meta_deseq))
                  colnames(meta_deseq) <- gsub("[[:punct:]]","_",colnames(meta_deseq))
                  desCol <- gsub(" ","_",desCol)
                  desCol <- gsub("[[:punct:]]","_",desCol)
                  covars <- gsub(" ","_",covars)
                  covars <- gsub("[[:punct:]]","_",covars)
                  
                  if (isTruthy(covars)) {
                    meta_deseq[,covars] <- lapply(meta_deseq[,covars,drop = F], factor)
                    desCol <- paste(c(covars, desCol),collapse = "+")
                  }
                  dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat,
                                                        colData = meta_deseq,
                                                        design = as.formula(paste0("~ ",desCol)))
                  dds <- DESeq(dds)
                  res <- results(dds)
                  resOrdered <- res[order(res$padj),]
                  top1 <- as.data.frame(resOrdered)
                }
              }
            }
          } else {
            req(input$comparisonA2_h)
            req(input$comparisonB2_h)
            meta <- meta_react()
            metacol <- input$DEGcolHeat
            expr <- expr_react()
            covars <- input$DEGCovarSelectHeat
            if (length(colnames(meta)) == 2) {
              metacol <- colnames(meta)[2]
            }
            metaSub <- meta[,c(colnames(meta)[1],metacol,covars)]
            metaSub <- metaSub[complete.cases(metaSub),]
            
            if (input$volcanoCompChoice4 == "Limma: Two groups") {
              A_choice <- input$comparisonA2_h            #Comparison group A
              B_choice <- input$comparisonB2_h            #Comparison group B
              A <- metaSub[which(metaSub[,metacol] == A_choice),1]
              B <- metaSub[which(metaSub[,metacol] == B_choice),1]
            } else if (input$volcanoCompChoice4 == "Limma: One group") {
              A_choice <- input$comparisonA2_h_one
              B_choice <- paste0("Not_",A_choice)
              metaSub[,paste0(metacol,"_Dichot")] <- NA
              metaSub[which(metaSub[,metacol] == A_choice),paste0(metacol,"_Dichot")] <- A_choice
              metaSub[which(metaSub[,metacol] != A_choice),paste0(metacol,"_Dichot")] <- paste0("Not_",A_choice)
              A <- metaSub[which(metaSub[,paste0(metacol,"_Dichot")] == A_choice),1]
              B <- metaSub[which(metaSub[,paste0(metacol,"_Dichot")] == paste0("Not_",A_choice)),1]
              metacol <- paste0(metacol,"_Dichot")
            }
            
            
            metaSub <- metaSub[which(metaSub[,1] %in% c(A,B)),]
            colnames(metaSub) <- gsub(" ","_",colnames(metaSub))
            colnames(metaSub) <- gsub("[[:punct:]]","_",colnames(metaSub))
            metacol <- gsub(" ","_",metacol)
            metacol <- gsub("[[:punct:]]","_",metacol)
            covars <- gsub(" ","_",covars)
            covars <- gsub("[[:punct:]]","_",covars)
            metaSub[,metacol] <- factor(metaSub[,metacol], levels = c(A_choice,B_choice))
            
            mat <- expr[,metaSub[,1]]
            mat <- log2(mat + 1.0)
            if (isTruthy(covars)) {
              metaSub[,covars] <- lapply(metaSub[,covars,drop = F], factor)
              form <- as.formula(paste0("~0 +",paste(c(metacol,covars),collapse = "+")))
            } else {
              #metaSub[,metacol] <- factor(metaSub[,metacol])
              form <- as.formula(paste0("~0 +",metacol))
            }
            
            if (nlevels(metaSub[,metacol]) >= 2) {
              
              designA <- eval(substitute(model.matrix(form, data = metaSub)))
              fit <- lmFit(mat, design = designA)
              colnames(designA) <- gsub(" ","_",colnames(designA))
              colnames(designA) <- gsub("[[:punct:]]","_",colnames(designA))
              contrast.matrix <- makeContrasts(contrasts = paste0(colnames(designA)[1],"-",colnames(designA)[2]), levels = designA)
              fit2 <- contrasts.fit(fit, contrast.matrix)
              fit2 <- eBayes(fit2)
              options(digits = 4)
              top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
              top2 <- top1
            }
          }
          
          incProgress(0.5, detail = "Complete!")
        })
        
        if (exists("top1")) {
          return(top1)
        }
        
      })
      
      ####----Data Exploration----####
      
      ####----MVG Heatmap----####
      
      heat_colAnn <- reactive({
        
        if (length(input$GroupColMulti) > 0) {
          #if (!is.null(input$GroupColMulti)) {
          #if (isTruthy(input$GroupColMulti)) {
          meta <- meta_react()
          rownames(meta) <- meta[,1]
          anno_cols <- input$GroupColMulti
          
          meta_sub <- meta[,anno_cols, drop = F]
          
          
          meta_sub[,anno_cols] <- as.data.frame(lapply(meta_sub[,anno_cols, drop = F], function(x) {
            if (!is.numeric(x) | all(suppressWarnings(as.numeric(x))%%1==0)) {
              return(factor(x))
            } else {
              return(x)
            }
          }))
          
          #meta_sub[,input$GroupColMulti] <- as.data.frame(lapply(meta_sub[,input$GroupColMulti, drop = F], factor))
          colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                      name = anno_cols,
                                                      which = 'col'
          )
          
          colAnn
        }  else {
          meta <- meta_react()
          if (ncol(meta) == 2) {
            rownames(meta) <- meta[,1]
            meta_sub <- meta[,2, drop = F]
            meta_sub[,1] <- as.factor(meta_sub[,1])
            colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                        name = colnames(meta_sub)[2],
                                                        which = 'col'
            )
            colAnn
          } else {
            colAnn <- NULL
            colAnn
          }
        }
        
      })
      
      MVG_heat_data <- reactive({
        
        req(expr_react())
        withProgress(message = "Processing", value = 0, {
          incProgress(0.25, detail = "Calculating Variance")
          top_probes <- input$NumFeatures
          clust_method <- input$ClusteringMethod
          var_type <- input$VarianceMeasure
          expr <- expr_react()
          exp <- expr
          mad <- NULL
          var <- NULL
          cv <- NULL
          if (var_type == "MAD"){
            mad <- apply(log2(exp + 1), 1, mad)
            mad <- sort(mad, decreasing = T)
            mad <- head(mad, n = (top_probes +1))
            out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
            colnames(out) <- c("Gene", "MAD", colnames(exp))
            dataset <- exp[names(mad),]
          }
          else if (var_type == "VAR"){
            var <- apply(log2(exp + 1), 1, var)
            var <- sort(var, decreasing = T)
            var <- head(var, n = (top_probes +1))
            out <- cbind(names(var), var[names(var)], exp[names(var),])
            colnames(out) <- c("Gene", "VAR", colnames(exp))
            dataset <- exp[names(var),]
          }
          else if (var_type == "CV"){
            cv <- apply(log2(exp + 1), 1, cv)
            cv <- sort(cv, decreasing = T)
            cv <- head(cv, n = (top_probes +1))
            out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
            colnames(out) <- c("Gene", "CV", colnames(exp))
            dataset <- exp[names(cv),]
          }
          incProgress(0.25, detail = "Calculating Z-Score")
          dataset <- log2(dataset + 1)
          zdataset <- apply(dataset, 1, scale)
          zdataset <- apply(zdataset, 1, rev)
          colnames(zdataset) <- colnames(dataset)
          dataset <- as.matrix(zdataset)
          dataset[is.na(dataset)] <- 0
          dataset = dataset[apply(dataset, 1, function(x) !all(x==0)),]
          incProgress(0.5, detail = "Complete!")
        })
        
        dataset
        
        
      })
      
      ## MVG heatmap reactive
      MVGheatmap_react <- reactive ({
        
        req(input$ClusteringMethod)
        withProgress(message = "Processing", value = 0, {
          incProgress(0.5, detail = "Generating Heatmap")
          row_names_choice <- input$ShowRowNames1
          col_names_choice <- input$ShowColNames1
          row_font <- input$heatmapFont2.r
          col_font <- input$heatmapFont2.c
          clust_method <- input$ClusteringMethod
          kmeans <- ifelse(is.na(input$NumClusters),1,input$NumClusters)
          clust_cols_opt <- input$clustcolsMVG
          clust_rows_opt <- input$clustrowsMVG
          colAnno <- heat_colAnn()
          dataset <- MVG_heat_data()
          p <- suppressMessages(ComplexHeatmap::Heatmap(dataset,
                                                        top_annotation = colAnno,
                                                        clustering_method_rows = clust_method, km = kmeans,
                                                        show_row_names = row_names_choice, show_column_names = col_names_choice,
                                                        cluster_rows = clust_rows_opt, cluster_columns = clust_cols_opt,
                                                        row_names_gp = gpar(fontsize = row_font), column_names_gp = gpar(fontsize = col_font),
                                                        heatmap_legend_param = list(title = "Expression"),
                                                        #width = ncol(dataset)*unit(3, "mm"), 
                                                        #height = nrow(dataset)*unit(2.5, "mm"),
                                                        border = F))
          incProgress(0.5, detail = "Complete!")
        })
        draw(p, padding = unit(c(50, 50, 2, 2), "mm")) # unit(c(bottom,left,right,top))
        
        
      })
      
      ####----Custom Heatmap----####
      
      Cutsom_heat_colAnn <- reactive({
        
        if (isTruthy(input$GroupColMulti)) {
          meta <- meta_react()
          usersamps <- input$userheatsamp2
          meta <- meta[which(meta[,1] %in% usersamps),]
          meta <- meta[match(usersamps,meta[,1]),]
          rownames(meta) <- meta[,1]
          anno_cols <- input$GroupColMulti
          meta_sub <- meta[,anno_cols, drop = F]
          meta_sub[,anno_cols] <- as.data.frame(lapply(meta_sub[,anno_cols, drop = F], function(x) {
            if (!is.numeric(x) | all(suppressWarnings(as.numeric(x))%%1==0)) {
              return(factor(x))
            } else {
              return(x)
            }
          }))
          #meta_sub[,input$GroupColMulti] <- as.data.frame(lapply(meta_sub[,input$GroupColMulti, drop = F], factor))
          colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                      name = anno_cols,
                                                      which = 'col'
          )
          colAnn
        } else {
          meta <- meta_react()
          if (ncol(meta) == 2) {
            rownames(meta) <- meta[,1]
            usersamps <- input$userheatsamp2
            meta <- meta[which(meta[,1] %in% usersamps),]
            meta <- meta[match(usersamps,meta[,1]),]
            rownames(meta) <- meta[,1]
            meta_sub <- meta[,2, drop = F]
            meta_sub[,1] <- as.factor(meta_sub[,1])
            colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                        name = colnames(meta_sub)[2],
                                                        which = 'col'
            )
            colAnn
          } else {
            colAnn <- NULL
            colAnn
          }
        }
        
      })
      
      Custom_heat_data <- reactive({
        genelist.uih <- unlist(strsplit(input$userheatgenes, " "))
        if (length(input$heatmapGeneSelec) >= 2 | length(genelist.uih) > 1) {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.25, detail = "Customizing")
            expr <- expr_react()
            meta <- meta_react()
            genelist.uih <- NULL
            genelist.ush <- NULL
            genelist.uih2 <- NULL
            genelist.ush <- input$heatmapGeneSelec
            genelist.uih <- unlist(strsplit(input$userheatgenes, " "))
            #genelist.uih2 <- unlist(strsplit(input$userheatgenes, "\t"))
            heatgenes <- c(genelist.ush,genelist.uih,genelist.uih2)
            heatgenes <- heatgenes[!is.na(heatgenes)]
            heatgenes <- heatgenes[!is.null(heatgenes)]
            usersamps <- input$userheatsamp2
            exp <- expr[heatgenes,which(colnames(expr) %in% usersamps)]
            exp <- exp[,match(usersamps,colnames(exp))]
            #meta <- meta[which(meta[,1] %in% usersamps),]
            #meta <- meta[match(colnames(usersamps),meta[,1]),]
            incProgress(0.25, detail = "Calculating Z-Score")
            dataset <- exp
            dataset <- log2(dataset + 1)
            zdataset <- apply(dataset, 1, scale)
            zdataset <- apply(zdataset, 1, rev)
            colnames(zdataset) <- colnames(dataset)
            dataset <- as.matrix(zdataset)
            dataset[is.na(dataset)] <- 0
            dataset = dataset[apply(dataset, 1, function(x) !all(x==0)),]
            dataset <- dataset[match(heatgenes, rownames(dataset)),]
            dataset <- dataset[complete.cases(dataset),]
            incProgress(0.5, detail = "Complete!")
          })
          
          dataset
        }
        
      })
      
      #Custom heat map reactive
      CustomHeatmap_react <- reactive({
        
        withProgress(message = "Processing", value = 0, {
          incProgress(0.5, detail = "Generating Heatmap")
          row_names_choice <- input$ShowRowNames2
          col_names_choice <- input$ShowColNames2
          clust_cols_opt <- input$ClusColOptCust
          clust_rows_opt <- input$ClusRowOptCust
          row_font <- input$heatmapFont3.r
          col_font <- input$heatmapFont3.c
          kmeans <- ifelse(is.na(input$NumClusters),1,input$NumClusters)
          clust_method <- input$ClusteringMethodCust
          colAnno <- Cutsom_heat_colAnn()
          dataset <- Custom_heat_data()
          
          p <- suppressMessages(ComplexHeatmap::Heatmap(dataset,
                                                        top_annotation = colAnno,
                                                        clustering_method_rows = clust_method, km = kmeans,
                                                        show_row_names = row_names_choice, show_column_names = col_names_choice,
                                                        cluster_rows = clust_rows_opt, cluster_columns = clust_cols_opt,
                                                        row_names_gp = gpar(fontsize = row_font), column_names_gp = gpar(fontsize = col_font),
                                                        heatmap_legend_param = list(title = "Expression"),
                                                        border = F))
          incProgress(0.5, detail = "Complete!")
        })
        draw(p, padding = unit(c(50, 50, 2, 2), "mm")) # unit(c(bottom,left,right,top))
        
      })
      
      ####----DEG Heatmap----####
      
      heat_colAnn_DEG <- reactive({
        
        if (length(input$GroupColMulti) > 0) {
          #if (!is.null(input$GroupColMulti)) {
          #if (isTruthy(input$GroupColMulti)) {
          meta <- meta_react()
          rownames(meta) <- meta[,1]
          
          
          metacol <- input$DEGcolHeat
          if (length(colnames(meta)) == 2) {
            metacol <- colnames(meta)[2]
          }
          
          if (input$volcanoCompChoice4 == "Limma: Two groups" | input$volcanoCompChoice4 == "DESeq2") {
            groupA <- input$comparisonA2_h
            groupB <- input$comparisonB2_h
            anno_cols <- input$GroupColMulti
            meta_sub <- meta[,anno_cols, drop = F]
            meta_sub[,anno_cols] <- as.data.frame(lapply(meta_sub[,anno_cols, drop = F], function(x) {
              if (!is.numeric(x) | all(suppressWarnings(as.numeric(x))%%1==0)) {
                return(factor(x))
              } else {
                return(x)
              }
            }))
            #meta_sub[,input$GroupColMulti] <- as.data.frame(lapply(meta_sub[,input$GroupColMulti, drop = F], factor))
            colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                        name = anno_cols,
                                                        which = 'col'
            )
          } else if (input$volcanoCompChoice4 == "Limma: One group") {
            groupA <- input$comparisonA2_h_one
            groupB <- paste0(unique(meta[which(meta[,metacol] != groupA),metacol]), collapse = " - ")
            #groupB <- paste0("Not_",groupA)
            meta[,paste0(metacol,"_Dichot")] <- NA
            meta[which(meta[,metacol] == groupA),paste0(metacol,"_Dichot")] <- groupA
            #metaSub_noNA[which(metaSub_noNA[,metacol] != groupA),paste0(metacol,"_Dichot")] <- paste0("Not_",groupA)
            metaSub_noNA[which(metaSub_noNA[,metacol] != groupA),paste0(metacol,"_Dichot")] <- groupB
            metacol <- paste0(metacol,"_Dichot")
            
            anno_cols <- c(metacol,input$GroupColMulti)
            meta_sub <- meta[,anno_cols, drop = F]
            meta_sub[,anno_cols] <- as.data.frame(lapply(meta_sub[,anno_cols, drop = F], function(x) {
              if (!is.numeric(x) | all(suppressWarnings(as.numeric(x))%%1==0)) {
                return(factor(x))
              } else {
                return(x)
              }
            }))
            #meta_sub[,input$GroupColMulti] <- as.data.frame(lapply(meta_sub[,input$GroupColMulti, drop = F], factor))
            colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                        name = anno_cols,
                                                        which = 'col'
            )
            
            
            
            
            
            #meta_sub <- meta[,c(metacol,input$GroupColMulti), drop = F]
            #meta_sub[,c(metacol,input$GroupColMulti)] <- as.data.frame(lapply(meta_sub[,c(metacol,input$GroupColMulti), drop = F], factor))
            #colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
            #                                            name = c(metacol,input$GroupColMulti),
            #                                            which = 'col'
            #)
          }
          
          colAnn
        }  else {
          meta <- meta_react()
          if (ncol(meta) == 2) {
            rownames(meta) <- meta[,1]
            meta_sub <- meta[,2, drop = F]
            meta_sub[,1] <- as.factor(meta_sub[,1])
            colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                        name = colnames(meta_sub)[2],
                                                        which = 'col'
            )
            colAnn
          } else {
            colAnn <- NULL
            colAnn
          }
        }
        
      })
      
      DEG_heat_data <- reactive({
        req(topgenereact2())
        withProgress(message = "Processing", value = 0, {
          incProgress(0.5, detail = "Calculating Z-Score")
          # UI Inputs
          meta <- meta_react()
          expr <- expr_react()
          FC_cutoff <- input$fc_cutoff_h              #FC cutoff for top gene selection
          P_cutoff <- input$p_cutoff_h                #P-value cutoff for top gene selections
          top_probes <- input$top_x_h                 #Number of top genes to show on heatmap
          top1 <- topgenereact2()
          if (input$volcanoCompChoice4 == "DESeq2") {
            top_above_cutoff <- top1[which(top1$log2FoldChange > abs(FC_cutoff) & top1$pvalue < P_cutoff),]
          } else {
            top_above_cutoff <- top1[which(top1$logFC > abs(FC_cutoff) & top1$P.Value < P_cutoff),]
          }
          
          # Get gene list from Top table
          genelist <- rownames(top_above_cutoff)[c(1:top_probes)]
          
          # Heatmap Calculations
          dataset <- expr[which(rownames(expr) %in% genelist),]
          if (length(dataset) > 0) {
            dataset <- log2(dataset + 1)
            zdataset <- apply(dataset, 1, scale)
            zdataset <- apply(zdataset, 1, rev)
            colnames(zdataset) <- colnames(dataset)
            dataset <- as.matrix(zdataset)
            dataset[is.na(dataset)] <- 0
            dataset = dataset[apply(dataset, 1, function(x) !all(x==0)),]
          }
          incProgress(0.5, detail = "Complete!")
        })
        dataset
      })
      
      # DEG heatmap reactive
      DEGHeatmap_react <- reactive({
        
        req(DEG_heat_data)
        withProgress(message = "Processing", value = 0, {
          incProgress(0.5, detail = "Generating Heatmap")
          ## UI Inputs
          row_names_choice <- input$ShowRowNames3     #choose to show row names or not
          col_names_choice <- input$ShowColNames3
          kmeans <- ifelse(is.na(input$NumClusters),1,input$NumClusters)
          clust_method <- input$ClusteringMethod_degh #Cluster method for heatmap
          col_font <- input$heatmapFont3.c.deg        #Heatmap column font size
          row_font <- input$heatmapFont3.r.deg        #Heatmap row font size
          #color_choice <- input$ColorPalette3
          clust_rows_opt <- input$ClusRowOptDEG
          clust_cols_opt <- input$ClusColOptDEG
          
          colAnno <- heat_colAnn_DEG()
          dataset <- DEG_heat_data()
          
          
          p <- suppressMessages(ComplexHeatmap::Heatmap(dataset,
                                                        top_annotation = colAnno,
                                                        clustering_method_rows = clust_method, km = kmeans,
                                                        show_row_names = row_names_choice, show_column_names = col_names_choice,
                                                        cluster_rows = clust_rows_opt, cluster_columns = clust_cols_opt,
                                                        row_names_gp = gpar(fontsize = row_font), column_names_gp = gpar(fontsize = col_font),
                                                        heatmap_legend_param = list(title = "Expression"),
                                                        border = F))
          incProgress(0.5, detail = "Complete!")
        })
        
        draw(p, padding = unit(c(50, 50, 2, 2), "mm")) # unit(c(bottom,left,right,top))
        
      })
      
      ####----Avg Expr Heatmap----####
      
      
      AvgExprHeatmap_react <- reactive({
        
        if (length(input$avgheatmapGeneSelec) >= 2 || length(input$avguserheatgenes) > 1) {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.25, detail = "Customizing")
            meta <- meta_react()
            metacol <- metacol_reactVal()
            if (length(colnames(meta)) == 2) {
              metacol <- colnames(meta)[2]
            }
            
            expr <- expr_react()
            row_names_choice <- input$ShowRowNames5      #Chose to show row names or not
            col_names_choice <- input$ShowColNames5
            clust_row_opt <- input$AvgClusRowOptCust
            clust_col_opt <- input$AvgClusColOptCust
            group_choices <- input$SampCondSelectionCust
            col_font <- input$heatmapFont3.c.deg.avg.cust     #Heatmap column font size
            row_font <- input$heatmapFont3.r.deg.avg.cust     #Heatmap row font size
            clust_method <- input$AvgClusteringMethodCust
            genelist.uih <- NULL
            genelist.ush <- NULL
            genelist.uih2 <- NULL
            genelist.ush <- input$avgheatmapGeneSelec
            genelist.uih <- unlist(strsplit(input$avguserheatgenes, " "))
            #genelist.uih2 <- unlist(strsplit(input$userheatgenes, "\t"))
            heatgenes <- c(genelist.ush,genelist.uih,genelist.uih2)
            heatgenes <- heatgenes[!is.na(heatgenes)]
            heatgenes <- heatgenes[!is.null(heatgenes)]
            req(length(heatgenes) > 2)
            incProgress(0.25, detail = "Averaging Groups")
            AvgExprDF <- data.frame(rownames(expr))
            for (i in group_choices) {
              samples <- meta[which(meta[,metacol] == i),1]
              if (length(samples) <= 1) {
                AvgExprDF[,paste("AvgExpr_",i, sep = "")] <- expr[,samples]
              }
              else if (length(samples) > 1) {
                AvgExprDF[,paste("AvgExpr_",i, sep = "")] <- rowMeans(expr[,samples], na.rm = T)
              }
            }
            rownames(AvgExprDF) <- AvgExprDF[,1]
            AvgExprDF <- AvgExprDF[,-1]
            getGenes <- intersect(rownames(AvgExprDF),heatgenes)
            AvgExprDF <- AvgExprDF[getGenes,]
            dataset <- AvgExprDF
            incProgress(0.25, detail = "Calculating Z-Score")
            #if (length(group_choices) > 0) {
            dataset <- log2(dataset + 1)
            zdataset <- apply(dataset, 1, scale)
            zdataset <- apply(zdataset, 1, rev)
            colnames(zdataset) <- colnames(dataset)
            dataset <- as.matrix(zdataset)
            dataset[is.na(dataset)] <- 0
            anno_meta <- data.frame(colnames(AvgExprDF))
            anno_meta$Type <- gsub("AvgExpr_","",anno_meta[,1])
            rownames(anno_meta) <- anno_meta[,1]
            anno_meta <- anno_meta[,-1,drop = F]
            dataset = dataset[apply(dataset, 1, function(x) !all(x==0)),]
            minimum = -5;
            maximum = 5;
            if (abs(min(dataset)) > abs(max(dataset))) {
              dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
            } else {
              dataset[dataset > abs(min(dataset))] = abs(min(dataset))
            }
            bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
            #Heatmap color
            color_choice <- input$ColorPalette5
            col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
            if (color_choice == "original") {
              HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
              hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
            }
            else if (color_choice %in% col_sets) {
              HeatMap_Colors <- brewer.pal(n = 5, color_choice)
              hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
            }
            else if (color_choice == "Inferno") {
              hmcols <- inferno(500)
            }
            else if (color_choice == "Viridis") {
              hmcols <- viridis(500)
            }
            else if (color_choice == "Plasma") {
              hmcols <- plasma(500)
            }
            else if (color_choice == "OmniBlueRed") {
              hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
            }
            else if (color_choice == "LightBlueBlackRed") {
              hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
            }
            else if (color_choice == "GreenBlackRed") {
              hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
            }
            dataset <- dataset[match(heatgenes, rownames(dataset)),]
            dataset <- dataset[complete.cases(dataset),]
            hm <- pheatmap::pheatmap(as.matrix(dataset),
                                     cluster_col = clust_col_opt,
                                     cluster_row = clust_row_opt,
                                     fontsize_row = row_font,
                                     fontsize_col = col_font,
                                     show_rownames = row_names_choice,
                                     show_colnames = col_names_choice,
                                     clustering_method = clust_method,
                                     color=hmcols,
                                     annotation_col = anno_meta,
                                     angle_col = 90,
                                     border_color = NA)
            hm
            #}
            incProgress(0.25, detail = "Complete!")
          })
          
          
        }
        
        
      })
      
      #meta_factored <- apply(meta,2,function(x) factor(x))
      #meta_factored <- as.data.frame(lapply(meta, factor))
      
      observe({
        
        req(meta_react())
        req(input$ssGSEAGroupColMulti)
        meta <- meta_react()
        metacol <- input$ssGSEAGroupColMulti[1]
        meta2 <- meta[,c(colnames(meta)[1],metacol)]
        meta2[,2] <- as.factor(meta2[,2])
        meta2 <- meta2[order(meta2[,2]),]
        
        updateSelectizeInput(session,"userheatsampSS",choices = meta2[,1], selected = meta2[,1])
      })
      
      ssgseaheatmap_df_react <- reactive({
        
        meta <- meta_react()
        dataset <- ssgsea_heat_data()
        if (isTruthy(input$ssGSEAGroupColMulti)) {
          meta <- meta[which(meta[,1] %in% colnames(dataset)),]
          meta_sub <- meta[,c(colnames(meta)[1],input$ssGSEAGroupColMulti), drop = F]
          dataset_t <- as.data.frame(t(dataset))
          dataset_t <- cbind(SampleName = rownames(dataset_t),dataset_t)
          heatdf <- merge(meta_sub,dataset_t, all.y = T, by.x = colnames(meta_sub)[1], by.y = "SampleName")
          heatdf <- heatdf[match(dataset_t[,1],heatdf[,1]),]
          heatdf
        } else {
          if (ncol(meta) == 2) {
            meta <- meta[which(meta[,1] %in% colnames(dataset)),]
            meta_sub <- meta
            dataset_t <- as.data.frame(t(dataset))
            dataset_t <- cbind(SampleName = rownames(dataset_t),dataset_t)
            heatdf <- merge(meta_sub,dataset_t, all.y = T, by.x = colnames(meta_sub)[1], by.y = "SampleName")
            heatdf <- heatdf[match(dataset_t[,1],heatdf[,1]),]
            heatdf
          } else {
            heatdf <- as.data.frame(t(dataset))
            heatdf <- cbind(SampleName = rownames(heatdf),heatdf)
            heatdf
          }
        }
        
      })
      
      output$ssgseaheatmap_df <- DT::renderDataTable({
        
        df <- ssgseaheatmap_df_react()
        dataset <- ssgsea_heat_data()
        DT::datatable(df, 
                      options = list(lengthMenu = c(5,10,20,50,100,1000),
                                     pageLength = 10,
                                     scrollX = TRUE),
                      rownames = FALSE)%>%
          formatRound(columns = rownames(dataset), digits = 4)
        
      })
      
      ssgsea_heat_anno <- reactive({
        
        if (isTruthy(input$ssGSEAGroupColMulti)) {
          meta <- meta_react()
          dataset <- ssgsea_heat_data()
          meta <- meta[which(meta[,1] %in% colnames(dataset)),]
          meta <- meta[match(colnames(dataset),meta[,1]),]
          rownames(meta) <- meta[,1]
          meta_sub <- meta[,input$ssGSEAGroupColMulti, drop = F]
          
          
          anno_cols <- input$ssGSEAGroupColMulti
          meta_sub <- meta[,anno_cols, drop = F]
          meta_sub[,anno_cols] <- as.data.frame(lapply(meta_sub[,anno_cols, drop = F], function(x) {
            if (!is.numeric(x) | all(suppressWarnings(as.numeric(x))%%1==0)) {
              return(factor(x))
            } else {
              return(x)
            }
          }))
          #meta_sub[,input$input$ssGSEAGroupColMulti] <- as.data.frame(lapply(meta_sub[,input$input$ssGSEAGroupColMulti, drop = F], factor))
          colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                      name = anno_cols,
                                                      which = 'col'
          )
          
          
          
          #meta_sub[,input$ssGSEAGroupColMulti] <- as.data.frame(lapply(meta_sub[,input$ssGSEAGroupColMulti, drop = F], factor))
          #
          #colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
          #                                            name = input$ssGSEAGroupColMulti,
          #                                            which = 'col'
          #)
          colAnn
        } else {
          meta <- meta_react()
          dataset <- ssgsea_heat_data()
          if (ncol(meta) == 2) {
            dataset <- ssgsea_heat_data()
            meta <- meta[which(meta[,1] %in% colnames(dataset)),]
            meta <- meta[match(colnames(dataset),meta[,1]),]
            rownames(meta) <- meta[,1]
            dataset <- ssgsea_heat_data()
            meta <- meta[which(meta[,1] %in% colnames(dataset)),]
            meta_sub <- meta[,2, drop = F]
            meta_sub[,1] <- as.factor(meta_sub[,1])
            colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                        name = colnames(meta_sub)[2],
                                                        which = 'col'
            )
            colAnn
          } else {
            colAnn <- NULL
            colAnn
          }
        }
        
      })
      
      
      ssgsea_heat_data <- reactive({
        
        req(meta_react())
        req(metacol_reactVal())
        req(input$ssGSEAtype)
        req(input$userheatsampSS)
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        A <- A_raw()
        scoreMethod <- input$ssGSEAtype
        
        samples_chosen <- input$userheatsampSS
        
        #A <- A[,samples_chosen]
        A <- A[,which(colnames(A) %in% samples_chosen)]
        meta <- meta[which(meta[,1] %in% samples_chosen),]
        GS <- gs_react_heat()
        geneset_names <- names(GS)
        
        
        if (as.numeric(tools::file_path_sans_ext(packageVersion("GSVA"))) >= 1.5) {
          if (scoreMethod == "ssgsea") {
            ssGSEA_param <- GSVA::ssgseaParam(A,GS)
          } else if (scoreMethod == "gsva") {
            ssGSEA_param <- GSVA::gsvaParam(A,GS)
          } else if (scoreMethod == "plage") {
            ssGSEA_param <- GSVA::plageParam(A,GS)
          } else if (scoreMethod == "zscore") {
            ssGSEA_param <- GSVA::zscoreParam(A,GS)
          }
          ssgsea <- GSVA::gsva(ssGSEA_param)
        } else {
          ssgsea <- gsva(A,GS,method = scoreMethod, verbose = FALSE)
        }
        
        ssgsea2 = t(ssgsea)
        ssgsea3 = apply(ssgsea2, 2, scale);
        ssgsea4 = apply(ssgsea3, 1, rev)
        colnames(ssgsea4) = rownames(ssgsea2)
        
        neworder_gs <- rownames(ssgsea4)
        final_gs <- intersect(geneset_names,neworder_gs)
        
        ssgsea5 <- ssgsea4[final_gs,]
        
        meta2 <- meta[,c(colnames(meta)[1],metacol)]
        meta2[,2] <- as.factor(meta2[,2])
        meta2 <- meta2[order(meta2[,2]),]
        rownames(meta2) <- meta2[,1]
        meta2 <- meta2[,-1,drop = F]
        ssgsea5 <- ssgsea5[,rownames(meta2)]
        #ssgsea5 <- ssgsea5[,match(colnames(ssgsea5),samples_chosen)]
        ssgsea5 <- ssgsea5[,samples_chosen]
        ssgsea5
        
        
        
      })
      
      
      ssgseaheatmap_react <- reactive({
        
        row_names_choice <- input$ShowRowNames1SSheat
        col_names_choice <- input$ShowColNamesSSheat
        row_font <- input$heatmapFont1.r
        col_font <- input$heatmapFont1.c
        clust_cols_opt <- input$clustcolsSSheat
        clust_rows_opt <- input$clustrowsSSheat
        kmeans <- ifelse(is.na(input$NumClusters),1,input$NumClusters)
        if (is.null(input$ClusterMethodSSheat) == TRUE) {
          clust_method <- 'complete'
        } else if (is.null(input$ClusterMethodSSheat) == FALSE) {
          clust_method <- input$ClusterMethodSSheat
        }
        
        colAnno <- ssgsea_heat_anno()
        dataset <- ssgsea_heat_data()
        
        p <- suppressMessages(ComplexHeatmap::Heatmap(dataset,
                                                      top_annotation = colAnno,
                                                      clustering_method_rows = clust_method,
                                                      show_row_names = row_names_choice, show_column_names = col_names_choice,
                                                      cluster_rows = clust_rows_opt, cluster_columns = clust_cols_opt,
                                                      row_names_gp = gpar(fontsize = row_font), column_names_gp = gpar(fontsize = col_font),
                                                      heatmap_legend_param = list(title = "ssGSEA Score"),
                                                      border = F))
        draw(p, padding = unit(c(25, 50, 2, 2), "mm")) # unit(c(bottom,left,right,top))
        
      })
      
      ssgseaheatmap2_df_react <- reactive({
        
        meta <- meta_react()
        dataset <- ssgsea_heat2_data()
        if (isTruthy(input$ssGSEAGroupColMulti)) {
          meta <- meta[which(meta[,1] %in% colnames(dataset)),]
          meta_sub <- meta[,c(colnames(meta)[1],input$ssGSEAGroupColMulti), drop = F]
          dataset_t <- as.data.frame(t(dataset))
          dataset_t <- cbind(SampleName = rownames(dataset_t),dataset_t)
          heatdf <- merge(meta_sub,dataset_t, all.y = T, by.x = colnames(meta_sub)[1], by.y = "SampleName")
          heatdf <- heatdf[match(dataset_t[,1],heatdf[,1]),]
          heatdf
        } else {
          if (ncol(meta) == 2) {
            meta <- meta[which(meta[,1] %in% colnames(dataset)),]
            meta_sub <- meta
            dataset_t <- as.data.frame(t(dataset))
            dataset_t <- cbind(SampleName = rownames(dataset_t),dataset_t)
            heatdf <- merge(meta_sub,dataset_t, all.y = T, by.x = colnames(meta_sub)[1], by.y = "SampleName")
            heatdf <- heatdf[match(dataset_t[,1],heatdf[,1]),]
            heatdf
          } else {
            heatdf <- as.data.frame(t(dataset))
            heatdf <- cbind(SampleName = rownames(heatdf),heatdf)
            heatdf
          }
        }
        
      })
      
      output$ssgseaheatmap2_df <- DT::renderDataTable({
        
        df <- ssgseaheatmap2_df_react()
        dataset <- ssgsea_heat2_data()
        DT::datatable(df, 
                      options = list(lengthMenu = c(5,10,20,50,100,1000),
                                     pageLength = 10,
                                     scrollX = TRUE),
                      rownames = FALSE)%>%
          formatRound(columns = rownames(dataset), digits = 4)
        
      })
      
      ssgsea_heat2_data <- reactive({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        A <- A_raw()
        scoreMethod <- input$ssGSEAtype
        
        samples_chosen <- input$userheatsampSS
        
        #A <- A[,samples_chosen]
        A <- A[,which(colnames(A) %in% samples_chosen)]
        meta <- meta[which(meta[,1] %in% samples_chosen),]
        GS <- gs_react_heat()
        geneset_names <- names(GS)
        if (as.numeric(tools::file_path_sans_ext(packageVersion("GSVA"))) >= 1.5) {
          if (scoreMethod == "ssgsea") {
            ssGSEA_param <- GSVA::ssgseaParam(A,GS)
          } else if (scoreMethod == "gsva") {
            ssGSEA_param <- GSVA::gsvaParam(A,GS)
          } else if (scoreMethod == "plage") {
            ssGSEA_param <- GSVA::plageParam(A,GS)
          } else if (scoreMethod == "zscore") {
            ssGSEA_param <- GSVA::zscoreParam(A,GS)
          }
          ssgsea <- GSVA::gsva(ssGSEA_param)
        } else {
          ssgsea <- gsva(A,GS,method = scoreMethod, verbose = FALSE, ssgsea.norm = F)
        }
        
        
        SD=apply(ssgsea,1, sd, na.rm = TRUE) #get SD
        
        ssgsea2 = t(ssgsea)
        ssgsea3 = apply(ssgsea2, 2, scale);
        ssgsea4 = apply(ssgsea3, 1, rev)
        colnames(ssgsea4) = rownames(ssgsea2)
        
        ssgsea5 = ssgsea4 * SD #multiply zscore matrix by SD
        
        neworder_gs <- rownames(ssgsea5)
        final_gs <- intersect(geneset_names,neworder_gs)
        
        ssgsea5 <- ssgsea5[final_gs,]
        meta2 <- meta[,c(colnames(meta)[1],metacol)]
        meta2[,2] <- as.factor(meta2[,2])
        meta2 <- meta2[order(meta2[,2]),]
        rownames(meta2) <- meta2[,1]
        meta2 <- meta2[,-1,drop = F]
        ssgsea5 <- ssgsea5[,rownames(meta2)]
        #ssgsea5 <- ssgsea5[,match(colnames(ssgsea5),samples_chosen)]
        ssgsea5 <- ssgsea5[,samples_chosen]
        ssgsea5
        
      })
      
      ssgsea_heat_anno2 <- reactive({
        
        if (isTruthy(input$ssGSEAGroupColMulti)) {
          meta <- meta_react()
          dataset <- ssgsea_heat2_data()
          meta <- meta[which(meta[,1] %in% colnames(dataset)),]
          meta <- meta[match(colnames(dataset),meta[,1]),]
          rownames(meta) <- meta[,1]
          anno_cols <- input$ssGSEAGroupColMulti
          meta_sub <- meta[,anno_cols, drop = F]
          meta_sub[,anno_cols] <- as.data.frame(lapply(meta_sub[,anno_cols, drop = F], function(x) {
            if (!is.numeric(x) | all(suppressWarnings(as.numeric(x))%%1==0)) {
              return(factor(x))
            } else {
              return(x)
            }
          }))
          #meta_sub[,input$input$ssGSEAGroupColMulti] <- as.data.frame(lapply(meta_sub[,input$input$ssGSEAGroupColMulti, drop = F], factor))
          colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                      name = anno_cols,
                                                      which = 'col'
          )
          colAnn
        } else {
          meta <- meta_react()
          if (ncol(meta) == 2) {
            dataset <- ssgsea_heat2_data()
            meta <- meta[which(meta[,1] %in% colnames(dataset)),]
            meta <- meta[match(colnames(dataset),meta[,1]),]
            rownames(meta) <- meta[,1]
            meta_sub <- meta[,2, drop = F]
            meta_sub[,1] <- as.factor(meta_sub[,1])
            colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                        name = colnames(meta_sub)[2],
                                                        which = 'col'
            )
            colAnn
          } else {
            colAnn <- NULL
            colAnn
          }
        }
        
      })
      
      
      ssgseaheatmap2_react <- reactive({
        
        row_names_choice <- input$ShowRowNames1SSheat
        col_names_choice <- input$ShowColNamesSSheat
        row_font <- input$heatmapFont1.r
        col_font <- input$heatmapFont1.c
        clust_cols_opt <- input$clustcolsSSheat
        clust_rows_opt <- input$clustrowsSSheat
        kmeans <- ifelse(is.na(input$NumClusters),1,input$NumClusters)
        if (is.null(input$ClusterMethodSSheat) == TRUE) {
          clust_method <- 'complete'
        }
        else if (is.null(input$ClusterMethodSSheat) == FALSE) {
          clust_method <- input$ClusterMethodSSheat
        }
        
        dataset <- ssgsea_heat2_data()
        colAnno <- ssgsea_heat_anno2()
        

        p <- suppressMessages(ComplexHeatmap::Heatmap(dataset,
                                                      top_annotation = colAnno,
                                                      clustering_method_rows = clust_method,
                                                      show_row_names = row_names_choice, show_column_names = col_names_choice,
                                                      cluster_rows = clust_rows_opt, cluster_columns = clust_cols_opt,
                                                      row_names_gp = gpar(fontsize = row_font), column_names_gp = gpar(fontsize = col_font),
                                                      heatmap_legend_param = list(title = "ssGSEA Score"),
                                                      border = F))
        draw(p, padding = unit(c(25, 50, 2, 2), "mm")) # unit(c(bottom,left,right,top))
        
        
        #meta <- meta_react()
        #metacol <- metacol_reactVal()
        #A <- expr_mat_react()
        ##if (ncol(meta) > 2) {
        ##  metacol <- input$GSEAmetaCol
        ##}
        ##else if (ncol(meta) == 2) {
        ##  metacol <- colnames(meta)[2]
        ##}
        #if (is.null(input$ClusterMethodSSheat) == TRUE) {
        #  clust_method <- 'complete'
        #}
        #else if (is.null(input$ClusterMethodSSheat) == FALSE) {
        #  clust_method <- input$ClusterMethodSSheat
        #}
        #color_choice <- input$ColorPalette_gseaHeat
        #clust_cols_opt <- input$clustcolsSSheat
        #clust_rows_opt <- input$clustrowsSSheat
        #scoreMethod <- input$ssGSEAtype
        #
        #if (is.null(input$ssgseaHeatGS) == TRUE) {
        #  #geneset_names <- gs_names_start
        #  geneset_names <- ssGSEA_Heat_GS()
        #}
        #else if (is.null(input$ssgseaHeatGS) == FALSE) {
        #  geneset_names <- input$ssgseaHeatGS
        #}
        #samples_chosen <- input$userheatsampSS
        #
        #A <- A[,samples_chosen]
        #meta <- meta[which(meta[,1] %in% samples_chosen),]
        #
        #if (input$tables == 1) {
        #  GS <- gs[geneset_names]
        #}
        #if (input$tables == 3) {
        #  GS <- gs2[geneset_names]
        #}
        #if (input$tables == 5) {
        #  GS <- RDataListGen()[geneset_names]
        #}
        #ssgsea <- gsva(A, GS, method = scoreMethod, verbose = F, ssgsea.norm = F)
        #
        #SD=apply(ssgsea,1, sd, na.rm = TRUE) #get SD
        #
        #ssgsea2 = t(ssgsea)
        #ssgsea3 = apply(ssgsea2, 2, scale);
        #ssgsea4 = apply(ssgsea3, 1, rev)
        #colnames(ssgsea4) = rownames(ssgsea2)
        #
        #ssgsea5 = ssgsea4 * SD #multiply zscore matrix by SD
        #
        #neworder_gs <- rownames(ssgsea5)
        #final_gs <- intersect(geneset_names,neworder_gs)
        #
        #ssgsea5 <- ssgsea5[final_gs,]
        #
        #minimum = min(ssgsea5)
        #maximum = max(ssgsea5)
        #if (abs(min(ssgsea5)) > abs(max(ssgsea5))) {
        #  ssgsea5[ssgsea5 < -abs(max(ssgsea5))] = -abs(max(ssgsea5))
        #} else {
        #  ssgsea5[ssgsea5 > abs(min(ssgsea5))] = abs(min(ssgsea5))
        #}
        #bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
        ##Heatmap color
        #col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
        #if (color_choice == "original") {
        #  HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
        #  hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
        #}
        #else if (color_choice %in% col_sets) {
        #  HeatMap_Colors <- brewer.pal(n = 5, color_choice)
        #  hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
        #}
        #else if (color_choice == "Inferno") {
        #  hmcols <- inferno(500)
        #}
        #else if (color_choice == "Viridis") {
        #  hmcols <- viridis(500)
        #}
        #else if (color_choice == "Plasma") {
        #  hmcols <- plasma(500)
        #}
        #else if (color_choice == "OmniBlueRed") {
        #  hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
        #}
        #else if (color_choice == "LightBlueBlackRed") {
        #  hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
        #}
        #else if (color_choice == "GreenBlackRed") {
        #  hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
        #}
        #
        #meta2 <- meta[,c(colnames(meta)[1],metacol)]
        #meta2[,2] <- as.factor(meta2[,2])
        #meta2 <- meta2[order(meta2[,2]),]
        #rownames(meta2) <- meta2[,1]
        #meta2 <- meta2[,-1,drop = F]
        #ssgsea5 <- ssgsea5[,rownames(meta2)]
        #
        ##meta <- meta[order(meta[,2]),]
        ##type <- meta[,2]
        ##meta2 <- as.data.frame(type)
        ##rownames(meta2) <- meta[,1]
        #
        #hm <- pheatmap::pheatmap(ssgsea5,
        #                         cluster_col = clust_cols_opt,
        #                         cluster_row = clust_rows_opt,
        #                         fontsize_row = row_font,
        #                         fontsize_col = col_font,
        #                         show_rownames = row_names_choice,
        #                         show_colnames = col_names_choice,
        #                         annotation_col = meta2,
        #                         clustering_method = clust_method,
        #                         color = hmcols)
        #
        #hm
        
        
      })
      
      ssgseaheatmap3_df_react <- reactive({
        
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        AvgssGSEADF <- ssgseaheatmap3_react_df()
        anno_meta <- data.frame(GroupName = colnames(AvgssGSEADF))
        anno_meta[,metacol] <- gsub("Avg_ssGSEA_","",anno_meta[,1])
        
        dataset_t <- as.data.frame(t(AvgssGSEADF))
        dataset_t <- cbind(GroupName = rownames(dataset_t),dataset_t)
        heatdf <- merge(anno_meta,dataset_t, all.y = T, by = "GroupName", sort = F)
        heatdf <- heatdf[match(dataset_t[,1],heatdf[,1]),]
        heatdf
        
      })
      
      output$ssgseaheatmap3_df <- DT::renderDataTable({
        
        df <- ssgseaheatmap3_df_react()
        dataset <- ssgseaheatmap3_react_df()
        DT::datatable(df, 
                      options = list(lengthMenu = c(5,10,20,50,100,1000),
                                     pageLength = 10,
                                     scrollX = TRUE),
                      rownames = FALSE)%>%
          formatRound(columns = rownames(dataset), digits = 4)
        
      })
      
      ssgseaheatmap3_react_df <- reactive({
        
        row_names_choice <- input$ShowRowNames1SSheat
        col_names_choice <- input$ShowColNamesSSheat
        row_font <- input$heatmapFont1.r
        col_font <- input$heatmapFont1.c
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        A <- A_raw()
        if (is.null(input$ClusterMethodSSheat) == TRUE) {
          clust_method <- 'complete'
        }
        else if (is.null(input$ClusterMethodSSheat) == FALSE) {
          clust_method <- input$ClusterMethodSSheat
        }
        color_choice <- input$ColorPalette_gseaHeat
        clust_cols_opt <- input$clustcolsSSheat
        clust_rows_opt <- input$clustrowsSSheat
        scoreMethod <- input$ssGSEAtype
        
        samples_chosen <- input$userheatsampSS
        
        #A <- A[,samples_chosen]
        A <- A[,which(colnames(A) %in% samples_chosen)]
        meta <- meta[which(meta[,1] %in% samples_chosen),]
        #GS <- gs_react()
        GS <- gs_react_heat()
        geneset_names <- names(GS)
        
        
        if (as.numeric(tools::file_path_sans_ext(packageVersion("GSVA"))) >= 1.5) {
          if (scoreMethod == "ssgsea") {
            ssGSEA_param <- GSVA::ssgseaParam(A,GS)
          } else if (scoreMethod == "gsva") {
            ssGSEA_param <- GSVA::gsvaParam(A,GS)
          } else if (scoreMethod == "plage") {
            ssGSEA_param <- GSVA::plageParam(A,GS)
          } else if (scoreMethod == "zscore") {
            ssGSEA_param <- GSVA::zscoreParam(A,GS)
          }
          ssgsea <- GSVA::gsva(ssGSEA_param)
        } else {
          ssgsea <- gsva(A,GS,method = scoreMethod, verbose = FALSE, ssgsea.norm = F)
        }
        
        
        SD=apply(ssgsea,1, sd, na.rm = TRUE) #get SD
        
        ssgsea2 = t(ssgsea)
        ssgsea3 = apply(ssgsea2, 2, scale);
        ssgsea4 = apply(ssgsea3, 1, rev)
        colnames(ssgsea4) = rownames(ssgsea2)
        
        ssgsea5 = ssgsea4 * SD #multiply zscore matrix by SD
        
        group_choices <- unique(meta[,metacol])
        AvgssGSEADF <- data.frame(rownames(ssgsea5))
        
        
        for (i in group_choices) {
          samples <- meta[which(meta[,metacol] == i),1]
          if (length(samples) <= 1) {
            AvgssGSEADF[,paste("Avg_ssGSEA_",i, sep = "")] <- ssgsea5[,samples]
          }
          else if (length(samples) > 1) {
            AvgssGSEADF[,paste("Avg_ssGSEA_",i, sep = "")] <- rowMeans(ssgsea5[,samples])
          }
        }
        
        rownames(AvgssGSEADF) <- AvgssGSEADF[,1]
        AvgssGSEADF <- AvgssGSEADF[,-1]
        
        neworder_gs <- rownames(AvgssGSEADF)
        final_gs <- intersect(geneset_names,neworder_gs)
        
        AvgssGSEADF <- AvgssGSEADF[final_gs,]
        
        minimum = min(AvgssGSEADF, na.rm = T)
        maximum = max(AvgssGSEADF, na.rm = T)
        if (abs(minimum) > abs(maximum)) {
          AvgssGSEADF[AvgssGSEADF < -abs(maximum)] = -abs(maximum)
        } else {
          AvgssGSEADF[AvgssGSEADF > abs(minimum)] = abs(minimum)
        }
        AvgssGSEADF
        
      })
      
      ssgseaheatmap3_react <- reactive({
        
        row_names_choice <- input$ShowRowNames1SSheat
        col_names_choice <- input$ShowColNamesSSheat
        row_font <- input$heatmapFont1.r
        col_font <- input$heatmapFont1.c
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        A <- A_raw()
        if (is.null(input$ClusterMethodSSheat) == TRUE) {
          clust_method <- 'complete'
        }
        else if (is.null(input$ClusterMethodSSheat) == FALSE) {
          clust_method <- input$ClusterMethodSSheat
        }
        color_choice <- input$ColorPalette_gseaHeat
        clust_cols_opt <- input$clustcolsSSheat
        clust_rows_opt <- input$clustrowsSSheat
        scoreMethod <- input$ssGSEAtype
        
        #samples_chosen <- input$userheatsampSS
        #
        ##A <- A[,samples_chosen]
        #A <- A[,which(colnames(A) %in% samples_chosen)]
        #meta <- meta[which(meta[,1] %in% samples_chosen),]
        ##GS <- gs_react()
        #GS <- gs_react_heat()
        #geneset_names <- names(GS)
        #
        #
        #if (as.numeric(tools::file_path_sans_ext(packageVersion("GSVA"))) >= 1.5) {
        #  if (scoreMethod == "ssgsea") {
        #    ssGSEA_param <- GSVA::ssgseaParam(A,GS)
        #  } else if (scoreMethod == "gsva") {
        #    ssGSEA_param <- GSVA::gsvaParam(A,GS)
        #  } else if (scoreMethod == "plage") {
        #    ssGSEA_param <- GSVA::plageParam(A,GS)
        #  } else if (scoreMethod == "zscore") {
        #    ssGSEA_param <- GSVA::zscoreParam(A,GS)
        #  }
        #  ssgsea <- GSVA::gsva(ssGSEA_param)
        #} else {
        #  ssgsea <- gsva(A,GS,method = scoreMethod, verbose = FALSE, ssgsea.norm = F)
        #}
        #
        #
        #SD=apply(ssgsea,1, sd, na.rm = TRUE) #get SD
        #
        #ssgsea2 = t(ssgsea)
        #ssgsea3 = apply(ssgsea2, 2, scale);
        #ssgsea4 = apply(ssgsea3, 1, rev)
        #colnames(ssgsea4) = rownames(ssgsea2)
        #
        #ssgsea5 = ssgsea4 * SD #multiply zscore matrix by SD
        #
        #group_choices <- unique(meta[,metacol])
        #AvgssGSEADF <- data.frame(rownames(ssgsea5))
        #
        #
        #for (i in group_choices) {
        #  samples <- meta[which(meta[,metacol] == i),1]
        #  if (length(samples) <= 1) {
        #    AvgssGSEADF[,paste("Avg_ssGSEA_",i, sep = "")] <- ssgsea5[,samples]
        #  }
        #  else if (length(samples) > 1) {
        #    AvgssGSEADF[,paste("Avg_ssGSEA_",i, sep = "")] <- rowMeans(ssgsea5[,samples])
        #  }
        #}
        #
        #rownames(AvgssGSEADF) <- AvgssGSEADF[,1]
        #AvgssGSEADF <- AvgssGSEADF[,-1]
        #
        #neworder_gs <- rownames(AvgssGSEADF)
        #final_gs <- intersect(geneset_names,neworder_gs)
        #
        #AvgssGSEADF <- AvgssGSEADF[final_gs,]
        #
        AvgssGSEADF <- ssgseaheatmap3_react_df()
        minimum = min(AvgssGSEADF, na.rm = T)
        maximum = max(AvgssGSEADF, na.rm = T)
        #if (abs(minimum) > abs(maximum)) {
        #  AvgssGSEADF[AvgssGSEADF < -abs(maximum)] = -abs(maximum)
        #} else {
        #  AvgssGSEADF[AvgssGSEADF > abs(minimum)] = abs(minimum)
        #}
        
        
        bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
        #Heatmap color
        col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
        if (color_choice == "original") {
          HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
          hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
        }
        else if (color_choice %in% col_sets) {
          HeatMap_Colors <- brewer.pal(n = 5, color_choice)
          hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
        }
        else if (color_choice == "Inferno") {
          hmcols <- inferno(500)
        }
        else if (color_choice == "Viridis") {
          hmcols <- viridis(500)
        }
        else if (color_choice == "Plasma") {
          hmcols <- plasma(500)
        }
        else if (color_choice == "OmniBlueRed") {
          hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
        }
        else if (color_choice == "LightBlueBlackRed") {
          hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
        }
        else if (color_choice == "GreenBlackRed") {
          hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
        }
        
        anno_meta <- data.frame(colnames(AvgssGSEADF))
        anno_meta[,metacol] <- gsub("Avg_ssGSEA_","",anno_meta[,1])
        rownames(anno_meta) <- anno_meta[,1]
        anno_meta <- anno_meta[,-1,drop = F]
        
        hm <- pheatmap::pheatmap(as.matrix(AvgssGSEADF),
                                 cluster_col = clust_cols_opt,
                                 cluster_row = clust_rows_opt,
                                 fontsize_row = row_font,
                                 fontsize_col = col_font,
                                 show_rownames = row_names_choice,
                                 show_colnames = col_names_choice,
                                 angle_col = 90,
                                 annotation_col = anno_meta,
                                 clustering_method = clust_method,
                                 color = hmcols)
        
        hm
        
        
      })
      
      ####----Gene Scatter----####
      
      #render gene expression comparison scatter plot
      output$geneScatter0 <- renderPlotly({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        expr <- expr_react()
        title_font <- input$GeneScatterTitleSize
        axis_font <- input$GeneScatterAxisSize
        #if (ncol(meta) > 2){
        #  req(input$ScatterPlotMetaCol)
        #  metacol <- input$ScatterPlotMetaCol
        #}
        #else if (ncol(meta) == 2){
        #  metacol <- colnames(meta)[2]
        #}
        
        
        #log if user designates
        if (input$logask == TRUE) {
          expr <- log2(expr + 1)
        }
        #transpose
        expr_t <- as.data.frame(t(expr))
        #reorder rowname to match meta for merging
        samporder <- meta[,1]
        expr_t2 <- as.data.frame(expr_t[samporder,])
        #add type
        expr_t3 <- expr_t2 %>%
          mutate(type = case_when(
            rownames(expr_t2) == meta[,1] ~ meta[,metacol],
          ))
        expr_t3 <- expr_t3 %>%
          relocate(type)
        expr_t3$type <- as.factor(expr_t3$type)
        colnames(expr_t3)[which(colnames(expr_t3) == "type")] <- metacol
        #user gene input
        gene1.u <- input$scatterG1
        gene2.u <- input$scatterG2
        #get columns and info based off user input
        gene1 <- expr_t3[,gene1.u]
        gene2 <- expr_t3[,gene2.u]
        if (input$logask == TRUE) {
          gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression (log2 +1)", sep = "")
          gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression (log2 +1)", sep = "")
        }
        else if (input$logask == FALSE) {
          gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression", sep = "")
          gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression", sep = "")
        }
        #plot
        p <- ggplot(expr_t3, aes(x = gene1, y = gene2,
                                 color = expr_t3[,metacol],
                                 text = paste("</br> Sample: ", rownames(expr_t3),
                                              "</br> ",metacol,": ", expr_t3[,metacol],
                                              sep = ""))) +
          geom_point() +
          theme_minimal() +
          labs(x = gene1.n, y = gene2.n,
               title = paste(gene1.n, " vs. ", gene2.n, sep = ''),
               color = metacol) +
          theme(axis.title = element_text(size = axis_font),
                plot.title = element_text(size = title_font))
        
        ggplotly(p, tooltip = 'text')
        
      })
      
      ####----Data Tables----####
      
      
      #render MSigDB gene set table
      output$msigdbTable <- DT::renderDataTable({
        DT::datatable(msigdb.gsea2(),
                      selection = list(mode = 'single', selected = 1),
                      options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T))
        
      })
      
      #create user input gene set table
      output$GStable.u <- DT::renderDataTable({
        gmt.u <- input$user.gmt.file
        ext <- tools::file_ext(gmt.u$datapath)
        req(gmt.u)
        validate(need(ext == c("gmt","tsv","txt"), "Please upload .gmt, .tsv, or .txt file"))
        if (ext == "gmt") {
          gmt.us <- read.gmt(gmt.u$datapath)
        }
        else {
          #gmt.us <- as.data.frame(read_delim(gmt.u$datapath, delim = '\t'))
          gmt.us <- as.data.frame(fread(gmt.u$datapath, sep = '\t'))
          #colnames(gmt.us) <- c("term","gene")
        }
        GeneSet <- as.data.frame(unique(gmt.us[1]))
        rownames(GeneSet) <- 1:nrow(GeneSet)
        colnames(GeneSet)[1] <- "Gene_Set"
        DT::datatable(GeneSet,
                      selection = list(mode = 'single', selected = 1),
                      options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100")))
      })
      
      #render tab2 gene set table
      output$tab2table <- DT::renderDataTable({
        DT::datatable(GeneSet2(),
                      selection = list(mode = 'single', selected = 1),
                      options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100")))
        
      })
      
      LeadingEdgeGenes_react <- reactive({
        req(datasetInput())
        if (length(input$GeneSetTable_rows_selected) > 0){
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          
          if (length(input$GeneSetTable_rows_selected) > 1) {
            GeneSymbol <- data.frame(GeneSet = character(0), GeneSymbol = character)
            for (row in seq(length(input$GeneSetTable_rows_selected))) {
              GS <- gsea.df[row,1]
              genes1 <- gsea.df[which(gsea.df$Description==GS),"core_enrichment"]
              genes2 <- sort(strsplit(genes1,"/")[[1]])
              term <- rep(GS, times = length(genes2))
              GeneSymbol <- rbind(GeneSymbol,
                                  data.frame(GeneSet = term,
                                             GeneSymbol = genes2))
            }
            colnames(GeneSymbol) <- c("Gene Set","Gene Symbol")
            GeneSymbol <- GeneSymbol[order(GeneSymbol[,2]),]
          } else {
            GS <- unique(gmt_react()[,1])
            genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
            genes2 <- sort(strsplit(genes1,"/")[[1]])
            GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
            colnames(GeneSymbol)[1] <- "Gene Symbol"
          }
          GeneSymbol
        }
        
      })
      
      #render leading edge genes list
      output$LeadingEdgeGenes <- DT::renderDataTable({
        GeneSymbol  <- LeadingEdgeGenes_react()
        DT::datatable(GeneSymbol, options = list(paging = F), rownames = F)
      })
      
      #render pre-loaded enriched signatures table
      output$enrich_sig_table <- DT::renderDataTable({
        gsea.df <- as_tibble(get(ES_Tab_List[which(ES_Tab_List == paste("ES_table",match(input$SigTableChoice, SigNames),sep = ""))]))
        DT::datatable(gsea.df,
                      extensions = c("KeyTable", "FixedHeader"),
                      caption = "Enriched Signatures",
                      options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"),scrollX = T)) %>%
          formatRound(columns = c(2:10), digits = 2)
      })
      
      #render user generated enriched signature table based off other gmt
      output$enrich_sig_table_gen <- DT::renderDataTable({
        #meta <- meta_react()
        #metacol <- metacol_reactVal()
        #if (length(colnames(meta)) == 2) {
        #  metacol <- colnames(meta)[2]
        #}
        #A <- A_raw()
        #gmt <- gmt_react()
        req(GeneratedMSigDBEST())
        gsea.df <- GeneratedMSigDBEST()
        gsea.df <- as_tibble(gsea.df)
        gsea.df <- gsea.df[,-2]
        ## displaying the GSEA results as interactive data table
        DT::datatable(gsea.df,
                      extensions = c("KeyTable", "FixedHeader"),
                      caption = "Enriched Signatures",
                      options = list(keys = T, searchHighlight = T, pageLength = 10,
                                     lengthMenu = c("10", "25", "50", "100"), scrollX = T)) %>%
          formatSignif(columns = c(5:7), digits = 3) %>%
          formatRound(columns = c(3,4), digits = 3)
        
        #if (input$tables == 3) {
        #  #if (ncol(meta) > 2) {
        #  #  req(input$GSEAmetaCol)
        #  #  metacol <- input$GSEAmetaCol
        #  #}
        #  #else if (ncol(meta) == 2) {
        #  #  metacol <- colnames(meta)[2]
        #  #}
        #  groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
        #  groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
        #  ##----Signal-to-Noise Calculation----##
        #  A <- A + 0.00000001
        #  P = as.matrix(as.numeric(colnames(A) %in% groupA))
        #  n1 <- sum(P[,1])
        #  M1 <- A %*% P
        #  M1 <- M1/n1
        #  A2 <- A*A
        #  S1 <- A2 %*% P
        #  S1 <- S1/n1 - M1*M1
        #  S1 <- sqrt(abs((n1/(n1-1)) * S1))
        #  P = as.matrix(as.numeric(colnames(A) %in% groupB))
        #  n2 <- sum(P[,1])
        #  M2 <- A %*% P
        #  M2 <- M2/n2
        #  A2 <- A*A
        #  S2 <- A2 %*% P
        #  S2 <- S2/n2 - M2*M2
        #  S2 <- sqrt(abs((n2/(n2-1)) * S2))
        #  rm(A2)
        #  # small sigma "fix" as used in GeneCluster
        #  S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
        #  S2 <- ifelse(S2 == 0, 0.2, S2)
        #  S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
        #  S1 <- ifelse(S1 == 0, 0.2, S1)
        #  M1 <- M1 - M2
        #  rm(M2)
        #  S1 <- S1 + S2
        #  rm(S2)
        #  s2n.matrix <- M1/S1
        #  ##----Reformatting----##
        #  s2n.df <- as.data.frame(s2n.matrix)
        #  s2n.df$GeneID <- rownames(s2n.df)
        #  rownames(s2n.df) <- NULL
        #  data <- dplyr::select(s2n.df, GeneID, V1)
        #  data.gsea <- data$V1
        #  names(data.gsea) <- as.character(data$GeneID)
        #  s2n.matrix.s <- sort(data.gsea, decreasing = T)
        #  ##----GSEA----##
        #  gmt.i <- tab2()
        #  gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
        #  gsea.df <- as_tibble(gsea.res@result)
        #  ## displaying the GSEA results as interactive data table
        #  DT::datatable(gsea.df,
        #                extensions = c("KeyTable", "FixedHeader"),
        #                caption = "Enriched Signatures",
        #                options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T)) %>%
        #    formatRound(columns = c(2:10), digits = 2)
        #}
        #else if (input$tables == 5) {
        #  #if (ncol(meta) > 2) {
        #  #  metacol <- input$GSEAmetaCol
        #  #}
        #  #else if (ncol(meta) == 2) {
        #  #  metacol <- colnames(meta)[2]
        #  #}
        #  groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
        #  groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
        #  ##----Signal-to-Noise Calculation----##
        #  A <- A + 0.00000001
        #  P = as.matrix(as.numeric(colnames(A) %in% groupA))
        #  n1 <- sum(P[,1])
        #  M1 <- A %*% P
        #  M1 <- M1/n1
        #  A2 <- A*A
        #  S1 <- A2 %*% P
        #  S1 <- S1/n1 - M1*M1
        #  S1 <- sqrt(abs((n1/(n1-1)) * S1))
        #  P = as.matrix(as.numeric(colnames(A) %in% groupB))
        #  n2 <- sum(P[,1])
        #  M2 <- A %*% P
        #  M2 <- M2/n2
        #  A2 <- A*A
        #  S2 <- A2 %*% P
        #  S2 <- S2/n2 - M2*M2
        #  S2 <- sqrt(abs((n2/(n2-1)) * S2))
        #  rm(A2)
        #  # small sigma "fix" as used in GeneCluster
        #  S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
        #  S2 <- ifelse(S2 == 0, 0.2, S2)
        #  S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
        #  S1 <- ifelse(S1 == 0, 0.2, S1)
        #  M1 <- M1 - M2
        #  rm(M2)
        #  S1 <- S1 + S2
        #  rm(S2)
        #  s2n.matrix <- M1/S1
        #  ##----Reformatting----##
        #  s2n.df <- as.data.frame(s2n.matrix)
        #  s2n.df$GeneID <- rownames(s2n.df)
        #  rownames(s2n.df) <- NULL
        #  data <- dplyr::select(s2n.df, GeneID, V1)
        #  data.gsea <- data$V1
        #  names(data.gsea) <- as.character(data$GeneID)
        #  s2n.matrix.s <- sort(data.gsea, decreasing = T)
        #  ##----GSEA----##
        #  gmt.i <- GStable.ubg()
        #  gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
        #  gsea.df <- as_tibble(gsea.res@result)
        #  ## displaying the GSEA results as interactive data table
        #  DT::datatable(gsea.df,
        #                extensions = c("KeyTable", "FixedHeader"),
        #                caption = "Enriched Signatures",
        #                options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T)) %>%
        #    formatRound(columns = c(2:10), digits = 2)
        #}
        #else if (input$tables == 1) {
        #  if (input$GenerateEST == TRUE) {
        #    #extract results and convert to tibble
        #    gsea.df <- GeneratedMSigDBEST()
        #    gsea.df <- as_tibble(gsea.df)
        #    ## displaying the GSEA results as interactive data table
        #    DT::datatable(gsea.df,
        #                  extensions = c("KeyTable", "FixedHeader"),
        #                  caption = "Enriched Signatures",
        #                  options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T)) %>%
        #      formatRound(columns = c(2:10), digits = 2)
        #  }
        #}
      })
      
      DEGtable1_react <- reactive({
        
        req(meta_react())
        req(metacol_reactVal())
        req(expr_react())
        req(input$volcanoCompChoice2)
        if (input$volcanoCompChoice2 == "DESeq2") {
          if (isTruthy(metacol_reactVal())) {
            desCol <- metacol_reactVal()
            expr <- expr_raw()
            meta <- meta_react()
            covars <- input$DEGCovarSelectTab
            meta <- meta[,c(colnames(meta)[1],desCol,covars)]
            meta <- meta[complete.cases(meta),]
            rownames(expr) <- expr[,1]
            expr <- expr[,-1]
            mat <- as.matrix(expr)
            sortedIndices <- match(colnames(mat), meta[,1])
            meta_deseq <- meta[sortedIndices,,drop=FALSE]
            rownames(meta_deseq) <- meta_deseq[,1]
            if (desCol != "Select column") {
              desRef <- input$DESeqDesignColRef_tab
              if (isTruthy(desRef)) {
                meta_deseq[,desCol] <- relevel(factor(meta_deseq[,desCol]), ref = desRef)
                colnames(meta_deseq) <- gsub(" ","_",colnames(meta_deseq))
                colnames(meta_deseq) <- gsub("[[:punct:]]","_",colnames(meta_deseq))
                desCol <- gsub(" ","_",desCol)
                desCol <- gsub("[[:punct:]]","_",desCol)
                covars <- gsub(" ","_",covars)
                covars <- gsub("[[:punct:]]","_",covars)
                
                if (isTruthy(covars)) {
                  meta_deseq[,covars] <- lapply(meta_deseq[,covars,drop = F], factor)
                  desCol <- paste(c(covars, desCol),collapse = "+")
                }
                dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat,
                                                      colData = meta_deseq,
                                                      design = as.formula(paste0("~ ",desCol)))
                dds <- DESeq(dds)
                res <- results(dds)
                resOrdered <- res[order(res$padj),]
                top1 <- as.data.frame(resOrdered)
                top1
              }
            }
          }
        } else {
          #req(input$comparisonA2.DEG)
          #req(input$comparisonB2.DEG)
          meta <- meta_react()
          metacol <- metacol_reactVal()
          if (length(colnames(meta)) == 2) {
            metacol <- colnames(meta)[2]
          }
          expr <- expr_react()
          covars <- input$DEGCovarSelectTab
          
          metaSub <- meta[,c(colnames(meta)[1],metacol,covars)]
          metaSub <- metaSub[complete.cases(metaSub),]
          
          
          if (input$volcanoCompChoice2 == "Limma: Two groups") {
            groupA <- input$comparisonA2.DEG
            groupB <- input$comparisonB2.DEG
            A <- metaSub[which(metaSub[,metacol] == groupA),1]
            B <- metaSub[which(metaSub[,metacol] == groupB),1]
          } else if (input$volcanoCompChoice2 == "Limma: One group") {
            groupA <- input$comparisonA2.DEG_one
            groupB <- paste0("Not_",groupA)
            metaSub[,paste0(metacol,"_Dichot")] <- NA
            metaSub[which(metaSub[,metacol] == groupA),paste0(metacol,"_Dichot")] <- groupA
            metaSub[which(metaSub[,metacol] != groupA),paste0(metacol,"_Dichot")] <- paste0("Not_",groupA)
            A <- metaSub[which(metaSub[,paste0(metacol,"_Dichot")] == groupA),1]
            B <- metaSub[which(metaSub[,paste0(metacol,"_Dichot")] == paste0("Not_",groupA)),1]
            metacol <- paste0(metacol,"_Dichot")
          }
          
          
          metaSub <- metaSub[which(metaSub[,1] %in% c(A,B)),]
          colnames(metaSub) <- gsub(" ","_",colnames(metaSub))
          colnames(metaSub) <- gsub("[[:punct:]]","_",colnames(metaSub))
          metacol <- gsub(" ","_",metacol)
          metacol <- gsub("[[:punct:]]","_",metacol)
          covars <- gsub(" ","_",covars)
          covars <- gsub("[[:punct:]]","_",covars)
          metaSub[,metacol] <- factor(metaSub[,metacol], levels = c(groupA,groupB))
          
          mat <- expr[,metaSub[,1]]
          mat <- log2(mat + 1.0)
          if (isTruthy(covars)) {
            metaSub[,covars] <- lapply(metaSub[,covars,drop = F], factor)
            form <- as.formula(paste0("~0 +",paste(c(metacol,covars),collapse = "+")))
          } else {
            #metaSub[,metacol] <- factor(metaSub[,metacol])
            form <- as.formula(paste0("~0 +",metacol))
          }
          
          if (nlevels(metaSub[,metacol]) >= 2) {
            designA <- eval(substitute(model.matrix(form, data = metaSub)))
            
            fit <- lmFit(mat, design = designA)
            colnames(designA) <- gsub(" ","_",colnames(designA))
            colnames(designA) <- gsub("[[:punct:]]","_",colnames(designA))
            contrast.matrix <- makeContrasts(contrasts = paste0(colnames(designA)[1],"-",colnames(designA)[2]), levels = designA)
            fit2 <- contrasts.fit(fit, contrast.matrix)
            fit2 <- eBayes(fit2)
            options(digits = 4)
            top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
            top1
          }
          
        }
        
        
      })
      
      #render DEG table
      output$DEGtable1 <- DT::renderDataTable({
        
        top1 <- DEGtable1_react()
        DT::datatable(top1, options = list(lengthMenu = c(50,100,1000, 5000, 10000), pageLength = 100, scrollX = TRUE),
                      selection=list(mode = "multiple"))
        
        
        
      })
      
      #render up regulated pathway enrichment data table
      output$UpRegPathwayTable1 <- DT::renderDataTable({
        adjp <- input$pathpval
        FC <- input$pathFC
        if (ncol(meta) > 2) {
          metacol <- input$EnrichRPathMetaCol
        }
        else if (ncol(meta) == 2) {
          metacol <- colnames(meta)[2]
        }
        A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
        B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
        mat <- expr[,c(A,B)]
        mat <- log2(mat + 1.0)
        groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
        designA <- model.matrix(~0 + groupAOther)
        fit <- lmFit(mat, design = designA)
        contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        options(digits = 4)
        top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
        genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
        dbs <- listEnrichrDbs()
        enrichRLive <- TRUE
        if (is.null(dbs)) {
          enrichRLive <- FALSE
        }
        dbs <- input$SelectedPathway
        enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
        DT::datatable(enriched[[1]], options = list(lengthMenu = c(10,20,50,100,1000), pageLength = 10, scrollX = TRUE),
                      selection=list(mode = "multiple"))
      })
      
      #render down regulated pathway enrichment data table
      output$DnRegPathwayTable1 <- DT::renderDataTable({
        adjp <- input$pathpval
        FC <- input$pathFC
        if (ncol(meta) > 2) {
          metacol <- input$EnrichRPathMetaCol
        }
        else if (ncol(meta) == 2) {
          metacol <- colnames(meta)[2]
        }
        A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
        B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
        mat <- expr[,c(A,B)]
        mat <- log2(mat + 1.0)
        groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
        designA <- model.matrix(~0 + groupAOther)
        fit <- lmFit(mat, design = designA)
        contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        options(digits = 4)
        top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
        genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
        dbs <- listEnrichrDbs()
        enrichRLive <- TRUE
        if (is.null(dbs)) {
          enrichRLive <- FALSE
        }
        dbs <- input$SelectedPathway
        enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
        DT::datatable(enriched[[1]], options = list(lengthMenu = c(10,20,50,100,1000), pageLength = 10, scrollX = TRUE),
                      selection=list(mode = "multiple"))
      })
      
      #render gene list table for boxplot selection
      output$GeneListTable <- DT::renderDataTable({
        geneList <- geneList_raw()
        DT::datatable(geneList,
                      selection = list(mode = 'single', selected = 1),
                      options = list(keys = TRUE,
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = c("10", "25", "50", "100")))
      })
      
      #render gene list table for boxplot selection
      output$GeneListTable2 <- DT::renderDataTable({
        geneList <- geneList_raw()
        DT::datatable(geneList,
                      selection = list(mode = 'single', selected = 1),
                      options = list(keys = TRUE,
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = c("10", "25", "50", "100")),
                      rownames = F)
      })
      
      #render gene list table for boxplot selection
      output$GeneListTableBarPlot <- DT::renderDataTable({
        geneList <- geneList_raw()
        DT::datatable(geneList,
                      selection = list(mode = 'single', selected = 1),
                      options = list(keys = TRUE,
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = c("10", "25", "50", "100")),
                      rownames = F)
      })
      
      #render gene scatter plot data table
      output$geneScatterTable <- DT::renderDataTable({
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        expr <- expr_react()
        #log if user designates
        if (input$logask == TRUE) {
          expr <- log2(expr + 1)
        }
        #transpose
        expr_t <- as.data.frame(t(expr))
        #if (ncol(meta) > 2) {
        #  metacol <- input$ScatterPlotMetaCol
        #}
        #else if (ncol(meta) == 2) {
        #  metacol <- colnames(meta)[2]
        #}
        #reorder rowname to match meta for merging
        samporder <- meta[,1]
        expr_t2 <- as.data.frame(expr_t[samporder,])
        #add type
        expr_t3 <- expr_t2 %>%
          mutate(type = case_when(
            rownames(expr_t2) == meta[,1] ~ meta[,metacol],
          ))
        expr_t3 <- expr_t3 %>%
          relocate(type)
        colnames(expr_t3)[which(colnames(expr_t3) == "type")] <- metacol
        #user gene input
        gene1.u <- input$scatterG1
        gene2.u <- input$scatterG2
        #get columns and info based off user input
        gene1 <- expr_t3[,gene1.u]
        gene2 <- expr_t3[,gene2.u]
        if (input$logask == TRUE) {
          gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression (log2 +1)", sep = "")
          gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression (log2 +1)", sep = "")
        }
        else if (input$logask == FALSE) {
          gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression", sep = "")
          gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression", sep = "")
        }
        #make table
        Sample <- rownames(expr_t3)
        Type <- expr_t3[,metacol]
        gene1col <- expr_t3[,gene1.u]
        gene2col <- expr_t3[,gene2.u]
        scatterTab <- data.frame(Sample, Type, gene1col, gene2col)
        colnames(scatterTab)[c(2,3,4)] <- c(metacol,gene1.n, gene2.n)
        #table output,
        DT::datatable(scatterTab,
                      options = list(keys = TRUE,
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = c("10", "25", "50", "100")))
        
      })
      
      AvgGeneScatterTable_react <- reactive({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        expr <- expr_react()
        A_choice <- input$comparisonA2.avg
        B_choice <- input$comparisonB2.avg
        log_choice <- input$AvgExpLogFC
        # Get A and B sample names
        A <- meta[which(meta[,metacol] == A_choice),1]
        B <- meta[which(meta[,metacol] == B_choice),1]
        
        # Make A and B expression Matrices
        mat_A <- expr[,A]
        mat_B <- expr[,B]
        
        logFC_text <- ""
        
        # Log if user designates
        if (log_choice == TRUE) {
          mat_A <- log2(mat_A + 1)
          mat_B <- log2(mat_B + 1)
          logFC_text <- " (log2 +1)"
        }
        
        mat_A <- as.data.frame(mat_A)
        mat_B <- as.data.frame(mat_B)
        
        # Get avg expression of each gene
        mat_A$AvgExpression_GroupA <- rowMeans(mat_A)
        mat_B$AvgExpression_GroupB <- rowMeans(mat_B)
        
        # Add Gene Names as column
        mat_A$GeneSymbol <- rownames(mat_A)
        mat_B$GeneSymbol <- rownames(mat_B)
        
        # Merge average columns
        AvgExpr_Table <- merge(mat_A[,c("AvgExpression_GroupA","GeneSymbol")],
                               mat_B[,c("AvgExpression_GroupB","GeneSymbol")],
                               by = "GeneSymbol",
                               all = T)
        
        colnames(AvgExpr_Table)[2] <- paste("Average Expression ", A_choice, sep = "")
        colnames(AvgExpr_Table)[3] <- paste("Average Expression ", B_choice, sep = "")
        AvgExpr_Table
        
      })
      
      output$AvggeneScatterTable <- renderDataTable({
        
        AvgExpr_Table <- AvgGeneScatterTable_react()
        
        #table output
        DT::datatable(AvgExpr_Table,
                      options = list(keys = TRUE,
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = c("10", "25", "50", "100")),
                      rownames = F)
        
      })
      
      #render ssGSEA table
      output$ssGSEAtable <- DT::renderDataTable({
        req(ssGSEAfunc())
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        
        if (length(input$GeneSetTable_rows_selected) > 0){
          ssgsea <- ssGSEAfunc()
          ssgsea2 <- as.data.frame(t(ssgsea))
          samporder <- meta[,1]
          ssgsea3 <- as.data.frame(ssgsea2[samporder,])
          colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
          rownames(ssgsea3) <- samporder
          ssgsea4 <- ssgsea3 %>%
            mutate(type = case_when(
              rownames(ssgsea3) == meta[,1] ~ meta[,metacol],
            ))
          ssgsea4 <- ssgsea4 %>%
            relocate(type)
          ssgsea4$type <- as.factor(ssgsea4$type)
          colnames(ssgsea4)[which(colnames(ssgsea4) == "type")] <- metacol
          #table output
          DT::datatable(ssgsea4,
                        options = list(keys = TRUE,
                                       searchHighlight = TRUE,
                                       pageLength = 10,
                                       lengthMenu = c("10", "25", "50", "100")))
        }
        #if (ncol(meta) > 2) {
        #  metacol <- input$GSEAmetaCol
        #}
        #else if (ncol(meta) == 2) {
        #  metacol <- colnames(meta)[2]
        #}
        #if (input$tables == 3) {
        #  if (length(input$tab2table_rows_selected) > 0){
        #    ssgsea <- ssGSEAfunc()
        #    ssgsea2 <- as.data.frame(t(ssgsea))
        #    samporder <- meta[,1]
        #    ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        #    colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        #    rownames(ssgsea3) <- samporder
        #    ssgsea4 <- ssgsea3 %>%
        #      mutate(type = case_when(
        #        rownames(ssgsea3) == meta[,1] ~ meta[,metacol],
        #      ))
        #    ssgsea4 <- ssgsea4 %>%
        #      relocate(type)
        #    ssgsea4$type <- as.factor(ssgsea4$type)
        #    colnames(ssgsea4)[which(colnames(ssgsea4) == "type")] <- metacol
        #    #table output
        #    DT::datatable(ssgsea4,
        #                  options = list(keys = TRUE,
        #                                 searchHighlight = TRUE,
        #                                 pageLength = 10,
        #                                 lengthMenu = c("10", "25", "50", "100")))
        #  }
        #}
        #else if (input$tables == 1) {
        #  if (length(input$msigdbTable_rows_selected) > 0){
        #    ssgsea <- ssGSEAfunc()
        #    ssgsea2 <- as.data.frame(t(ssgsea))
        #    samporder <- meta[,1]
        #    ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        #    colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        #    rownames(ssgsea3) <- samporder
        #    ssgsea4 <- ssgsea3 %>%
        #      mutate(type = case_when(
        #        rownames(ssgsea3) == meta[,1] ~ meta[,metacol],
        #      ))
        #    ssgsea4 <- ssgsea4 %>%
        #      relocate(type)
        #    ssgsea4$type <- as.factor(ssgsea4$type)
        #    colnames(ssgsea4)[which(colnames(ssgsea4) == "type")] <- metacol
        #    #table output
        #    DT::datatable(ssgsea4,
        #                  options = list(keys = TRUE,
        #                                 searchHighlight = TRUE,
        #                                 pageLength = 10,
        #                                 lengthMenu = c("10", "25", "50", "100")))
        #  }
        #}
        #else if (input$tables == 5) {
        #  if (length(input$GStable.u_rows_selected) > 0){
        #    ssgsea <- ssGSEAfunc()
        #    ssgsea2 <- as.data.frame(t(ssgsea))
        #    samporder <- meta[,1]
        #    ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        #    colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        #    rownames(ssgsea3) <- samporder
        #    ssgsea4 <- ssgsea3 %>%
        #      mutate(type = case_when(
        #        rownames(ssgsea3) == meta[,1] ~ meta[,metacol],
        #      ))
        #    ssgsea4 <- ssgsea4 %>%
        #      relocate(type)
        #    ssgsea4$type <- as.factor(ssgsea4$type)
        #    colnames(ssgsea4)[which(colnames(ssgsea4) == "type")] <- metacol
        #    #table output
        #    DT::datatable(ssgsea4,
        #                  options = list(keys = TRUE,
        #                                 searchHighlight = TRUE,
        #                                 pageLength = 10,
        #                                 lengthMenu = c("10", "25", "50", "100")))
        #  }
        #}
      })
      
      
      ####----Plots----####
      
      enrichplot0_react <- reactive({
        
        req(datasetInput())
        gmt <- gmt_react()
        res <- datasetInput()
        gseaplot2(res,
                  geneSetID = 1:length(res@result[["ID"]]),
                  title = paste(unique(gmt[,1]),collapse = "\n"),
                  pvalue_table = F)
        
      })
      
      #render GSEA plot
      output$enrichplot0 <- renderPlot({
        p <- enrichplot0_react()
        p
      })
      
      GeneSetName_React <- reactive({
        
        names(gs_react())
        #if (input$tables == 1) {
        #  if (length(input$msigdbTable_rows_selected) > 0) {
        #    GS <- as.character(msigdb.gsea2()[input$msigdbTable_rows_selected,3])
        #    GS
        #  }
        #} else if (input$tables == 3) {
        #  if (length(input$tab2table_rows_selected) > 0) {
        #    GS <- as.character(GeneSet2()[input$tab2table_rows_selected,1])
        #    GS
        #  }
        #} else if (input$tables == 5) {
        #  if (length(input$GStable.u_rows_selected) > 0) {
        #    GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        #    GS
        #  }
        #}
        
      })
      
      gsea_heat_anno <- reactive({
        
        if (isTruthy(input$GSEAGroupColMulti)) {
          meta <- meta_react()
          dataset <- gsea_heat_data()
          meta <- meta[which(meta[,1] %in% colnames(dataset)),]
          meta <- meta[match(colnames(dataset),meta[,1]),]
          rownames(meta) <- meta[,1]
          metacol <- metacol_reactVal()
          if (length(colnames(meta)) == 2) {
            metacol <- colnames(meta)[2]
          }
          
          if (input$volcanoCompChoice3 == "Two groups") {
            groupA <- input$comparisonA
            groupB <- input$comparisonB
            
            
            anno_cols <- input$GSEAGroupColMulti
            meta_sub <- meta[,anno_cols, drop = F]
            meta_sub[,anno_cols] <- as.data.frame(lapply(meta_sub[,anno_cols, drop = F], function(x) {
              if (!is.numeric(x) | all(suppressWarnings(as.numeric(x))%%1==0)) {
                return(factor(x))
              } else {
                return(x)
              }
            }))
            #meta_sub[,input$input$GSEAGroupColMulti] <- as.data.frame(lapply(meta_sub[,input$input$GSEAGroupColMulti, drop = F], factor))
            colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                        name = anno_cols,
                                                        which = 'col'
            )
            
            #meta_sub <- meta[,input$GSEAGroupColMulti, drop = F]
            #meta_sub[,input$GSEAGroupColMulti] <- as.data.frame(lapply(meta_sub[,input$GSEAGroupColMulti, drop = F], factor))
            #colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
            #                                            name = input$GSEAGroupColMulti,
            #                                            which = 'col'
            #)
          } else if (input$volcanoCompChoice3 == "One group") {
            groupA <- input$comparisonA_one
            groupB <- paste0("Not_",groupA)
            meta[,paste0(metacol,"_Dichot")] <- NA
            meta[which(meta[,metacol] == groupA),paste0(metacol,"_Dichot")] <- groupA
            meta[which(meta[,metacol] != groupA),paste0(metacol,"_Dichot")] <- paste0("Not_",groupA)
            metacol <- paste0(metacol,"_Dichot")
            
            anno_cols <- c(metacol,input$GSEAGroupColMulti)
            meta_sub <- meta[,anno_cols, drop = F]
            meta_sub[,anno_cols] <- as.data.frame(lapply(meta_sub[,anno_cols, drop = F], function(x) {
              if (!is.numeric(x) | all(suppressWarnings(as.numeric(x))%%1==0)) {
                return(factor(x))
              } else {
                return(x)
              }
            }))
            #meta_sub[,input$input$GSEAGroupColMulti] <- as.data.frame(lapply(meta_sub[,input$input$GSEAGroupColMulti, drop = F], factor))
            colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                        name = anno_cols,
                                                        which = 'col'
            )
            
            
            #meta_sub <- meta[,c(metacol,input$GSEAGroupColMulti), drop = F]
            #meta_sub[,c(metacol,input$GSEAGroupColMulti)] <- as.data.frame(lapply(meta_sub[,c(metacol,input$GSEAGroupColMulti), drop = F], factor))
            #colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
            #                                            name = c(metacol,input$GSEAGroupColMulti),
            #                                            which = 'col'
            #)
          }
          
          colAnn
        } else {
          meta <- meta_react()
          if (ncol(meta) == 2) {
            rownames(meta) <- meta[,1]
            dataset <- gsea_heat_data()
            meta <- meta[which(meta[,1] %in% colnames(dataset)),]
            meta <- meta[match(colnames(dataset),meta[,1]),]
            rownames(meta) <- meta[,1]
            meta_sub <- meta[,2, drop = F]
            meta_sub[,1] <- as.factor(meta_sub[,1])
            colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                        name = colnames(meta_sub)[2],
                                                        which = 'col'
            )
            colAnn
          } else {
            colAnn <- NULL
            colAnn
          }
        }
        
      })
      
      gsea_heat_data <- reactive({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        A <- A_raw()
        GS <- GeneSetName_React()
        #groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
        #groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
        
        if (input$volcanoCompChoice3 == "Two groups") {
          groupA <- input$comparisonA
          groupB <- input$comparisonB
        } else if (input$volcanoCompChoice3 == "One group") {
          groupA <- input$comparisonA_one
          groupB <- paste0("Not_",groupA)
          meta[,paste0(metacol,"_Dichot")] <- NA
          meta[which(meta[,metacol] == groupA),paste0(metacol,"_Dichot")] <- groupA
          meta[which(meta[,metacol] != groupA),paste0(metacol,"_Dichot")] <- paste0("Not_",groupA)
          metacol <- paste0(metacol,"_Dichot")
        }
        
        groupA <- meta[which(meta[,metacol] == groupA),1]
        groupB <- meta[which(meta[,metacol] == groupB),1]
        
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
        genes2 <- strsplit(genes1,"/")
        genes3 <- as.data.frame(genes2, col.names = "genes")
        gene_symbol <- genes3$genes
        ## convert expression matrix to numeric
        class(A) <- "numeric"
        ## Transforming data
        A <- A[,c(groupA,groupB)]
        exp.mat1 = log2(A + 1) # log
        exp.mat2 = apply(exp.mat1, 1, scale); # z score
        exp.mat3 = apply(exp.mat2, 1, rev); # transpose
        colnames(exp.mat3) = colnames(A) # set the column name
        exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
        exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
        # reassign data
        dataset <- exp.mat5
        dataset
        
      })
      
      heatmap0_react <- reactive({
        
        #meta <- meta_react()
        #metacol <- metacol_reactVal()
        #A <- expr_mat_react()
        colAnno <- gsea_heat_anno()
        dataset <- gsea_heat_data()
        
        p <- suppressMessages(ComplexHeatmap::Heatmap(dataset,
                                                      top_annotation = colAnno,
                                                      #clustering_method_rows = clust_method,
                                                      #show_row_names = row_names_choice, show_column_names = col_names_choice,
                                                      cluster_rows = F, cluster_columns = F,
                                                      #row_names_gp = gpar(fontsize = row_font), column_names_gp = gpar(fontsize = col_font),
                                                      heatmap_legend_param = list(title = "Expression"),
                                                      border = F))
        draw(p, padding = unit(c(100, 50, 2, 2), "mm")) # unit(c(bottom,left,right,top))
        
      })
      
      #render heatmap
      output$heatmap0 <- renderPlot({
        
        p <- heatmap0_react()
        p
        
      })
      
      #render MVG heatmap - per sample - 1
      output$heatmap1 <- renderPlot({
        
        heat <- MVGheatmap_react()
        heat
        
      })
      
      #render custom heatmap - per sample - 2
      output$heatmap2 <- renderPlot({
        
        heat <- CustomHeatmap_react()
        heat
        
      })
      
      # Render DEG Expression heatmap - #3
      output$heatmap3 <- renderPlot({
        
        heat <- DEGHeatmap_react()
        heat
        
        
      })
      
      # Render Average Expression heatmap - custom
      output$avgheatmap1Cust <- renderPlot({
        
        heat <- AvgExprHeatmap_react()
        heat
        
        
      })
      
      output$ssgseaheatmap <- renderPlot({
        
        heat <- ssgseaheatmap_react()
        heat
        
      })
      
      output$ssgseaheatmap2 <- renderPlot({
        
        heat <- ssgseaheatmap2_react()
        heat
        
      })
      
      output$ssgseaheatmap3 <- renderPlot({
        
        heat <- ssgseaheatmap3_react()
        heat
        
      })
      
      #render MA plot
      output$MAPlot1 <- renderPlot({
        #title_font <- input$VolMATitleSize
        axis_font <- input$VolMAAxisSize
        top2 <- topgenereact()
        if (input$volcanoCompChoice == "DESeq2") {
          top2$AveExpr <- log2(top2$AveExpr+1)
        }
        #add color categories based on FC and pval
        top2['threshold'] <- "none"
        top2[which(top2$logFC > abs(input$fc_cutoff)), "threshold"] <- "up"
        top2[which(top2$logFC < -abs(input$fc_cutoff)), "threshold"] <- "down"
        upRed <- "lightcoral"
        dnBlue <- "cadetblue3"
        mdGray <- "gray70"
        top2u <- top2[order(top2[,1], decreasing = TRUE),]
        top2d <- top2[order(top2[,1]),]
        top_hits_up <- top2u[head(which(top2u$logFC > abs(input$fc_cutoff)), n = input$top_x),]
        top_hits_dn <- top2d[head(which(top2d$logFC < -abs(input$fc_cutoff)), n = input$top_x),]
        #select genes to label based on user selection
        genesel.s <- NULL
        genesel.t <- NULL
        genesel.u <- NULL
        genesel.s <- unlist(strsplit(input$userGeneSelec2, " "))
        genesel.t <- unlist(strsplit(input$userGeneSelec2, "\t"))
        genesel.u <- input$userGeneSelec
        genesel.text <- c(genesel.s,genesel.t,genesel.u)
        top2_selec <- top2 %>%
          filter(GeneName %in% genesel.text)
        
        
        
        x <- ggplot(data = top2, aes(x = AveExpr, y = logFC)) +
          geom_point(size = 2, shape = 16) +
          theme_light(base_size = 16)
        #colors
        x <- x + aes(color = threshold) +
          scale_color_manual(values = c("up" = upRed,"down" = dnBlue, "none" = mdGray))
        x <- x + geom_hline(yintercept = c(-abs(input$fc_cutoff),abs(input$fc_cutoff)), linetype="dashed", color="gray20")
        x <- x + geom_text_repel(
          data =  top_hits_up,
          aes(label = rownames(top_hits_up)),
          color="gray20",
          size = 6,
          #nudge_x = 0.2,
          #nudge_y=0.2,
          box.padding = unit(0.9, "lines"),
          point.padding = unit(.3+4*0.1, "lines"),
          max.overlaps = 50)
        x <- x + geom_text_repel(
          data =  top_hits_dn,
          aes(label = rownames(top_hits_dn)),
          color="gray20",
          size = 6,
          #nudge_x = 0.2,
          #nudge_y=0.2,
          box.padding = unit(0.9, "lines"),
          point.padding = unit(.3+4*0.1, "lines"),
          max.overlaps = 50)
        x <- x + geom_text_repel(
          data =  top2_selec,
          aes(label = rownames(top2_selec)),
          color="gray20",
          size = 6,
          #nudge_x = 0.2,
          #nudge_y=0.2,
          box.padding = unit(0.9, "lines"),
          point.padding = unit(.3+4*0.1, "lines"),
          max.overlaps = 50)
        #coloring selected points
        x <- x + geom_point(data = top2_selec,
                            aes(x = AveExpr, y = logFC),
                            pch = 21,
                            color = "black",
                            size = 2)
        x <- x + geom_point(data = top_hits_dn,
                            aes(x = AveExpr, y = logFC),
                            pch = 21,
                            color = "black",
                            size = 2)
        x <- x + geom_point(data = top_hits_up,
                            aes(x = AveExpr, y = logFC),
                            pch = 21,
                            color = "black",
                            size = 2)
        x <- x + theme(legend.position="none")
        x <- x + labs(x = "Average Expression", y = "log2FC")
        x <- x + theme(axis.text = element_text(size=18))
        x <- x + theme(axis.title = element_text(size = axis_font))
        x
      })
      
      # DGE --------------------------------------------------------------------
      
      ## Boxplot ---------------------------------------------------------------
      
      output$rendBPgroupCriteria2 <- renderUI({
        
        meta <- meta_react()
        GroupChoices <- colnames(meta_react())[-1]
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        if (isTruthy(input$SubsetCol)) {
          if (input$SubsetCol != "Select All Samples") {
            GroupChoices <- GroupChoices[which(GroupChoices!=input$SubsetCol)]
          }
        }
        selectInput("BPgroupCriteria2","Grouping Criteria:",choices = GroupChoices, selected = metacol)
        
      })
      
      output$rendBPgroupSelection2 <- renderUI({
        
        req(input$BPgroupCriteria2)
        req(meta_react())
        meta <- meta_react()
        groupCrit <- input$BPgroupCriteria2
        if (input$BPremoveSingles2 == T) {
          tab <- table(meta[,groupCrit])
          meta <- meta[meta[,groupCrit] %in% names(tab[tab >= 2]), ]
        }
        GroupSelec <- unique(meta[,groupCrit])
        selectInput("BPgroupSelection2","Select Groups:",choices = GroupSelec, selected = GroupSelec, multiple = T)
        
      })
      
      BP_Feature_Choices2 <- reactive({
        
        req(expr_react())
        req(meta_react())
        Features <- NULL
        FeatCat <- input$BPFeatureCategory2
        if (FeatCat == "Matrix Features") {
          Features <- rownames(expr_react())
          Features
        } else if (FeatCat == "Meta Features") {
          #Features <- colnames(aligned_meta_file())[-1]
          Features <- colnames(meta_react())[-1]
          Features
        } else if (FeatCat == "Immune Deconvolution Features") {
          Features <- ImmDeconv_uncorr_react()[,1]
          Features
        } else {
          Features <- NULL
          Features
        }
        
      })
      
      shiny::observe({
        req(BP_Feature_Choices2())
        updateSelectizeInput(session = session, inputId = "BPFeatSelection2",
                             choices = BP_Feature_Choices2(),
                             server = T)
        
      })
      
      CohortBPPlot_df_react2 <- reactive({
        
        meta <- meta_react()
        mat <- expr_react()
        req(input$BPgroupCriteria2)
        groupCrit <- input$BPgroupCriteria2
        FeatCat <- input$BPFeatureCategory2
        Feature <- input$BPFeatSelection2
        NameCol <- colnames(meta)[1]
        removeSingles <- input$BPremoveSingles2
        BPlog <- input$BPlogOpt2
        if (FeatCat == "Matrix Features") {
          feat <- unlist(mat[Feature,])
          featdf <- data.frame(
            NameCol = names(feat),
            Feature = unname(feat))
          colnames(featdf) <- c(NameCol,Feature)
          meta <- merge(meta,featdf)
        } #else if (FeatCat == "Immune Deconvolution Features") {
        #mat <- ImmDeconv_uncorr_react()
        #featdf <- mat[which(mat[,1] == Feature),]
        #rownames(featdf) <- featdf[,1]
        #featdf <- featdf[,-1]
        #featdf <- as.data.frame(t(featdf))
        #featdf[,NameCol] <- rownames(featdf)
        #meta <- merge(meta,featdf)
        #} 
        
        if (isTruthy(Feature)) {
          if (Feature %in% colnames(meta)) {
            meta <- meta %>% select(any_of(c(NameCol,groupCrit,Feature)))
            if (removeSingles == T) {
              tab <- table(meta[,groupCrit])
              meta <- meta[meta[,groupCrit] %in% names(tab[tab >= 2]), ]
            }
            meta[,Feature] <- as.numeric(meta[,Feature])
            meta <- meta[which(!is.na(meta[,Feature])),]
            meta <- meta[which(meta[,Feature]!=Inf & meta[,Feature]!=-Inf),]
            meta[,groupCrit] <- as.factor(meta[,groupCrit])
            if (BPlog == "Log2") {
              meta[,Feature] <- log2(meta[,Feature])
            } else if (BPlog == "Log2+1") {
              meta[,Feature] <- log2(meta[,Feature]+1)
            } else if (BPlog == "Log10") {
              meta[,Feature] <- log10(meta[,Feature])
            } else if (BPlog == "Log10+1") {
              meta[,Feature] <- log10(meta[,Feature]+1)
            }
            #meta
          }
          meta <- meta %>%
            group_by(!!sym(groupCrit)) %>%
            mutate(Outlier = ifelse(rstatix::is_outlier(!!sym(Feature)),TRUE,FALSE)) %>%
            as.data.frame()
          meta
        }
      })
      
      output$DataExplor_Box_plot_df2 <- DT::renderDataTable({
        
        df <- CohortBPPlot_df_react2()
        DT::datatable(df,
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 10,
                                     scrollX = T),
                      rownames = F)
        
      })
      
      CohortBPPlot_react2 <- reactive({
        
        plotdf_full <- CohortBPPlot_df_react2()
        groupCrit <- input$BPgroupCriteria2
        FeatCat <- input$BPFeatureCategory2
        Feature <- input$BPFeatSelection2
        NameCol <- colnames(plotdf_full)[1]
        BPplottheme <- input$BPTheme2
        StatMethod <- input$BPplotstatComp2
        dotChoice <- input$BPplotsampledots2
        dotSize <- input$BPplotDotSize2
        bpFlip <- input$BPflipBP2
        BPorVI <- input$BPorViolin2
        Xaxis_font <- input$BPplot1XAxisSize2              # Axis font size
        Yaxis_font <- input$BPplot1YAxisSize2              # Axis font size
        Yaxis_lim <- input$BPplot1YAxisLim2
        hjust_orient <- 1                                # Initial hjust
        axis_orient <- as.numeric(input$BPxAxisOrient2)  # X-axis label orientation
        if (axis_orient == 0) {                          # Adjust hjust if orientation is 0
          hjust_orient <- 0.5
        }
        BPorder <- input$BPplotXaxOrder2
        BPGroupSelect <- input$BPgroupSelection2
        
        if (isTruthy(Feature)) {
          plotdf <- plotdf_full[,c(NameCol,groupCrit,Feature)]
          plotdf <- plotdf[which(plotdf[,groupCrit] %in% BPGroupSelect),]
          colnames(plotdf) <- c("SampleName","Group","Feature")
          
          if (BPorder == "Descending"){
            barp <- ggplot(data = plotdf, aes(x=reorder(Group,-Feature, FUN = median),y=Feature, fill=Group))
            plotdf_dots <- plotdf
            plotdf_dots$Group <- reorder(plotdf_dots$Group,-plotdf_dots$Feature, FUN = median)
            plotdf_dots$xj <- jitter(as.numeric(factor(plotdf_dots$Group)))
          }
          if (BPorder == "Ascending"){
            barp <- ggplot(data = plotdf, aes(x=reorder(Group,Feature, FUN = median),y=Feature, fill=Group))
            plotdf_dots <- plotdf
            plotdf_dots$Group <- reorder(plotdf_dots$Group,plotdf_dots$Feature, FUN = median)
            plotdf_dots$xj <- jitter(as.numeric(factor(plotdf_dots$Group)))
          }
          if (BPorder == "Not Specificed"){
            barp <- ggplot(data = plotdf, aes(x=Group,y=Feature, fill=Group))
            plotdf_dots <- plotdf
            plotdf_dots$xj <- jitter(as.numeric(factor(plotdf_dots$Group)))
          }
          if (BPorVI == "Box Plot") {
            barp <- barp + geom_boxplot(width = 0.5, lwd = 1)
          }
          if (BPorVI == "Violin Plot") {
            barp <- barp + geom_violin() +
              stat_summary(fun=median, geom="crossbar", width=0.5, color="black")
          }
          if (isTruthy(Yaxis_lim)) {
            barp <- barp +
              ylim(paste0(as.numeric(strsplit(Yaxis_lim,",")[[1]][1]),as.numeric(strsplit(Yaxis_lim,",")[[1]][2])))
          }
          
          barp <- barp +
            get(BPplottheme)() +
            labs(#title = BPTitle_in,
              x = groupCrit, y = Feature,
              fill = groupCrit)
          if (StatMethod != "none") {
            barp <- barp + ggpubr::stat_compare_means(method = StatMethod)
          }
          if (dotChoice) {
            barp <- barp + geom_point(data = plotdf_dots, aes(x=xj), col="grey14", size=dotSize)
          }
          barp <- barp + theme(axis.text.x = element_text(size = Xaxis_font,angle = axis_orient, hjust = hjust_orient),
                               axis.title.x = element_text(size = Xaxis_font),
                               axis.text.y = element_text(size = Yaxis_font),
                               axis.title.y = element_text(size = Yaxis_font),
                               legend.position = "none")
          if (bpFlip) {
            barp <- barp + coord_flip()
          }
          barp
        }
        
      })
      
      output$uncorrected_Box_plot2 <- renderPlot({
        
        barp <- CohortBPPlot_react2()
        barp
        
      })
      
      
      #render boxplot
      output$boxplot1 <- renderPlot({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        expr <- expr_react()
        geneList <- geneList_raw()
        title_font <- input$boxplot2TitleSize
        axis_font <- input$boxplot2AxisSize
        #if (ncol(meta) > 2) {
        #  metacol <- input$BoxPlotMetaColSelec
        #}
        #else if (ncol(meta) == 2) {
        #  metacol <- colnames(meta)[2]
        #}
        
        if (length(input$GeneListTable_rows_selected) > 0){
          gene <- geneList[input$GeneListTable_rows_selected, 1]
          min <- min(log2(expr[gene,] + 1.0))
          max <- max(log2(expr[gene,] + 1.0))
          meta_temp <- meta
          rownames(meta_temp) <- meta[,1]
          meta_temp <- meta_temp[,metacol,drop = F]
          meta_temp[,metacol] <- as.factor(meta_temp[,metacol])
          #meta_temp <- meta_temp %>%
          #  select(Group)
          data = merge(t(expr[gene,]), meta_temp, by=0)
          colnames(data) = c("SampleName", "GeneExpr", "Cluster")
          ggplot(data, aes(x=Cluster, y=log2(GeneExpr + 1.0))) +
            geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=1) +
            geom_dotplot(binaxis='y', stackdir='center', dotsize=input$boxplotDot) +
            stat_compare_means(label = "p.signif") + ylim(min * 0.9, max * 1.3) +
            theme_bw() +
            labs(title= paste(gene, "Expression (log2)")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  axis.title = element_text(size = axis_font),
                  plot.title = element_text(size = title_font))
        }
      })
      
      
      # Data Exploration -------------------------------------------------------
      
      ## Boxplot ---------------------------------------------------------------
      
      output$rendBPgroupCriteria <- renderUI({
        
        meta <- meta_react()
        GroupChoices <- colnames(meta_react())[-1]
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        if (isTruthy(input$SubsetCol)) {
          if (input$SubsetCol != "Select All Samples") {
            GroupChoices <- GroupChoices[which(GroupChoices!=input$SubsetCol)]
          }
        }
        selectInput("BPgroupCriteria","Grouping Criteria:",choices = GroupChoices, selected = metacol)
        
      })
      
      output$rendBPgroupSelection <- renderUI({
        
        req(input$BPgroupCriteria)
        req(meta_react())
        meta <- meta_react()
        groupCrit <- input$BPgroupCriteria
        if (input$BPremoveSingles == T) {
          tab <- table(meta[,groupCrit])
          meta <- meta[meta[,groupCrit] %in% names(tab[tab >= 2]), ]
        }
        GroupSelec <- unique(meta[,groupCrit])
        selectInput("BPgroupSelection","Select Groups:",choices = GroupSelec, selected = GroupSelec, multiple = T)
        
      })
      
      BP_Feature_Choices <- reactive({
        
        req(expr_react())
        req(meta_react())
        Features <- NULL
        FeatCat <- input$BPFeatureCategory
        if (FeatCat == "Matrix Features") {
          Features <- rownames(expr_react())
          Features
        } else if (FeatCat == "Meta Features") {
          #Features <- colnames(aligned_meta_file())[-1]
          Features <- colnames(meta_react())[-1]
          Features
        } else if (FeatCat == "Immune Deconvolution Features") {
          Features <- ImmDeconv_uncorr_react()[,1]
          Features
        } else {
          Features <- NULL
          Features
        }
        
      })
      
      shiny::observe({
        
        updateSelectizeInput(session = session, inputId = "BPFeatSelection",
                             choices = BP_Feature_Choices(),
                             server = T)
        
      })
      
      CohortBPPlot_df_react <- reactive({
        
        meta <- meta_react()
        mat <- expr_react()
        req(input$BPgroupCriteria)
        groupCrit <- input$BPgroupCriteria
        FeatCat <- input$BPFeatureCategory
        Feature <- input$BPFeatSelection
        NameCol <- colnames(meta)[1]
        removeSingles <- input$BPremoveSingles
        BPlog <- input$BPlogOpt
        if (FeatCat == "Matrix Features") {
          feat <- unlist(mat[Feature,])
          featdf <- data.frame(
            NameCol = names(feat),
            Feature = unname(feat))
          colnames(featdf) <- c(NameCol,Feature)
          meta <- merge(meta,featdf)
        } #else if (FeatCat == "Immune Deconvolution Features") {
        #mat <- ImmDeconv_uncorr_react()
        #featdf <- mat[which(mat[,1] == Feature),]
        #rownames(featdf) <- featdf[,1]
        #featdf <- featdf[,-1]
        #featdf <- as.data.frame(t(featdf))
        #featdf[,NameCol] <- rownames(featdf)
        #meta <- merge(meta,featdf)
        #} 
        if (isTruthy(Feature)) {
          if (Feature %in% colnames(meta)) {
            meta <- meta %>% select(any_of(c(NameCol,groupCrit,Feature)))
            if (removeSingles == T) {
              tab <- table(meta[,groupCrit])
              meta <- meta[meta[,groupCrit] %in% names(tab[tab >= 2]), ]
            }
            meta[,Feature] <- as.numeric(meta[,Feature])
            meta <- meta[which(!is.na(meta[,Feature])),]
            meta <- meta[which(meta[,Feature]!=Inf & meta[,Feature]!=-Inf),]
            meta[,groupCrit] <- as.factor(meta[,groupCrit])
            if (BPlog == "Log2") {
              meta[,Feature] <- log2(meta[,Feature])
            } else if (BPlog == "Log2+1") {
              meta[,Feature] <- log2(meta[,Feature]+1)
            } else if (BPlog == "Log10") {
              meta[,Feature] <- log10(meta[,Feature])
            } else if (BPlog == "Log10+1") {
              meta[,Feature] <- log10(meta[,Feature]+1)
            }
            #meta
          }
          meta <- meta %>%
            group_by(!!sym(groupCrit)) %>%
            mutate(Outlier = ifelse(rstatix::is_outlier(!!sym(Feature)),TRUE,FALSE)) %>%
            as.data.frame()
          meta
        }
      })
      
      output$DataExplor_Box_plot_df <- DT::renderDataTable({
        
        df <- CohortBPPlot_df_react()
        DT::datatable(df,
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 10,
                                     scrollX = T),
                      rownames = F)
        
      })
      
      CohortBPPlot_react <- reactive({
        
        plotdf_full <- CohortBPPlot_df_react()
        groupCrit <- input$BPgroupCriteria
        FeatCat <- input$BPFeatureCategory
        Feature <- input$BPFeatSelection
        NameCol <- colnames(plotdf_full)[1]
        BPplottheme <- input$BPTheme
        StatMethod <- input$BPplotstatComp
        dotChoice <- input$BPplotsampledots
        dotSize <- input$BPplotDotSize
        bpFlip <- input$BPflipBP
        BPorVI <- input$BPorViolin
        Xaxis_font <- input$BPplot1XAxisSize              # Axis font size
        Yaxis_font <- input$BPplot1YAxisSize              # Axis font size
        Yaxis_lim <- input$BPplot1YAxisLim
        hjust_orient <- 1                                # Initial hjust
        axis_orient <- as.numeric(input$BPxAxisOrient)  # X-axis label orientation
        if (axis_orient == 0) {                          # Adjust hjust if orientation is 0
          hjust_orient <- 0.5
        }
        BPorder <- input$BPplotXaxOrder
        BPGroupSelect <- input$BPgroupSelection
        
        if (isTruthy(Feature)) {
          plotdf <- plotdf_full[,c(NameCol,groupCrit,Feature)]
          plotdf <- plotdf[which(plotdf[,groupCrit] %in% BPGroupSelect),]
          colnames(plotdf) <- c("SampleName","Group","Feature")
          
          if (BPorder == "Descending"){
            barp <- ggplot(data = plotdf, aes(x=reorder(Group,-Feature, FUN = median),y=Feature, fill=Group))
            plotdf_dots <- plotdf
            plotdf_dots$Group <- reorder(plotdf_dots$Group,-plotdf_dots$Feature, FUN = median)
            plotdf_dots$xj <- jitter(as.numeric(factor(plotdf_dots$Group)))
          }
          if (BPorder == "Ascending"){
            barp <- ggplot(data = plotdf, aes(x=reorder(Group,Feature, FUN = median),y=Feature, fill=Group))
            plotdf_dots <- plotdf
            plotdf_dots$Group <- reorder(plotdf_dots$Group,plotdf_dots$Feature, FUN = median)
            plotdf_dots$xj <- jitter(as.numeric(factor(plotdf_dots$Group)))
          }
          if (BPorder == "Not Specificed"){
            barp <- ggplot(data = plotdf, aes(x=Group,y=Feature, fill=Group))
            plotdf_dots <- plotdf
            plotdf_dots$xj <- jitter(as.numeric(factor(plotdf_dots$Group)))
          }
          if (BPorVI == "Box Plot") {
            barp <- barp + geom_boxplot(width = 0.5, lwd = 1)
          }
          if (BPorVI == "Violin Plot") {
            barp <- barp + geom_violin() +
              stat_summary(fun=median, geom="crossbar", width=0.5, color="black")
          }
          if (isTruthy(Yaxis_lim)) {
            barp <- barp +
              ylim(paste0(as.numeric(strsplit(Yaxis_lim,",")[[1]][1]),as.numeric(strsplit(Yaxis_lim,",")[[1]][2])))
          }
          
          barp <- barp +
            get(BPplottheme)() +
            labs(#title = BPTitle_in,
              x = groupCrit, y = Feature,
              fill = groupCrit)
          if (StatMethod != "none") {
            barp <- barp + ggpubr::stat_compare_means(method = StatMethod)
          }
          if (dotChoice) {
            barp <- barp + geom_point(data = plotdf_dots, aes(x=xj), col="grey14", size=dotSize)
          }
          barp <- barp + theme(axis.text.x = element_text(size = Xaxis_font,angle = axis_orient, hjust = hjust_orient),
                               axis.title.x = element_text(size = Xaxis_font),
                               axis.text.y = element_text(size = Yaxis_font),
                               axis.title.y = element_text(size = Yaxis_font),
                               legend.position = "none")
          if (bpFlip) {
            barp <- barp + coord_flip()
          }
          barp
        }
        
      })
      
      output$uncorrected_Box_plot <- renderPlot({
        
        barp <- CohortBPPlot_react()
        barp
        
      })
      
      
      #render boxplot
      #output$boxplot3 <- renderPlot({
      #  
      #  meta <- meta_react()
      #  metacol <- metacol_reactVal()
      #  if (length(colnames(meta)) == 2) {
      #    metacol <- colnames(meta)[2]
      #  }
      #  expr <- expr_react()
      #  geneList <- geneList_raw()
      #  title_font <- input$boxplot1TitleSize
      #  axis_font <- input$boxplot1AxisSize
      #  #if (ncol(meta) > 2) {
      #  #  metacol <- input$BoxPlot1MetaCol
      #  #}
      #  #else if (ncol(meta) == 2) {
      #  #  metacol <- colnames(meta)[2]
      #  #}
      #  
      #  if (length(input$GeneListTable2_rows_selected) > 0){
      #    gene <- geneList[input$GeneListTable2_rows_selected, 1]
      #    min <- min(log2(expr[gene,] + 1.0))
      #    max <- max(log2(expr[gene,] + 1.0))
      #    meta_temp <- meta
      #    rownames(meta_temp) <- meta[,1]
      #    meta_temp <- meta_temp[,metacol,drop = F]
      #    meta_temp[,metacol] <- as.factor(meta_temp[,metacol])
      #    #meta_temp <- meta_temp %>%
      #    #  select(Group)
      #    data = merge(t(expr[gene,]), meta_temp, by=0)
      #    colnames(data) = c("SampleName", "GeneExpr", "Cluster")
      #    ggplot(data, aes(x=Cluster, y=log2(GeneExpr + 1.0))) +
      #      geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=1) +
      #      geom_dotplot(binaxis='y', stackdir='center', dotsize=input$boxplotDot2) +
      #      stat_compare_means(label = "p.signif") + ylim(min * 0.9, max * 1.3) +
      #      theme_bw() +
      #      labs(title= paste(gene, "Expression (log2)")) +
      #      theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
      #            axis.title = element_text(size = axis_font),
      #            plot.title = element_text(size = title_font))
      #  }
      #})
      
      ####----Bar Plot----####
      barplot_react <- reactive({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        expr <- expr_react()
        geneList <- geneList_raw()
        title_font <- input$barplot1TitleSize            # Title font size
        axis_font <- input$barplot1AxisSize              # Axis font size
        hjust_orient <- 1                                # Initial hjust
        axis_orient <- as.numeric(input$barxAxisOrient)  # X-axis label orientation
        if (axis_orient == 0) {                          # Adjust hjust if orientation is 0
          hjust_orient <- 0.5
        }
        bar_type <- input$errorbarplot                   # Error bar type
        logchoice <- input$log2barplot                   # Log expression data option
        colorin <- input$barplotColoCodes
        colorout <- input$barplotColoCodesOut
        bpylim <- input$barPlotYlim                      # Y-limit
        bpybreaks <- input$barplotYbreaks                # Y-axis breaks
        dotsizein <- input$barplotDotSize
        
        #if (ncol(meta) > 2) {
        #  metacol <- input$BarPlot1MetaCol
        #}
        #else if (ncol(meta) == 2) {
        #  metacol <- colnames(meta)[2]
        #}
        meta_temp <- meta
        rownames(meta_temp) <- meta[,1]
        meta_temp <- meta_temp[,metacol,drop = F]
        meta_temp[,metacol] <- as.factor(meta_temp[,metacol])
        
        if (length(input$GeneListTableBarPlot_rows_selected) > 0){
          
          gene <- geneList[input$GeneListTableBarPlot_rows_selected, 1]
          expr_gene <- data.frame(row.names = names(expr[gene,]),
                                  GeneName = unname(unlist(expr[gene,])))
          #colnames(expr_gene)[1] <- gene
          expr_gene <- merge(expr_gene, meta_temp, by=0)
          colnames(expr_gene)[c(1,3)] <- c("SampleName","Type")
          
          expr_gene2 <- merge(expr_gene,meta, all = T)
          plottitle <- paste(gene,"Average Gene Expression Across",metacol)
          genetitle <- paste(gene,"Average Expression")
          
          if (logchoice == T) {
            expr_gene2[,"GeneName"] <- log2(expr_gene2[,"GeneName"] + 1)
            plottitle <- paste(gene,"Average Gene Expression (Log2) Across",metacol)
            genetitle <- paste(gene,"Average Expression (Log2)")
          }
          
          
          se <- function(x) sd(x)/sqrt(length(x))
          expr_gene_stats <- expr_gene2 %>%
            group_by(Type) %>%
            summarise_at("GeneName",funs(mean,sd,se))
          
          y_min <- 0
          y_max <- round(max((expr_gene2[,2]) + 1))
          if (bpylim != "") {
            y_min <- as.numeric(gsub(" ","",strsplit(bpylim,",")[[1]][1]))
            y_max <- as.numeric(gsub(" ","",strsplit(bpylim,",")[[1]][2]))
          }
          
          if (input$barplotXaxOrder == "Descending"){
            barp <- ggplot(expr_gene_stats, aes(reorder(Type,-mean, na.rm = TRUE), mean, fill=Type))
          }
          if (input$barplotXaxOrder == "Ascending"){
            barp <- ggplot(expr_gene_stats, aes(reorder(Type, mean, na.rm = TRUE), mean, fill=Type))
          }
          if (bar_type != "None") {
            if (bar_type == "Standard Deviation") {
              barp <- barp + geom_errorbar(aes(ymin=0,ymax = mean+sd), size = 0.75, width = 0.5)
            }
            else if (bar_type == "Standard Error") {
              barp <- barp + geom_errorbar(aes(ymin=0,ymax = mean+se), size = 0.75, width = 0.5)
            }
          }
          barp <- barp + geom_bar(stat = "identity",
                                  width=0.75,
                                  size = 0.75,
                                  color="black",
                                  show.legend = FALSE) +
            theme_minimal()
          
          barp <- barp + labs(x = metacol,
                              y = genetitle,
                              title = plottitle)
          barp <- barp + theme(axis.text.x = element_text(angle = axis_orient, hjust = hjust_orient),
                               axis.text = element_text(size = axis_font),
                               axis.title = element_text(size = axis_font),
                               plot.title = element_text(size = title_font))
          if (bpylim == "") {
            if (is.na(bpybreaks)) {
              barp <- barp + scale_y_continuous(expand = c(0, 0))
            }
            else if (!is.na(bpybreaks)) {
              barp <- barp + scale_y_continuous(limits=c(0,y_max),expand = c(0, 0),breaks=seq(0,y_max,bpybreaks),oob = rescale_none)
            }
          }
          else if (bpylim != "") {
            if (is.na(bpybreaks)) {
              barp <- barp + scale_y_continuous(limits=c(y_min,y_max),expand = c(0, 0),oob = rescale_none)
            }
            else if (!is.na(bpybreaks)) {
              barp <- barp + scale_y_continuous(limits=c(y_min,y_max),expand = c(0, 0),breaks=seq(y_min,y_max,bpybreaks),oob = rescale_none)
            }
          }
          
          if (input$barplotsampledots == T) {
            barp <- barp + geom_dotplot(data = expr_gene2, aes(x=Type,y=GeneName),
                                        binaxis='y', stackdir='center',
                                        stackratio=1, dotsize=dotsizein, fill = "black")
          }
          
          if (colorin != "") {
            colorin <- strsplit(colorin," ")[[1]]
            if (length(colorin) == 1) {
              colorin <- rep(colorin,length(unique(meta[,metacol])))
            }
            barp <- barp + scale_fill_manual(values=colorin)
          }
          barp <- barp + coord_cartesian(clip = "off")
          
          
          
          barp
        }
      })
      
      output$barplot <- renderPlot({
        
        bpplot <- barplot_react()
        bpplot
        
      })
      
      #render up regulated pathway enrichment plot
      output$UpRegPathway1 <- renderPlot({
        title_font <- input$PathwayTitleSize
        axis_font <- input$PathwayAxisSize
        adjp <- input$pathpval
        FC <- input$pathFC
        if (ncol(meta) > 2) {
          metacol <- input$EnrichRPathMetaCol
        }
        else if (ncol(meta) == 2) {
          metacol <- colnames(meta)[2]
        }
        A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
        B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
        mat <- expr[,c(A,B)]
        mat <- log2(mat + 1.0)
        groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
        designA <- model.matrix(~0 + groupAOther)
        fit <- lmFit(mat, design = designA)
        contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        options(digits = 4)
        top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
        genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
        dbs <- listEnrichrDbs()
        enrichRLive <- TRUE
        if (is.null(dbs)) {
          enrichRLive <- FALSE
        }
        dbs <- input$SelectedPathway
        enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
        plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value") +
          theme(axis.title = element_text(size = axis_font),
                plot.title = element_text(size = title_font))
      })
      
      #render up regulated pathway enrichment plot
      output$DnRegPathway1 <- renderPlot({
        title_font <- input$PathwayTitleSize
        axis_font <- input$PathwayAxisSize
        adjp <- input$pathpval
        FC <- input$pathFC
        if (ncol(meta) > 2) {
          metacol <- input$EnrichRPathMetaCol
        }
        else if (ncol(meta) == 2) {
          metacol <- colnames(meta)[2]
        }
        A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
        B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
        mat <- expr[,c(A,B)]
        mat <- log2(mat + 1.0)
        groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
        designA <- model.matrix(~0 + groupAOther)
        fit <- lmFit(mat, design = designA)
        contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        options(digits = 4)
        top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
        genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
        dbs <- listEnrichrDbs()
        enrichRLive <- TRUE
        if (is.null(dbs)) {
          enrichRLive <- FALSE
        }
        dbs <- input$SelectedPathway
        enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
        plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value") +
          theme(axis.title = element_text(size = axis_font),
                plot.title = element_text(size = title_font))
      })
      
      Volcano3_react <- reactive({
        req(topgenereact())
        #title_font <- input$VolMATitleSize
        axis_font <- input$VolMAAxisSize
        tick_font <- input$VolMATickSize
        anno_font <- input$VolMAAnnoSize
        anno_face <- input$VolMAAnnoFont
        top2 <- topgenereact()
        pvalChoice <- input$PorAdjPval
        upRed <- "lightcoral"
        dnBlue <- "cadetblue3"
        mdGray <- "gray70"
        log_fc <- input$fc_cutoff
        pval_cut <- input$p_cutoff
        top_hits <- input$top_x
        #add color categories based on FC and pval
        top2['threshold'] <- "none"
        if (pvalChoice == "-log10(P.value)") {
          top2[which(top2$logFC > abs(log_fc) & top2$P.Value < pval_cut), "threshold"] <- "up"
          top2[which(top2$logFC < -abs(log_fc) & top2$P.Value < pval_cut), "threshold"] <- "down"
          #select number of top hits based on input
          top_hits_up <- top2[head(which(top2$logFC > abs(log_fc) & top2$P.Value < pval_cut), n = top_hits),]
          top_hits_dn <- top2[head(which(top2$logFC < -abs(log_fc) & top2$P.Value < pval_cut), n = top_hits),]
          #create plot
          x <- ggplot(data = top2, aes(x = logFC, y = -log10(P.Value))) +
            geom_point(size = 2, shape = 16) +
            theme_light(base_size = 16)
        } else if (pvalChoice == "-log10(Adjusted P.value)") {
          top2[which(top2$logFC > abs(log_fc) & top2$adj.P.Val < pval_cut), "threshold"] <- "up"
          top2[which(top2$logFC < -abs(log_fc) & top2$adj.P.Val < pval_cut), "threshold"] <- "down"
          #select number of top hits based on input
          top_hits_up <- top2[head(which(top2$logFC > abs(log_fc) & top2$adj.P.Val < pval_cut), n = top_hits),]
          top_hits_dn <- top2[head(which(top2$logFC < -abs(log_fc) & top2$adj.P.Val < pval_cut), n = top_hits),]
          #create plot
          x <- ggplot(data = top2, aes(x = logFC, y = -log10(adj.P.Val))) +
            geom_point(size = 2, shape = 16) +
            theme_light(base_size = 16)
        }
        #select genes to label based on user selection
        genesel.s <- NULL
        genesel.t <- NULL
        genesel.u <- NULL
        genesel.s <- unlist(strsplit(input$userGeneSelec2, " "))
        genesel.t <- unlist(strsplit(input$userGeneSelec2, "\t"))
        genesel.u <- input$userGeneSelec
        genesel.text <- c(genesel.s,genesel.t,genesel.u)
        top2_selec <- top2 %>%
          filter(GeneName %in% genesel.text)
        #colors
        x <- x + aes(color = threshold) +
          scale_color_manual(values = c("up" = upRed,"down" = dnBlue, "none" = mdGray))
        #FC and pval lines
        x <- x + geom_vline(xintercept = c(-abs(log_fc),abs(log_fc)), linetype="dashed", color="gray20")
        x <- x + geom_hline(yintercept = -log10(pval_cut), linetype="dashed", color="gray20")
        #label top hits if needed
        if (top_hits > 0) {
          x <- x + geom_text_repel(
            data =  top_hits_up,
            #aes(label = rownames(top_hits_up), fontface = anno_face),
            aes(label = rownames(top_hits_up)),
            size = anno_font,
            color="gray20",
            min.segment.length = 0,
            #nudge_x = 0.2,
            #nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = Inf)
          x <- x + geom_text_repel(
            data =  top_hits_dn,
            #aes(label = rownames(top_hits_dn), fontface = anno_face),
            aes(label = rownames(top_hits_dn)),
            size = anno_font,
            color="gray20",
            min.segment.length = 0,
            #nudge_x = 0.2,
            #nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = Inf)
        }
        x <- x + geom_text_repel(
          data =  top2_selec,
          #aes(label = rownames(top2_selec), fontface = anno_face),
          aes(label = rownames(top2_selec)),
          size = anno_font,
          color="gray20",
          min.segment.length = 0,
          #nudge_x = 0.2,
          #nudge_y=0.2,
          box.padding = unit(0.9, "lines"),
          point.padding = unit(.3+4*0.1, "lines"),
          max.overlaps = Inf)
        #coloring selected points
        if (pvalChoice == "-log10(P.value)") {
          x <- x + geom_point(data = top2_selec,
                              aes(x = logFC, y = -log10(P.Value)),
                              pch = 21,
                              color = "black",
                              size = 2)
          x <- x + geom_point(data = top_hits_dn,
                              aes(x = logFC, y = -log10(P.Value)),
                              pch = 21,
                              color = "black",
                              size = 2)
          x <- x + geom_point(data = top_hits_up,
                              aes(x = logFC, y = -log10(P.Value)),
                              pch = 21,
                              color = "black",
                              size = 2)
        } else if (pvalChoice == "-log10(Adjusted P.value)") {
          x <- x + geom_point(data = top2_selec,
                              aes(x = logFC, y = -log10(adj.P.Val)),
                              pch = 21,
                              color = "black",
                              size = 2)
          x <- x + geom_point(data = top_hits_dn,
                              aes(x = logFC, y = -log10(adj.P.Val)),
                              pch = 21,
                              color = "black",
                              size = 2)
          x <- x + geom_point(data = top_hits_up,
                              aes(x = logFC, y = -log10(adj.P.Val)),
                              pch = 21,
                              color = "black",
                              size = 2)
        }
        #axis parameters
        x <- x + theme(legend.position="none")
        x <- x + theme(axis.text = element_text(size=tick_font),
                       axis.title = element_text(size = axis_font))
        #x <- x + theme(axis.title = element_text(size = axis_font))
        x <- x + labs(x = "log2FC", y = pvalChoice)
        #save(list = ls(), file = "volcano_env.RData", envir = environment())
        x
        
        
      })
      
      #render volcano plot
      output$Volcano3 <- renderPlot({
        plot <- Volcano3_react()
        plot
      })
      
      boxplot2_react <- reactive({
        req(ssGSEAfunc())
        title_font <- input$gseaBoxTitleSize
        axis_font <- input$gseaBoxAxisSize
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        #if (ncol(meta) > 2) {
        #  metacol <- input$GSEAmetaCol
        #}
        #else if (ncol(meta) == 2) {
        #  metacol <- colnames(meta)[2]
        #}
        
        ssgsea <- ssGSEAfunc()
        ssgsea2 <- as.data.frame(t(ssgsea))
        samporder <- meta[,1]
        ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        rownames(ssgsea3) <- samporder
        ssgsea4 <- ssgsea3 %>%
          mutate(type = case_when(
            rownames(ssgsea3) == meta[,1] ~ meta[,metacol],
          ))
        ssgsea4 <- ssgsea4 %>%
          relocate(type)
        ssgsea4$type <- as.factor(ssgsea4$type)
        colnames(ssgsea4)[which(colnames(ssgsea4) == "type")] <- metacol
        if (input$boxplotcompare == "none") {
          p <- ggplot(ssgsea4, aes(factor(ssgsea4[,metacol]), ssgsea4[,2], fill = ssgsea4[,metacol])) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                 title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = ""),
                 fill = metacol) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = axis_font),
                  plot.title = element_text(size = title_font),
                  legend.position = "none")
        } else {
          p <- ggplot(ssgsea4, aes(factor(ssgsea4[,metacol]), ssgsea4[,2], fill = ssgsea4[,metacol])) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                 title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = ""),
                 fill = metacol) +
            theme_bw() +
            stat_compare_means(method = input$boxplotcompare) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = axis_font),
                  plot.title = element_text(size = title_font),
                  legend.position = "none")
        }
        p
        
        
        
      })
      
      #render ssGSEA boxplot
      output$boxplot2 <- renderPlot({
        
        p <- boxplot2_react()
        p
        
        
      })
      
      
      AvgGeneScatter2_react <- reactive({
        
        A_choice <- input$comparisonA2.avg
        B_choice <- input$comparisonB2.avg
        log_choice <- input$AvgExpLogFC
        title_font <- input$AvgExprScatterTitleSize
        axis_font <- input$AvgExprScatterAxisSize
        
        #select genes to label based on user selection
        genesel.s <- NULL
        genesel.t <- NULL
        genesel.u <- NULL
        genesel.s <- unlist(strsplit(input$gsSelection2, " "))
        genesel.t <- unlist(strsplit(input$gsSelection2, "\t"))
        genesel.u <- input$scatterGeneSelec
        genesel.text <- c(genesel.s,genesel.t,genesel.u)
        
        AvgExpr_Table <- AvgExprReact()
        
        logFC_text <- ""
        
        # Log if user designates
        if (log_choice == TRUE) {
          logFC_text <- " (log2 +1)"
        }
        
        # Add group for above or below abline
        AvgExpr_Table <- AvgExpr_Table %>%
          mutate(ColorCol = case_when(
            AvgExpr_Table$AvgExpression_GroupB > AvgExpr_Table$AvgExpression_GroupA ~ paste(B_choice," Expr > ",A_choice,' Expr', sep = ''),
            AvgExpr_Table$AvgExpression_GroupA > AvgExpr_Table$AvgExpression_GroupB ~ paste(A_choice," Expr > ",B_choice,' Expr', sep = '')
          ))
        
        AvgExpr_Table_selec <- AvgExpr_Table %>%
          filter(GeneSymbol %in% genesel.text)
        
        # plot
        p <- ggplot(AvgExpr_Table, aes(x = AvgExpression_GroupA, y = AvgExpression_GroupB,
                                       color = ColorCol)) +
          geom_point() +
          theme_minimal() +
          theme(legend.position="none",
                axis.title = element_text(size = axis_font),
                plot.title = element_text(size = title_font)) +
          xlab(paste("Average Expression: ",A_choice,logFC_text,sep = ""))+
          ylab(paste("Average Expression: ",B_choice,logFC_text,sep = "")) +
          labs(title = paste("Average Expression",logFC_text,": ",A_choice," vs. ",B_choice,sep = "")) +
          scale_color_manual(values = c("cadetblue3","lightcoral"))
        
        p <- p + geom_text_repel(
          data =  AvgExpr_Table_selec,
          aes(label = GeneSymbol),
          size = 6,
          color="gray20",
          #nudge_x = 0.2,
          #nudge_y=0.2,
          box.padding = unit(0.9, "lines"),
          point.padding = unit(.3+4*0.1, "lines"),
          max.overlaps = 50)
        
        #coloring selected points
        p <- p + geom_point(data = AvgExpr_Table_selec,
                            aes(x = AvgExpression_GroupA, y = AvgExpression_GroupB),
                            pch = 21,
                            color = "black",
                            size = 2)
        
        p
        
      })
      
      
      # Render Average Gene Expression comparison plot
      output$AvggeneScatter2 <- renderPlot({
        
        p <- AvgGeneScatter2_react()
        p
        
      })
      
      
      ####----Download Handlers----####
      
      output$dnldssgseaheatmap_df <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_ssGSEA_zScore_Heatmap_Table.txt", sep = '')
        },
        content = function(file) {
          df <- ssgseaheatmap_df_react()
          write.table(df,file, sep = '\t', row.names = F)
        }
      )
      output$dnldssgseaheatmap2_df <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_ssGSEA_RawDifference_Heatmap_Table.txt", sep = '')
        },
        content = function(file) {
          df <- ssgseaheatmap2_df_react()
          write.table(df,file, sep = '\t', row.names = F)
        }
      )
      output$dnldssgseaheatmap3_df <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_ssGSEA_AvgRawDifference_Heatmap_Table.txt", sep = '')
        },
        content = function(file) {
          df <- ssgseaheatmap3_df_react()
          write.table(df,file, sep = '\t', row.names = F)
        }
      )
      
      output$dnldPlotSVG_exprBar <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_AvgExpression_Barplot.svg", sep = '')
        },
        content = function(file) {
          bpplot <- barplot_react()
          ggsave(file,bpplot, width = 10, height = 8)
        }
      )
      
      output$dnldPlotPDF_exprBar <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_AvgExpression_Barplot.pdf", sep = '')
        },
        content = function(file) {
          bpplot <- barplot_react()
          ggsave(file,bpplot, width = 10, height = 8)
        }
      )
      
      output$dnldPlotSVG_ssgseaHeat <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_ssGSEA_zScore_Heatmap.svg", sep = '')
        },
        content = function(file) {
          svg(filename = file, height = input$ssgseaheatmapHeight, width = input$ssgseaheatmapWidth)
          ComplexHeatmap::draw(ssgseaheatmap_react())
          dev.off()
        }
      )
      
      
      output$dnldPlotPDF_ssgseaHeat <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_ssGSEA_zScore_Heatmap.pdf", sep = '')
        },
        content = function(file) {
          pdf(filename = file, height = input$ssgseaheatmapHeight, width = input$ssgseaheatmapWidth)
          ComplexHeatmap::draw(ssgseaheatmap_react())
          dev.off()
        }
      )
      
      
      output$dnldPlotSVG_ssgseaHeat2 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_ssGSEA_RawDiff_Heatmap.svg", sep = '')
        },
        content = function(file) {
          svg(filename = file, height = input$ssgseaheatmapHeight, width = input$ssgseaheatmapWidth)
          ComplexHeatmap::draw(ssgseaheatmap2_react())
          dev.off()
        }
      )
      
      
      output$dnldPlotPDF_ssgseaHeat2 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_ssGSEA_RawDiff_Heatmap.pdf", sep = '')
        },
        content = function(file) {
          pdf(filename = file, height = input$ssgseaheatmapHeight, width = input$ssgseaheatmapWidth)
          ComplexHeatmap::draw(ssgseaheatmap2_react())
          dev.off()
        }
      )
      
      
      
      output$dnldPlotSVG_ssgseaHeat3 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_ssGSEA_Heatmap.svg", sep = '')
        },
        content = function(file) {
          
          heat <- ssgseaheatmap3_react()
          
          ggsave(file,heat, width = input$ssgseaheatmapWidth, height = input$ssgseaheatmapHeight)
        }
      )
      
      output$dnldPlotPDF_ssgseaHeat3 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_ssGSEA_Heatmap.pdf", sep = '')
        },
        content = function(file) {
          
          heat <- ssgseaheatmap3_react()
          ggsave(file,heat, width = 15, height = 20)
        }
      )
      
      
      
      
      output$dnldPlotSVG_heat1 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_MVG_Heatmap.svg", sep = '')
        },
        content = function(file) {
          
          #heat <- MVGheatmap_react()
          
          #ggsave(file,heat, width = 10, height = 20)
          svg(filename = file, height = input$heatmapHeight1, width = input$heatmapWidth1)
          ComplexHeatmap::draw(MVGheatmap_react())
          dev.off()
        }
      )
      
      output$dnldPlotSVG_heat2 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_Custom_Heatmap.svg", sep = '')
        },
        content = function(file) {
          
          #heat <- CustomHeatmap_react()
          #ggsave(file,heat, width = 10, height = 20)
          svg(filename = file, height = input$heatmapHeight2, width = input$heatmapWidth2)
          ComplexHeatmap::draw(CustomHeatmap_react())
          dev.off()
        }
      )
      
      output$dnldPlotSVG_heat3 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_DEGExpr_Heatmap.svg", sep = '')
        },
        content = function(file) {
          
          #heat <- DEGHeatmap_react()
          
          #ggsave(file,heat, width = 10, height = 20)
          svg(filename = file, height = input$heatmapHeight3, width = input$heatmapWidth3)
          ComplexHeatmap::draw(DEGHeatmap_react())
          dev.off()
        }
      )
      
      output$dnldPlotSVG_heat5 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_CustomAvgExpr_Heatmap.svg", sep = '')
        },
        content = function(file) {
          
          heat <- AvgExprHeatmap_react()
          ggplot2::ggsave(file,heat, width = 10, height = 20, unit = "in")
          
        }
      )
      
      output$dnldPlotSVG_scatter <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_2GeneExpression_ScatterPlot.svg", sep = '')
        },
        content = function(file) {
          meta <- meta_react()
          metacol <- metacol_reactVal()
          if (length(colnames(meta)) == 2) {
            metacol <- colnames(meta)[2]
          }
          expr <- expr_react()
          #log if user designates
          if (input$logask == TRUE) {
            expr <- log2(expr + 1)
          }
          #transpose
          expr_t <- as.data.frame(t(expr))
          #reorder rowname to match meta for merging
          samporder <- meta[,1]
          expr_t2 <- as.data.frame(expr_t[samporder,])
          #add type
          expr_t3 <- expr_t2 %>%
            mutate(type = case_when(
              rownames(expr_t2) == meta[,1] ~ meta[,metacol],
            ))
          expr_t3 <- expr_t3 %>%
            relocate(type)
          #user gene input
          gene1.u <- input$scatterG1
          gene2.u <- input$scatterG2
          #get columns and info based off user input
          gene1 <- expr_t3[,gene1.u]
          gene2 <- expr_t3[,gene2.u]
          if (input$logask == TRUE) {
            gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression (log2 +1)", sep = "")
            gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression (log2 +1)", sep = "")
          }
          else if (input$logask == FALSE) {
            gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression", sep = "")
            gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression", sep = "")
          }
          #plot
          p <- ggplot(expr_t3, aes(x = gene1, y = gene2,
                                   color = type)) +
            geom_point() +
            theme_minimal() +
            labs(x = gene1.n, y = gene2.n,
                 title = paste(gene1.n, " vs. ", gene2.n, sep = ''))
          ggsave(file,p, width = 10, height = 8)
        }
      )
      
      output$dnldPlotSVG_exprBox <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_Boxplot.svg", sep = '')
        },
        content = function(file) {
          p <- CohortBPPlot_react()
          ggsave(file,p, width = 10, height = 8)
        }
      )
      
      output$dnldDataExplor_Box_plot_df <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_Boxplot_Data.txt", sep = '')
        },
        content = function(file) {
          df <- CohortBPPlot_df_react()
          write.table(df,file,sep = '\t', row.names = F)
        }
      )
      output$dnldDataExplor_Box_plot_df2 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_Boxplot_Data.txt", sep = '')
        },
        content = function(file) {
          df <- CohortBPPlot_df_react2()
          write.table(df,file,sep = '\t', row.names = F)
        }
      )
      
      output$dnldPlotSVG_vol <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_VolcanoPlot.svg", sep = '')
        },
        content = function(file) {
          
          x <- Volcano3_react()
          PlotH <- input$VolMAdnldHeight
          PlotW <- input$VolMAdnldWidth
          PlotU <- input$VolMAdnldSizeUnits
          ggsave(file,x, width = PlotW, height = PlotH, units = PlotU)
        }
      )
      
      output$dnldPlotSVG_MA <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_MA_Plot.svg", sep = '')
        },
        content = function(file) {
          top2 <- topgenereact()
          #add color categories based on FC and pval
          top2['threshold'] <- "none"
          top2[which(top2$logFC > abs(input$fc_cutoff)), "threshold"] <- "up"
          top2[which(top2$logFC < -abs(input$fc_cutoff)), "threshold"] <- "down"
          upRed <- "lightcoral"
          dnBlue <- "cadetblue3"
          mdGray <- "gray70"
          top2u <- top2[order(top2[,1], decreasing = TRUE),]
          top2d <- top2[order(top2[,1]),]
          top_hits_up <- top2u[head(which(top2u$logFC > abs(input$fc_cutoff)), n = input$top_x),]
          top_hits_dn <- top2d[head(which(top2d$logFC < -abs(input$fc_cutoff)), n = input$top_x),]
          #select genes to label based on user selection
          genesel.s <- NULL
          genesel.t <- NULL
          genesel.u <- NULL
          genesel.s <- unlist(strsplit(input$userGeneSelec2, " "))
          genesel.t <- unlist(strsplit(input$userGeneSelec2, "\t"))
          genesel.u <- input$userGeneSelec
          genesel.text <- c(genesel.s,genesel.t,genesel.u)
          top2_selec <- top2 %>%
            filter(GeneName %in% genesel.text)
          x <- ggplot(data = top2, aes(x = AveExpr, y = logFC)) +
            geom_point(size = 2, shape = 16) +
            theme_light(base_size = 16)
          #colors
          x <- x + aes(color = threshold) +
            scale_color_manual(values = c("up" = upRed,"down" = dnBlue, "none" = mdGray))
          x <- x + geom_hline(yintercept = c(-abs(input$fc_cutoff),abs(input$fc_cutoff)), linetype="dashed", color="gray20")
          x <- x + geom_text_repel(
            data =  top_hits_up,
            aes(label = rownames(top_hits_up)),
            color="gray20",
            size = 6,
            #nudge_x = 0.2,
            #nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = 50)
          x <- x + geom_text_repel(
            data =  top_hits_dn,
            aes(label = rownames(top_hits_dn)),
            color="gray20",
            size = 6,
            #nudge_x = 0.2,
            #nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = 50)
          x <- x + geom_text_repel(
            data =  top2_selec,
            aes(label = rownames(top2_selec)),
            color="gray20",
            size = 6,
            #nudge_x = 0.2,
            #nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = 50)
          #coloring selected points
          x <- x + geom_point(data = top2_selec,
                              aes(x = AveExpr, y = logFC),
                              pch = 21,
                              color = "black",
                              size = 2)
          x <- x + geom_point(data = top_hits_dn,
                              aes(x = AveExpr, y = logFC),
                              pch = 21,
                              color = "black",
                              size = 2)
          x <- x + geom_point(data = top_hits_up,
                              aes(x = AveExpr, y = logFC),
                              pch = 21,
                              color = "black",
                              size = 2)
          x <- x + theme(legend.position="none")
          x <- x + labs(x = "Average Expression", y = "log2FC")
          x <- x + theme(axis.text = element_text(size=18))
          x <- x + theme(axis.title = element_text(size=24))
          ggsave(file,x, width = 10, height = 8)
        }
      )
      
      output$dnldPlotSVG_exprBox2 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_Boxplot.svg", sep = '')
        },
        content = function(file) {
          p <- CohortBPPlot_react2()
          ggsave(file,p, width = 10, height = 8)
        }
      )
      
      output$dnldPlotSVG_AvgScatter <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_AvgExpression_ScatterPlot.svg", sep = '')
        },
        content = function(file) {
          p <- AvgGeneScatter2_react()
          
          ggsave(file,p, width = 10, height = 8)
        }
      )
      
      output$dnldPlotSVG_upPath <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_UpRegulated_Pathways.svg", sep = '')
        },
        content = function(file) {
          adjp <- input$pathpval
          FC <- input$pathFC
          if (ncol(meta) > 2) {
            metacol <- input$EnrichRPathMetaCol
          }
          else if (ncol(meta) == 2) {
            metacol <- colnames(meta)[2]
          }
          A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
          B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
          mat <- expr[,c(A,B)]
          mat <- log2(mat + 1.0)
          groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
          designA <- model.matrix(~0 + groupAOther)
          fit <- lmFit(mat, design = designA)
          contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
          fit2 <- contrasts.fit(fit, contrast.matrix)
          fit2 <- eBayes(fit2)
          options(digits = 4)
          top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
          genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
          dbs <- listEnrichrDbs()
          enrichRLive <- TRUE
          if (is.null(dbs)) {
            enrichRLive <- FALSE
          }
          dbs <- input$SelectedPathway
          enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
          x <- plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
          ggsave(file,x, width = 10, height = 8)
        }
      )
      
      output$dnldPlotSVG_dnPath <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_DownRegulated_Pathways.svg", sep = '')
        },
        content = function(file) {
          adjp <- input$pathpval
          FC <- input$pathFC
          if (ncol(meta) > 2) {
            metacol <- input$EnrichRPathMetaCol
          }
          else if (ncol(meta) == 2) {
            metacol <- colnames(meta)[2]
          }
          A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
          B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
          mat <- expr[,c(A,B)]
          mat <- log2(mat + 1.0)
          groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
          designA <- model.matrix(~0 + groupAOther)
          fit <- lmFit(mat, design = designA)
          contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
          fit2 <- contrasts.fit(fit, contrast.matrix)
          fit2 <- eBayes(fit2)
          options(digits = 4)
          top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
          genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
          dbs <- listEnrichrDbs()
          enrichRLive <- TRUE
          if (is.null(dbs)) {
            enrichRLive <- FALSE
          }
          dbs <- input$SelectedPathway
          enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
          x <- plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
          ggsave(file,x, width = 10, height = 8)
        }
      )
      
      output$dnldPlotSVG_gsea <- downloadHandler(
        filename = function() {
          gmt <- gmt_react()
          GS <- paste(unique(gmt[,1]), collapse = "_")
          paste(gsub(" ","",ProjectName_react()),"_",GS,"_GSEA_EnrichmentPlot.svg", sep = '')
        },
        content = function(file) {
          x <- enrichplot0_react()
          ggsave(file,x, width = 10, height = 8)
        }
      )
      
      output$dnldPlotSVG_gseaHeat <- downloadHandler(
        filename = function() {
          gmt <- gmt_react()
          GS <- paste(unique(gmt[,1]), collapse = "_")
          paste(gsub(" ","",ProjectName_react()),"_",GS,"_GSEA_Heatmap.svg", sep = '')
        },
        content = function(file) {
          
          svg(filename = file, height = input$GSEAheatmapHeight1, width = input$GSEAheatmapWidth1)
          ComplexHeatmap::draw(heatmap0_react())
          dev.off()
        }
      )
      
      
      output$dnldPlotSVG_gseaBox <- downloadHandler(
        filename = function() {
          gs <- gs_react()
          GS <- names(gs)
          paste(gsub(" ","",ProjectName_react()),"_",GS,"_ssGSEA_Boxplot.svg", sep = '')
        },
        content = function(file) {
          x <- boxplot2_react()
          ggsave(file,x, width = 10, height = 8)
        }
      )
      
      output$dnldPlotPDF_heat1 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_MVG_Heatmap.pdf", sep = '')
        },
        content = function(file) {
          
          heat <- MVGheatmap_react()
          
          ggsave(file,heat, width = 10, height = 20)
        }
      )
      
      output$dnldPlotPDF_heat2 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_Custom_Heatmap.pdf", sep = '')
        },
        content = function(file) {
          
          heat <- CustomHeatmap_react()
          ggsave(file,heat, width = 10, height = 20)
          
        }
      )
      
      output$dnldPlotPDF_heat3 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_DEGExpr_Heatmap.pdf", sep = '')
        },
        content = function(file) {
          
          heat <- DEGHeatmap_react()
          
          ggsave(file,heat, width = 10, height = 20)
        }
      )
      
      output$dnldPlotPDF_heat5 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_CustomAvgExpr_Heatmap.pdf", sep = '')
        },
        content = function(file) {
          
          heat <- AvgExprHeatmap_react()
          ggsave(file,heat, width = 10, height = 20)
          
        }
      )
      
      output$dnldPlotPDF_scatter <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_2GeneExpression_ScatterPlot.pdf", sep = '')
        },
        content = function(file) {
          meta <- meta_react()
          metacol <- metacol_reactVal()
          if (length(colnames(meta)) == 2) {
            metacol <- colnames(meta)[2]
          }
          expr <- expr_react()
          #log if user designates
          if (input$logask == TRUE) {
            expr <- log2(expr + 1)
          }
          #transpose
          expr_t <- as.data.frame(t(expr))
          #reorder rowname to match meta for merging
          samporder <- meta[,1]
          expr_t2 <- as.data.frame(expr_t[samporder,])
          #add type
          expr_t3 <- expr_t2 %>%
            mutate(type = case_when(
              rownames(expr_t2) == meta[,1] ~ meta[,metacol],
            ))
          expr_t3 <- expr_t3 %>%
            relocate(type)
          #user gene input
          gene1.u <- input$scatterG1
          gene2.u <- input$scatterG2
          #get columns and info based off user input
          gene1 <- expr_t3[,gene1.u]
          gene2 <- expr_t3[,gene2.u]
          if (input$logask == TRUE) {
            gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression (log2 +1)", sep = "")
            gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression (log2 +1)", sep = "")
          }
          else if (input$logask == FALSE) {
            gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression", sep = "")
            gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression", sep = "")
          }
          #plot
          p <- ggplot(expr_t3, aes(x = gene1, y = gene2,
                                   color = type)) +
            geom_point() +
            theme_minimal() +
            labs(x = gene1.n, y = gene2.n,
                 title = paste(gene1.n, " vs. ", gene2.n, sep = ''))
          ggsave(file,p, width = 10, height = 8)
        }
      )
      
      output$dnldPlotPDF_exprBox <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_Boxplot.pdf", sep = '')
        },
        content = function(file) {
          p <- CohortBPPlot_react()
          ggsave(file,p, width = 10, height = 8)
        }
      )
      
      output$dnldPlotPDF_vol <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_VolcanoPlot.pdf", sep = '')
        },
        content = function(file) {
          
          x <- Volcano3_react()
          PlotH <- input$VolMAdnldHeight
          PlotW <- input$VolMAdnldWidth
          PlotU <- input$VolMAdnldSizeUnits
          ggsave(file,x, width = PlotW, height = PlotH, units = PlotU)
        }
      )
      
      output$dnldPlotPDF_MA <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_MA_Plot.pdf", sep = '')
        },
        content = function(file) {
          top2 <- topgenereact()
          #add color categories based on FC and pval
          top2['threshold'] <- "none"
          top2[which(top2$logFC > abs(input$fc_cutoff)), "threshold"] <- "up"
          top2[which(top2$logFC < -abs(input$fc_cutoff)), "threshold"] <- "down"
          upRed <- "lightcoral"
          dnBlue <- "cadetblue3"
          mdGray <- "gray70"
          top2u <- top2[order(top2[,1], decreasing = TRUE),]
          top2d <- top2[order(top2[,1]),]
          top_hits_up <- top2u[head(which(top2u$logFC > abs(input$fc_cutoff)), n = input$top_x),]
          top_hits_dn <- top2d[head(which(top2d$logFC < -abs(input$fc_cutoff)), n = input$top_x),]
          #select genes to label based on user selection
          genesel.s <- NULL
          genesel.t <- NULL
          genesel.u <- NULL
          genesel.s <- unlist(strsplit(input$userGeneSelec2, " "))
          genesel.t <- unlist(strsplit(input$userGeneSelec2, "\t"))
          genesel.u <- input$userGeneSelec
          genesel.text <- c(genesel.s,genesel.t,genesel.u)
          top2_selec <- top2 %>%
            filter(GeneName %in% genesel.text)
          x <- ggplot(data = top2, aes(x = AveExpr, y = logFC)) +
            geom_point(size = 2, shape = 16) +
            theme_light(base_size = 16)
          #colors
          x <- x + aes(color = threshold) +
            scale_color_manual(values = c("up" = upRed,"down" = dnBlue, "none" = mdGray))
          x <- x + geom_hline(yintercept = c(-abs(input$fc_cutoff),abs(input$fc_cutoff)), linetype="dashed", color="gray20")
          x <- x + geom_text_repel(
            data =  top_hits_up,
            aes(label = rownames(top_hits_up)),
            color="gray20",
            size = 6,
            #nudge_x = 0.2,
            #nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = 50)
          x <- x + geom_text_repel(
            data =  top_hits_dn,
            aes(label = rownames(top_hits_dn)),
            color="gray20",
            size = 6,
            #nudge_x = 0.2,
            #nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = 50)
          x <- x + geom_text_repel(
            data =  top2_selec,
            aes(label = rownames(top2_selec)),
            color="gray20",
            size = 6,
            #nudge_x = 0.2,
            #nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = 50)
          #coloring selected points
          x <- x + geom_point(data = top2_selec,
                              aes(x = AveExpr, y = logFC),
                              pch = 21,
                              color = "black",
                              size = 2)
          x <- x + geom_point(data = top_hits_dn,
                              aes(x = AveExpr, y = logFC),
                              pch = 21,
                              color = "black",
                              size = 2)
          x <- x + geom_point(data = top_hits_up,
                              aes(x = AveExpr, y = logFC),
                              pch = 21,
                              color = "black",
                              size = 2)
          x <- x + theme(legend.position="none")
          x <- x + labs(x = "Average Expression", y = "log2FC")
          x <- x + theme(axis.text = element_text(size=18))
          x <- x + theme(axis.title = element_text(size=24))
          ggsave(file,x, width = 10, height = 8)
        }
      )
      
      output$dnldPlotPDF_exprBox2 <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_Boxplot.pdf", sep = '')
        },
        content = function(file) {
          p <- CohortBPPlot_react2()
          ggsave(file,p, width = 10, height = 8)
        }
      )
      
      output$dnldPlotPDF_AvgScatter <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_AvgExpression_ScatterPlot.pdf", sep = '')
        },
        content = function(file) {
          p <- AvgGeneScatter2_react()
          
          ggsave(file,p, width = 10, height = 8)
        }
      )
      
      output$dnldPlotPDF_upPath <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_UpRegulated_Pathways.pdf", sep = '')
        },
        content = function(file) {
          adjp <- input$pathpval
          FC <- input$pathFC
          if (ncol(meta) > 2) {
            metacol <- input$EnrichRPathMetaCol
          }
          else if (ncol(meta) == 2) {
            metacol <- colnames(meta)[2]
          }
          A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
          B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
          mat <- expr[,c(A,B)]
          mat <- log2(mat + 1.0)
          groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
          designA <- model.matrix(~0 + groupAOther)
          fit <- lmFit(mat, design = designA)
          contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
          fit2 <- contrasts.fit(fit, contrast.matrix)
          fit2 <- eBayes(fit2)
          options(digits = 4)
          top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
          genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
          dbs <- listEnrichrDbs()
          enrichRLive <- TRUE
          if (is.null(dbs)) {
            enrichRLive <- FALSE
          }
          dbs <- input$SelectedPathway
          enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
          x <- plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
          ggsave(file,x, width = 15, height = 8)
        }
      )
      
      output$dnldPlotPDF_dnPath <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_DownRegulated_Pathways.pdf", sep = '')
        },
        content = function(file) {
          adjp <- input$pathpval
          FC <- input$pathFC
          if (ncol(meta) > 2) {
            metacol <- input$EnrichRPathMetaCol
          }
          else if (ncol(meta) == 2) {
            metacol <- colnames(meta)[2]
          }
          A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
          B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
          mat <- expr[,c(A,B)]
          mat <- log2(mat + 1.0)
          groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
          designA <- model.matrix(~0 + groupAOther)
          fit <- lmFit(mat, design = designA)
          contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
          fit2 <- contrasts.fit(fit, contrast.matrix)
          fit2 <- eBayes(fit2)
          options(digits = 4)
          top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
          genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
          dbs <- listEnrichrDbs()
          enrichRLive <- TRUE
          if (is.null(dbs)) {
            enrichRLive <- FALSE
          }
          dbs <- input$SelectedPathway
          enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
          x <- plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
          ggsave(file,x, width = 15, height = 8)
        }
      )
      
      output$dnldPlotPDF_gsea <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_GSEA_EnrichmentPlot.pdf", sep = '')
        },
        content = function(file) {
          x <- enrichplot0_react()
          ggsave(file,x, width = 10, height = 8)
        }
      )
      
      output$dnldPlotPDF_gseaHeat <- downloadHandler(
        filename = function() {
          paste(gsub(" ","",ProjectName_react()),"_GSEA_Heatmap.pdf", sep = '')
        },
        content = function(file) {
          meta <- meta_react()
          metacol <- metacol_reactVal()
          if (length(colnames(meta)) == 2) {
            metacol <- colnames(meta)[2]
          }
          A <- A_raw()
          if (input$tables == 1) {
            if (length(input$msigdbTable_rows_selected) > 0){
              groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
              groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
              res <- datasetInput()
              gsea.df <- as.data.frame(res@result)
              GS <- as.character(msigdb.gsea2()[input$msigdbTable_rows_selected,3])
              genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
              genes2 <- strsplit(genes1,"/")
              genes3 <- as.data.frame(genes2, col.names = "genes")
              gene_symbol <- genes3$genes
              ## convert expression matrix to numeric
              class(A) <- "numeric"
              ## Transforming data
              A <- A[,c(groupA,groupB)]
              exp.mat1 = log2(A + 1) # log
              exp.mat2 = apply(exp.mat1, 1, scale); # z score
              exp.mat3 = apply(exp.mat2, 1, rev); # transpose
              colnames(exp.mat3) = colnames(A) # set the column name
              exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
              exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
              # reassign data
              dataset <- exp.mat5
              ## generate color for pheatmap
              if (abs(min(dataset)) > abs(max(dataset))) {
                dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
              } else {
                dataset[dataset > abs(min(dataset))] = abs(min(dataset))
              }
              meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
              meta2 <- meta2[order(meta2[,2]),]
              type <- meta2[,2]
              meta3 <- as.data.frame(type)
              rownames(meta3) <- meta2[,1]
              zscore_range = 10;
              minimum = -zscore_range;
              maximum = zscore_range;
              bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
              #hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
              #Heatmap color
              color_choice <- input$ColorPalette_gseaHeat
              col_sets <- c("OrRd","PuBu","Greens","YlGnBu","Spectral")
              if (color_choice == "original") {
                HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
                hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
              }
              else if (color_choice %in% col_sets) {
                HeatMap_Colors <- brewer.pal(n = 5, color_choice)
                hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
              }
              else if (color_choice == "Inferno") {
                hmcols <- inferno(500)
              }
              else if (color_choice == "Viridis") {
                hmcols <- viridis(500)
              }
              else if (color_choice == "Plasma") {
                hmcols <- plasma(500)
              }
              x <- pheatmap(dataset,            #data
                            cluster_cols = F,    #cluster columns - NO
                            cluster_row = F,     #cluster rows - YES
                            fontsize_col = input$heatmapFont1.c,   #column fontsize
                            fontsize_row = input$heatmapFont1.r,
                            show_rownames = T,
                            show_colnames = T,
                            color=hmcols,
                            annotation_col = meta3)
            }
          }
          else if (input$tables == 3) {
            if (length(input$tab2table_rows_selected) > 0){
              groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
              groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
              res <- datasetInput()
              gsea.df <- as.data.frame(res@result)
              GS <- as.character(GeneSet2()[input$tab2table_rows_selected,1])
              genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
              genes2 <- strsplit(genes1,"/")
              genes3 <- as.data.frame(genes2, col.names = "genes")
              gene_symbol <- genes3$genes
              ## convert expression matrix to numeric
              class(A) <- "numeric"
              ## Transforming data
              A <- A[,c(groupA,groupB)]
              exp.mat1 = log2(A + 1) # log
              exp.mat2 = apply(exp.mat1, 1, scale); # z score
              exp.mat3 = apply(exp.mat2, 1, rev); # transpose
              colnames(exp.mat3) = colnames(A) # set the column name
              exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
              exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
              # reassign data
              dataset <- exp.mat5
              ## generate color for pheatmap
              if (abs(min(dataset)) > abs(max(dataset))) {
                dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
              } else {
                dataset[dataset > abs(min(dataset))] = abs(min(dataset))
              }
              meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
              meta2 <- meta2[order(meta2[,2]),]
              type <- meta2[,2]
              meta3 <- as.data.frame(type)
              rownames(meta3) <- meta2[,1]
              zscore_range = 10;
              minimum = -zscore_range;
              maximum = zscore_range;
              bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
              #hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
              #Heatmap color
              color_choice <- input$ColorPalette_gseaHeat
              col_sets <- c("OrRd","PuBu","Greens","YlGnBu","Spectral")
              if (color_choice == "original") {
                HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
                hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
              }
              else if (color_choice %in% col_sets) {
                HeatMap_Colors <- brewer.pal(n = 5, color_choice)
                hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
              }
              else if (color_choice == "Inferno") {
                hmcols <- inferno(500)
              }
              else if (color_choice == "Viridis") {
                hmcols <- viridis(500)
              }
              else if (color_choice == "Plasma") {
                hmcols <- plasma(500)
              }
              x <- pheatmap(dataset,            #data
                            cluster_cols = F,    #cluster columns - NO
                            cluster_row = F,     #cluster rows - YES
                            fontsize_col = input$heatmapFont1.c,   #column fontsize
                            fontsize_row = input$heatmapFont1.r,
                            show_rownames = T,
                            show_colnames = T,
                            color=hmcols,
                            annotation_col = meta3)
            }
          }
          else if (input$tables == 5) {
            if (length(input$GStable.u_rows_selected) > 0){
              groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
              groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
              res <- datasetInput()
              gsea.df <- as.data.frame(res@result)
              GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
              genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
              genes2 <- strsplit(genes1,"/")
              genes3 <- as.data.frame(genes2, col.names = "genes")
              gene_symbol <- genes3$genes
              ## convert expression matrix to numeric
              class(A) <- "numeric"
              ## Transforming data
              A <- A[,c(groupA,groupB)]
              exp.mat1 = log2(A + 1) # log
              exp.mat2 = apply(exp.mat1, 1, scale); # z score
              exp.mat3 = apply(exp.mat2, 1, rev); # transpose
              colnames(exp.mat3) = colnames(A) # set the column name
              exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
              exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
              # reassign data
              dataset <- exp.mat5
              ## generate color for pheatmap
              if (abs(min(dataset)) > abs(max(dataset))) {
                dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
              } else {
                dataset[dataset > abs(min(dataset))] = abs(min(dataset))
              }
              meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
              meta2 <- meta2[order(meta2[,2]),]
              type <- meta2[,2]
              meta3 <- as.data.frame(type)
              rownames(meta3) <- meta2[,1]
              zscore_range = 10;
              minimum = -zscore_range;
              maximum = zscore_range;
              bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
              #hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
              #Heatmap color
              color_choice <- input$ColorPalette_gseaHeat
              col_sets <- c("OrRd","PuBu","Greens","YlGnBu","Spectral")
              if (color_choice == "original") {
                HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
                hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
              }
              else if (color_choice %in% col_sets) {
                HeatMap_Colors <- brewer.pal(n = 5, color_choice)
                hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
              }
              else if (color_choice == "Inferno") {
                hmcols <- inferno(500)
              }
              else if (color_choice == "Viridis") {
                hmcols <- viridis(500)
              }
              else if (color_choice == "Plasma") {
                hmcols <- plasma(500)
              }
              x <- pheatmap(dataset,            #data
                            cluster_cols = F,    #cluster columns - NO
                            cluster_row = F,     #cluster rows - YES
                            fontsize_col = input$heatmapFont1.c,   #column fontsize
                            fontsize_row = input$heatmapFont1.r,
                            show_rownames = T,
                            show_colnames = T,
                            color=hmcols,
                            annotation_col = meta3)
            }
          }
          ggsave(file,x, width = 10, height = 20)
        }
      )
      
      output$dnldPlotPDF_gseaBox <- downloadHandler(
        filename = function() {
          gs <- gs_react()
          GS <- names(gs)
          paste(gsub(" ","",ProjectName_react()),"_",GS,"_ssGSEA_Boxplot.pdf", sep = '')
        },
        content = function(file) {
          x <- boxplot2_react()
          ggsave(file,x, width = 10, height = 8)
        }
      )
      
      
      output$downloadClusters <- downloadHandler(
        filename = function() {
          paste(ProjectName_react(),"_ClusterResults.tsv", sep = '')
        },
        content = function(file) {
          top_probes <- input$NumFeatures
          expr <- expr_react()
          exp <- expr
          col_labels <- colnames(expr)
          #isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
          #exp <- expr[isexpr,]
          mad <- NULL
          var <- NULL
          cv <- NULL
          var_type <- input$VarianceMeasure
          if (var_type == "MAD"){
            mad <- apply(log2(exp + 1), 1, mad)
            mad <- sort(mad, decreasing = T)
            mad <- head(mad, n = (top_probes +1))
            out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
            colnames(out) <- c("Gene", "MAD", colnames(exp))
            dataset <- exp[names(mad),]
            variable_gene_list <- names(mad)
          }
          if (var_type == "VAR"){
            var <- apply(log2(exp + 1), 1, var)
            var <- sort(var, decreasing = T)
            var <- head(var, n = (top_probes +1))
            out <- cbind(names(var), var[names(var)], exp[names(var),])
            colnames(out) <- c("Gene", "VAR", colnames(exp))
            dataset <- exp[names(var),]
            variable_gene_list <- names(var)
          }
          if (var_type == "CV"){
            cv <- apply(log2(exp + 1), 1, cv)
            cv <- sort(cv, decreasing = T)
            cv <- head(cv, n = (top_probes +1))
            out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
            colnames(out) <- c("Gene", "VAR", colnames(exp))
            dataset <- exp[names(cv),]
            variable_gene_list <- names(cv)
          }
          dataset <- log2(dataset + 1)
          zdataset <- apply(dataset, 1, scale)
          zdataset <- apply(zdataset, 1, rev)
          colnames(zdataset) <- colnames(dataset)
          dataset <- as.matrix(zdataset)
          dataset[is.na(dataset)] <- 0
          dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
          minimum = -5;
          maximum = 5;
          if (abs(min(dataset)) > abs(max(dataset))) {
            dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
          } else {
            dataset[dataset > abs(min(dataset))] = abs(min(dataset))
          }
          results2 = hclust(dist(t(dataset)), method = input$ClusteringMethod)
          m = sort(cutree(results2, k=input$NumClusters))
          output = cbind(colnames(m), as.matrix(m))
          colnames(output) = c("Cluster")
          output <- as.data.frame(output)
          output$SampleName <- rownames(output)
          output <- output %>%
            relocate(SampleName)
          write_tsv(as.data.frame(output),file)
        }
      )
      
      #download gene scatter plot expression data
      output$geneScatterDownload <- downloadHandler(
        filename = function() {
          #user gene input
          gene1.u <- input$scatterG1
          gene2.u <- input$scatterG2
          paste(gene1.u, "_vs_", gene2.u, "_Expression.tsv", sep = "")
        },
        content = function(file) {
          meta <- meta_react()
          metacol <- metacol_reactVal()
          if (length(colnames(meta)) == 2) {
            metacol <- colnames(meta)[2]
          }
          expr <- expr_react()
          #transpose
          expr_t <- as.data.frame(t(expr))
          #reorder rowname to match meta for merging
          samporder <- meta[,1]
          expr_t2 <- as.data.frame(expr_t[samporder,])
          #add type
          expr_t3 <- expr_t2 %>%
            mutate(type = case_when(
              rownames(expr_t2) == meta[,1] ~ meta[,metacol],
            ))
          expr_t3 <- expr_t3 %>%
            relocate(type)
          #user gene input
          gene1.u <- input$scatterG1
          gene2.u <- input$scatterG2
          #get columns and info based off user input
          gene1 <- expr_t3[,gene1.u]
          gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression", sep = "")
          gene2 <- expr_t3[,gene2.u]
          gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression", sep = "")
          #make table
          Sample <- rownames(expr_t3)
          Type <- expr_t3[,'type']
          gene1col <- expr_t3[,gene1.u]
          gene2col <- expr_t3[,gene2.u]
          scatterTab <- data.frame(Sample, Type, gene1col, gene2col)
          colnames(scatterTab)[c(3,4)] <- c(gene1.n, gene2.n)
          write_tsv(scatterTab, file)
        }
      )
      
      output$AvggeneScatterDownload <- downloadHandler(
        filename = function() {
          A_choice <- gsub(" ","",input$comparisonA2.avg)
          B_choice <- gsub(" ","",input$comparisonB2.avg)
          log_choice <- input$AvgExpLogFC
          paste("AvgExpr_",A_choice,"_vs_",B_choice,".tsv",sep = "")
        },
        content = function(file) {
          AvgExpr_Table <- AvgGeneScatterTable_react()
          write.table(AvgExpr_Table,file, sep = '\t', row.names = F)
        }
      )
      
      #download ssGSEA score table
      output$ssGSEAdownload <- downloadHandler(
        filename = function() {
          GS <- gs_react()
          paste('ssGSEAscore_',names(GS),'.tsv',sep = '')
          #if (input$tables == 1) {
          #  paste('ssGSEAscore_',names(gs()[(msigdb.gsea2()[input$msigdbTable_rows_selected,3])]),'.tsv',sep = '')
          #}
          #else if (input$tables == 3) {
          #  paste('ssGSEAscore_',names(gs2()[(GeneSet2()[input$tab2table_rows_selected,1])]),'.tsv',sep = '')
          #}
          #else if (input$tables == 5) {
          #  paste('ssGSEAscore_',user_gs_mirror()[input$GStable.u_rows_selected,1],'.tsv',sep = '')
          #  #paste('ssGSEAscore_',names(RDataListGen()[(GStable.ubg()[input$GStable.u_rows_selected,1])]),'.tsv',sep = '')
          #}
        },
        content = function(file) {
          meta <- meta_react()
          metacol <- metacol_reactVal()
          GS <- gs_react()
          if (length(colnames(meta)) == 2) {
            metacol <- colnames(meta)[2]
          }
          #if (input$tables == 1) {
          #  GS <- gs()[(msigdb.gsea2()[input$msigdbTable_rows_selected,3])]
          #}
          #else if (input$tables == 3) {
          #  GS <- gs2()[(GeneSet2()[input$tab2table_rows_selected,1])]
          #}
          #else if (input$tables == 5) {
          #  gmt <- GStable.ubg()
          #  colnames(gmt) <- c("term","gene")
          #  gs_name <- user_gs_mirror()[input$GStable.u_rows_selected,1]
          #  gmt_sub <- gmt[which(gmt$term == gs_name),]
          #  GS <- list()
          #  for (i in unique(gmt_sub[,1])){
          #    GS[[i]] <- gmt_sub[gmt_sub[,1] == i,]$gene
          #  }
          #  #GS <- RDataListGen()[geneset_names]()[(user_gs_mirror()[input$GStable.u_rows_selected,1])]
          #}
          ssgsea <- ssGSEAfunc()
          ssgsea2 <- as.data.frame(t(ssgsea))
          samporder <- meta[,1]
          ssgsea3 <- as.data.frame(ssgsea2[samporder,])
          colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
          rownames(ssgsea3) <- samporder
          ssgsea4 <- ssgsea3 %>%
            mutate(type = case_when(
              rownames(ssgsea3) == meta[,1] ~ meta[,metacol],
            ))
          ssgsea4 <- ssgsea4 %>%
            relocate(type)
          ssgsea4$sample <- rownames(ssgsea4)
          ssgsea5 <- ssgsea4[,c(3,1,2)]
          rownames(ssgsea5) <- 1:nrow(ssgsea5)
          write_tsv(ssgsea5, file)
        }
      )
      
      observe({
        req(input$volcanoCompChoice2)
        req(input$top_x2)
        req(metacol_reactVal())
        req(input$UpDnChoice)
        
        if (input$volcanoCompChoice2 == "Limma: Two groups") {
          feature_labs <- paste0("_",input$comparisonA2.DEG,"_",input$comparisonB2.DEG)
        } else if (input$volcanoCompChoice2 == "Limma: One group") {
          feature_labs <- paste0("_",input$comparisonA2.DEG_one)
        } else if (input$volcanoCompChoice2 == "DESeq2") {
          feature_labs <- paste0("_",input$DESeqDesignColRef_tab)
        }
        
        if (input$UpDnChoice == "UpAndDown_Regulated") {
          GeneSetName_Up <- paste0("Top",input$top_x2,"Up_",ProjectName_react(),"_",metacol_reactVal(),feature_labs)
          GeneSetName_Dn <- paste0("Top",input$top_x2,"Down_",ProjectName_react(),"_",metacol_reactVal(),feature_labs)
          updateTextInput(session,"DEGGeneSetName1", label = "Upregulated Gene Set Name:",value = GeneSetName_Up)
          updateTextInput(session,"DEGGeneSetName2", label = "Downregulated Gene Set Name:",value = GeneSetName_Dn)
        } else if (input$UpDnChoice == "Up_Regulated") {
          GeneSetName_Up <- paste0("Top",input$top_x2,"Up_",ProjectName_react(),"_",metacol_reactVal(),feature_labs)
          updateTextInput(session,"DEGGeneSetName1", label = "Upregulated Gene Set Name:",value = GeneSetName_Up)
        } else if (input$UpDnChoice == "Down_Regulated") {
          GeneSetName_Dn <- paste0("Top",input$top_x2,"Down_",ProjectName_react(),"_",metacol_reactVal(),feature_labs)
          updateTextInput(session,"DEGGeneSetName1", label = "Downregulated Gene Set Name:",value = GeneSetName_Dn)
        }
        
      })
      
      observe({
        req(input$volcanoCompChoice2)
        req(input$top_x2)
        req(metacol_reactVal())
        req(input$UpDnChoice)
        
        if (input$volcanoCompChoice2 == "Limma: Two groups") {
          feature_labs <- paste0("_",input$comparisonA2.DEG,"_",input$comparisonB2.DEG)
        } else if (input$volcanoCompChoice2 == "Limma: One group") {
          feature_labs <- paste0("_",input$comparisonA2.DEG_one)
        } else if (input$volcanoCompChoice2 == "DESeq2") {
          feature_labs <- paste0("_",input$DESeqDesignColRef_tab)
        }
        
        if (input$UpDnChoice == "UpAndDown_Regulated") {
          NewFileName <- paste0("Top",input$top_x2,"UpAndDown_",ProjectName_react(),"_",metacol_reactVal(),feature_labs,"_Geneset")
        } else if (input$UpDnChoice == "Up_Regulated") {
          NewFileName <- paste0("Top",input$top_x2,"Up_",ProjectName_react(),"_",metacol_reactVal(),feature_labs,"_Geneset")
        } else if (input$UpDnChoice == "Down_Regulated") {
          NewFileName <- paste0("Top",input$top_x2,"Down_",ProjectName_react(),"_",metacol_reactVal(),feature_labs,"_Geneset")
        }
        updateTextInput(session,"DEGfileName",value = NewFileName)
      })
      
      #render download button for DEG GMT
      output$DEGgmtDownload <- downloadHandler(
        filename = function() {
          paste(input$DEGfileName,".gmt",sep = "")
        },
        content = function(file) {
          top1 <- DEGtable1_react()
          
          if (input$UpDnChoice == "UpAndDown_Regulated") {
            GeneSetName_Up <- input$DEGGeneSetName1
            GeneSetName_Dn <- input$DEGGeneSetName2
          } else if (input$UpDnChoice == "Up_Regulated") {
            GeneSetName_Up <- input$DEGGeneSetName1
          } else if (input$UpDnChoice == "Down_Regulated") {
            GeneSetName_Dn <- input$DEGGeneSetName1
          }
          
          description <- paste0("FCCutoff_",abs(input$fc_cutoff2),"_PvalCutoff_",input$p_cutoff2)
          
          if (input$UpDnChoice == "UpAndDown_Regulated"){
            GeneSetName_Up <- ifelse(is.null(GeneSetName_Up),paste0("Top",input$top_x2,"_UpReg"),GeneSetName_Up)
            GeneSetName_Dn <- ifelse(is.null(GeneSetName_Dn),paste0("Top",input$top_x2,"_DownReg"),GeneSetName_Dn)
            if (input$fc_cutoff2 != 0) {
              genes_up <- rownames(top1)[which(top1$logFC > abs(input$fc_cutoff2) & top1$P.Value < input$p_cutoff2)]
              genes_dn <- rownames(top1)[which(top1$logFC > -abs(input$fc_cutoff2) & top1$P.Value < input$p_cutoff2)]
            }  else {
              genes_up <- rownames(top1)[which(top1$logFC > input$fc_cutoff2 & top1$P.Value < input$p_cutoff2)]
              genes_dn <- rownames(top1)[which(top1$logFC < input$fc_cutoff2 & top1$P.Value < input$p_cutoff2)]
            }
            genes_up_h <- t(as.data.frame(head(genes_up, n=input$top_x2)))
            genes_dn_h <- t(as.data.frame(head(genes_dn, n=input$top_x2)))
            genes_up_h_gmt <- data.frame(name = GeneSetName_Up,
                                         description = description,
                                         genes_up_h)
            genes_dn_h_gmt <- data.frame(name = GeneSetName_Dn,
                                         description = description,
                                         genes_dn_h)
            genes_h_gmt <- rbindlist(list(genes_up_h_gmt,genes_dn_h_gmt), fill = T)
          }
          else if (input$UpDnChoice == "Up_Regulated"){
            GeneSetName_Up <- ifelse(is.null(GeneSetName_Up),paste0("Top",input$top_x2,"_UpReg"),GeneSetName_Up)
            genes_up <- rownames(top1)[which(top1$logFC > abs(input$fc_cutoff2) & top1$P.Value < input$p_cutoff2)]
            genes_up_h <- t(as.data.frame(head(genes_up, n=input$top_x2)))
            genes_up_h_gmt <- data.frame(GeneSetName_Up,
                                         description,
                                         genes_up_h)
            genes_h_gmt <- genes_up_h_gmt
          }
          else if (input$UpDnChoice == "Down_Regulated"){
            GeneSetName_Dn <- ifelse(is.null(GeneSetName_Dn),paste0("Top",input$top_x2,"_DownReg"),GeneSetName_Dn)
            if (input$fc_cutoff2 != 0) {
              genes_dn <- rownames(top1)[which(top1$logFC > -abs(input$fc_cutoff2) & top1$P.Value < input$p_cutoff2)]
            }  else {
              genes_dn <- rownames(top1)[which(top1$logFC < input$fc_cutoff2 & top1$P.Value < input$p_cutoff2)]
            }
            #genes_dn <- rownames(top1)[which(top1$logF < -abs(input$fc_cutoff2) & top1$P.Value < input$p_cutoff2)]
            genes_dn_h <- t(as.data.frame(head(genes_dn, n=input$top_x2)))
            genes_dn_h_gmt <- data.frame(GeneSetName_Dn,
                                         description,
                                         genes_dn_h)
            genes_h_gmt <- genes_dn_h_gmt
          }
          write_delim(genes_h_gmt, file, delim = '\t', col_names = F)
        }
      )
      
      #render download button for DEG GMT
      output$DEGtsvDownload <- downloadHandler(
        filename = function() {
          paste(input$DEGfileName,".tsv",sep = "")
        },
        content = function(file) {
          top1 <- DEGtable1_react()
          
          if (input$UpDnChoice == "UpAndDown_Regulated") {
            GeneSetName_Up <- input$DEGGeneSetName1
            GeneSetName_Dn <- input$DEGGeneSetName2
          } else if (input$UpDnChoice == "Up_Regulated") {
            GeneSetName_Up <- input$DEGGeneSetName1
          } else if (input$UpDnChoice == "Down_Regulated") {
            GeneSetName_Dn <- input$DEGGeneSetName1
          }
          
          description <- paste0("FCCutoff_",abs(input$fc_cutoff2),"_PvalCutoff_",input$p_cutoff2)
          
          if (input$UpDnChoice == "UpAndDown_Regulated"){
            GeneSetName_Up <- ifelse(is.null(GeneSetName_Up),paste0("Top",input$top_x2,"_UpReg"),GeneSetName_Up)
            GeneSetName_Dn <- ifelse(is.null(GeneSetName_Dn),paste0("Top",input$top_x2,"_DownReg"),GeneSetName_Dn)
            if (input$fc_cutoff2 != 0) {
              genes_up <- rownames(top1)[which(top1$logFC > abs(input$fc_cutoff2) & top1$P.Value < input$p_cutoff2)]
              genes_dn <- rownames(top1)[which(top1$logFC > -abs(input$fc_cutoff2) & top1$P.Value < input$p_cutoff2)]
            }  else {
              genes_up <- rownames(top1)[which(top1$logFC > input$fc_cutoff2 & top1$P.Value < input$p_cutoff2)]
              genes_dn <- rownames(top1)[which(top1$logFC < input$fc_cutoff2 & top1$P.Value < input$p_cutoff2)]
            }
            genes_up_h <- head(genes_up, n=input$top_x2)
            genes_dn_h <- head(genes_dn, n=input$top_x2)
            genes_up_h_tsv <- data.frame(term = rep(GeneSetName_Up,length(genes_up_h)),
                                         gene = genes_up_h)
            genes_dn_h_tsv <- data.frame(term = rep(GeneSetName_Dn,length(genes_dn_h)),
                                         gene = genes_dn_h)
            genes_h_tsv <- rbindlist(list(genes_up_h_tsv,genes_dn_h_tsv), fill = T)
          }
          else if (input$UpDnChoice == "Up_Regulated"){
            GeneSetName_Up <- ifelse(is.null(GeneSetName_Up),paste0("Top",input$top_x2,"_UpReg"),GeneSetName_Up)
            genes_up <- rownames(top1)[which(top1$logFC > abs(input$fc_cutoff2) & top1$P.Value < input$p_cutoff2)]
            genes_up_h <- head(genes_up, n=input$top_x2)
            genes_up_h_tsv <- data.frame(term = rep(GeneSetName_Up,length(genes_up_h)),
                                         gene = genes_up_h)
            genes_h_tsv <- genes_up_h_tsv
          }
          else if (input$UpDnChoice == "Down_Regulated"){
            GeneSetName_Dn <- ifelse(is.null(GeneSetName_Dn),paste0("Top",input$top_x2,"_DownReg"),GeneSetName_Dn)
            if (input$fc_cutoff2 != 0) {
              genes_dn <- rownames(top1)[which(top1$logFC > -abs(input$fc_cutoff2) & top1$P.Value < input$p_cutoff2)]
            }  else {
              genes_dn <- rownames(top1)[which(top1$logFC < input$fc_cutoff2 & top1$P.Value < input$p_cutoff2)]
            }
            genes_dn_h <- head(genes_dn, n=input$top_x2)
            genes_dn_h_tsv <- data.frame(term = rep(GeneSetName_Dn,length(genes_dn_h)),
                                         gene = genes_dn_h)
            genes_h_tsv <- genes_dn_h_tsv
          }
          write_delim(genes_h_tsv, file, delim = '\t')
        }
      )
      
      #render download button for upreg gene pathways GMT
      output$UpRegPathDownloadgmt <- downloadHandler(
        filename = function() {
          paste("UpReg_",input$SelectedPathway,"_pathway.gmt",sep = "")
        },
        content = function(file) {
          top1 <- topgenereact()
          genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
          dbs <- listEnrichrDbs()
          enrichRLive <- TRUE
          if (is.null(dbs)) {
            enrichRLive <- FALSE
          }
          dbs <- input$SelectedPathway
          enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
          df <- enriched[[1]]
          path.genes <- as.data.frame(df[,'Genes'])
          genes.n <- c()
          for (i in path.genes[,1]) {
            g <- strsplit(i,";")
            for (j in g){
              genes.n <- c(genes.n, j)
            }
          }
          rm(g)
          genes.n <- t(as.data.frame(unique(genes.n)))
          gmt.df <- data.frame(paste("UpReg_",input$SelectedPathway,"_PathwayGenes",sep = ""),
                               input$SelectedPathway,genes.n)
          rownames(gmt.df) <- NULL
          write_delim(gmt.df, file, delim = '\t', col_names = F)
        }
      )
      
      #render download button for dnreg gene pathways GMT
      output$DnRegPathDownloadgmt <- downloadHandler(
        filename = function() {
          paste("DnReg_",input$SelectedPathway,"_pathway.gmt",sep = "")
        },
        content = function(file) {
          top1 <- topgenereact()
          genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
          dbs <- listEnrichrDbs()
          enrichRLive <- TRUE
          if (is.null(dbs)) {
            enrichRLive <- FALSE
          }
          dbs <- input$SelectedPathway
          enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
          df <- enriched[[1]]
          path.genes <- as.data.frame(df[,'Genes'])
          genes.n <- c()
          for (i in path.genes[,1]) {
            g <- strsplit(i,";")
            for (j in g){
              genes.n <- c(genes.n, j)
            }
          }
          rm(g)
          genes <- t(as.data.frame(unique(genes.n)))
          gmt.df <- data.frame(paste("DnReg_",input$SelectedPathway,"_PathwayGenes",sep = ""),
                               input$SelectedPathway,genes.n)
          rownames(gmt.df) <- NULL
          write_delim(gmt.df, file, delim = '\t', col_names = F)
        }
      )
      
      #render download button for upreg gene pathways
      output$UpRegPathDownload <- downloadHandler(
        filename = function() {
          paste("UpReg_",input$SelectedPathway,"_pathway.tsv",sep = "")
        },
        content = function(file) {
          top1 <- topgenereact()
          genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
          dbs <- listEnrichrDbs()
          enrichRLive <- TRUE
          if (is.null(dbs)) {
            enrichRLive <- FALSE
          }
          dbs <- input$SelectedPathway
          enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
          df <- enriched[[1]]
          write_delim(df, file, delim = '\t')
        }
      )
      
      #render download button for dnreg gene pathways
      output$DnRegPathDownload <- downloadHandler(
        filename = function() {
          paste("DnReg_",input$SelectedPathway,"_pathway.tsv",sep = "")
        },
        content = function(file) {
          top1 <- topgenereact()
          genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
          dbs <- listEnrichrDbs()
          enrichRLive <- TRUE
          if (is.null(dbs)) {
            enrichRLive <- FALSE
          }
          dbs <- input$SelectedPathway
          enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
          df <- enriched[[1]]
          write_delim(df, file, delim = '\t')
        }
      )
      
      #render Most Variable Gene download button
      output$MVGdownload <- downloadHandler(
        filename = function() {
          top_probes <- input$NumFeatures
          paste(top_probes,"_Most_Variable_Genes", ".tsv", sep = "")
        },
        content = function(file){
          top_probes <- input$NumFeatures
          expr <- expr_react()
          exp <- expr
          col_labels <- colnames(expr)
          #isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
          #exp <- expr[isexpr,]
          mad <- NULL
          var <- NULL
          cv <- NULL
          var_type <- input$VarianceMeasure
          if (var_type == "MAD"){
            mad <- apply(log2(exp + 1), 1, mad)
            mad <- sort(mad, decreasing = T)
            mad <- head(mad, n = (top_probes +1))
            out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
            colnames(out) <- c("Gene", "MAD", colnames(exp))
            dataset <- exp[names(mad),]
            variable_gene_list <- names(mad)
          }
          if (var_type == "VAR"){
            var <- apply(log2(exp + 1), 1, var)
            var <- sort(var, decreasing = T)
            var <- head(var, n = (top_probes +1))
            out <- cbind(names(var), var[names(var)], exp[names(var),])
            colnames(out) <- c("Gene", "VAR", colnames(exp))
            dataset <- exp[names(var),]
            variable_gene_list <- names(var)
          }
          if (var_type == "CV"){
            cv <- apply(log2(exp + 1), 1, cv)
            cv <- sort(cv, decreasing = T)
            cv <- head(cv, n = (top_probes +1))
            out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
            colnames(out) <- c("Gene", "VAR", colnames(exp))
            dataset <- exp[names(cv),]
            variable_gene_list <- names(cv)
          }
          variable_gene_list <- as.data.frame(variable_gene_list)
          colnames(variable_gene_list)[1] <- "Genes"
          variable_gene_list$Rank <- rownames(variable_gene_list)
          variable_gene_list <- variable_gene_list %>%
            select(Rank, Genes)
          write_delim(variable_gene_list, file, delim = '\t')
        }
      )
      
      #render Most Variable Gene GMT download button
      output$MVGdownloadgmt <- downloadHandler(
        filename = function() {
          top_probes <- input$NumFeatures
          paste(top_probes,"_Most_Variable_Genes", ".gmt", sep = "")
        },
        content = function(file){
          top_probes <- input$NumFeatures
          expr <- expr_react()
          exp <- expr
          col_labels <- colnames(expr)
          #isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
          #exp <- expr[isexpr,]
          mad <- NULL
          var <- NULL
          cv <- NULL
          var_type <- input$VarianceMeasure
          if (var_type == "MAD"){
            mad <- apply(log2(exp + 1), 1, mad)
            mad <- sort(mad, decreasing = T)
            mad <- head(mad, n = (top_probes +1))
            out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
            colnames(out) <- c("Gene", "MAD", colnames(exp))
            dataset <- exp[names(mad),]
            variable_gene_list <- names(mad)
          }
          if (var_type == "VAR"){
            var <- apply(log2(exp + 1), 1, var)
            var <- sort(var, decreasing = T)
            var <- head(var, n = (top_probes +1))
            out <- cbind(names(var), var[names(var)], exp[names(var),])
            colnames(out) <- c("Gene", "VAR", colnames(exp))
            dataset <- exp[names(var),]
            variable_gene_list <- names(var)
          }
          if (var_type == "CV"){
            cv <- apply(log2(exp + 1), 1, cv)
            cv <- sort(cv, decreasing = T)
            cv <- head(cv, n = (top_probes +1))
            out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
            colnames(out) <- c("Gene", "VAR", colnames(exp))
            dataset <- exp[names(cv),]
            variable_gene_list <- names(cv)
          }
          variable_gene_list <- as.data.frame(variable_gene_list)
          colnames(variable_gene_list)[1] <- "Genes"
          genes.n <- c()
          for (i in variable_gene_list[,1]) {
            g <- strsplit(i,";")
            for (j in g){
              genes.n <- c(genes.n, j)
            }
          }
          rm(g)
          genes <- t(as.data.frame(unique(genes.n)))
          gmt.df <- data.frame(paste(top_probes,"_Most_Variable_Genes", sep = ""),
                               "Most_Variable_Genes",genes)
          rownames(gmt.df) <- NULL
          write_delim(gmt.df, file, delim = '\t', col_names = F)
        }
      )
      
      #download button for leading edge genes
      output$LEGdownload <- downloadHandler(
        filename = function() {
          gmt <- gmt_react()
          GS <- unique(gmt[,1])
          paste(GS,"_leading_edge_genes", ".tsv", sep = "")
        },
        content = function(file){
          GeneSymbol <- LeadingEdgeGenes_react()
          write.table(GeneSymbol, file, sep = '\t', row.names = F)
        })
      
      #render DEG table download button
      output$DEGtableDownload <- downloadHandler(
        filename = function() {
          paste("DEG_Table_",Sys.Date(), ".txt", sep = "")
        },
        content = function(file){
          top1 <- DEGtable1_react()
          top1 <- cbind(rownames(top1),top1)
          colnames(top1)[1] <- "Gene"
          write.table(top1,file,sep = '\t', row.names = F)
        }
      )
      
      
      
      #download button for Enrich Sig Table
      output$enrich_sig_download <- downloadHandler(
        filename = function() {
          groupA <- input$comparisonA
          groupB <- input$comparisonB
          paste(input$SigTableChoice,".tsv", sep = "")
        },
        content = function(file){
          gsea.df <- as_tibble(get(ES_Tab_List[which(ES_Tab_List == paste("ES_table",match(input$SigTableChoice, SigNames),sep = ""))]))
          write_delim(gsea.df, file, delim = '\t')
        })
      
      #download button for user enriched signature table
      output$enrich_sig_download.u <- downloadHandler(
        filename = function() {
          groupA <- input$comparisonA
          groupB <- input$comparisonB
          paste("Enrich_Sig_Table_",groupA,"vs",groupB,".tsv", sep = "")
        },
        content = function(file) {
          
          req(GeneratedMSigDBEST())
          gsea.df <- GeneratedMSigDBEST()
          gsea.df <- as_tibble(gsea.df)
          gsea.df <- gsea.df[,-2]
          write.table(gsea.df, file, sep = '\t', row.names = F)
          
          #meta <- meta_react()
          #metacol <- metacol_reactVal()
          #if (length(colnames(meta)) == 2) {
          #  metacol <- colnames(meta)[2]
          #}
          #expr <- expr_react()
          #if (input$tables == 3) {
          #  groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
          #  groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
          #  ##----Signal-to-Noise Calculation----##
          #  A <- A + 0.00000001
          #  P = as.matrix(as.numeric(colnames(A) %in% groupA))
          #  n1 <- sum(P[,1])
          #  M1 <- A %*% P
          #  M1 <- M1/n1
          #  A2 <- A*A
          #  S1 <- A2 %*% P
          #  S1 <- S1/n1 - M1*M1
          #  S1 <- sqrt(abs((n1/(n1-1)) * S1))
          #  P = as.matrix(as.numeric(colnames(A) %in% groupB))
          #  n2 <- sum(P[,1])
          #  M2 <- A %*% P
          #  M2 <- M2/n2
          #  A2 <- A*A
          #  S2 <- A2 %*% P
          #  S2 <- S2/n2 - M2*M2
          #  S2 <- sqrt(abs((n2/(n2-1)) * S2))
          #  rm(A2)
          #  # small sigma "fix" as used in GeneCluster
          #  S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
          #  S2 <- ifelse(S2 == 0, 0.2, S2)
          #  S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
          #  S1 <- ifelse(S1 == 0, 0.2, S1)
          #  M1 <- M1 - M2
          #  rm(M2)
          #  S1 <- S1 + S2
          #  rm(S2)
          #  s2n.matrix <- M1/S1
          #  ##----Reformatting----##
          #  s2n.df <- as.data.frame(s2n.matrix)
          #  s2n.df$GeneID <- rownames(s2n.df)
          #  rownames(s2n.df) <- NULL
          #  data <- dplyr::select(s2n.df, GeneID, V1)
          #  data.gsea <- data$V1
          #  names(data.gsea) <- as.character(data$GeneID)
          #  s2n.matrix.s <- sort(data.gsea, decreasing = T)
          #  ##----GSEA----##
          #  gmt.i <- tab2()
          #  gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
          #  gsea.df <- as.data.frame(gsea.res@result)
          #  write_delim(gsea.df, file, delim = '\t')
          #}
          #else if (input$tables == 5) {
          #  groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
          #  groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
          #  ##----Signal-to-Noise Calculation----##
          #  A <- A + 0.00000001
          #  P = as.matrix(as.numeric(colnames(A) %in% groupA))
          #  n1 <- sum(P[,1])
          #  M1 <- A %*% P
          #  M1 <- M1/n1
          #  A2 <- A*A
          #  S1 <- A2 %*% P
          #  S1 <- S1/n1 - M1*M1
          #  S1 <- sqrt(abs((n1/(n1-1)) * S1))
          #  P = as.matrix(as.numeric(colnames(A) %in% groupB))
          #  n2 <- sum(P[,1])
          #  M2 <- A %*% P
          #  M2 <- M2/n2
          #  A2 <- A*A
          #  S2 <- A2 %*% P
          #  S2 <- S2/n2 - M2*M2
          #  S2 <- sqrt(abs((n2/(n2-1)) * S2))
          #  rm(A2)
          #  # small sigma "fix" as used in GeneCluster
          #  S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
          #  S2 <- ifelse(S2 == 0, 0.2, S2)
          #  S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
          #  S1 <- ifelse(S1 == 0, 0.2, S1)
          #  M1 <- M1 - M2
          #  rm(M2)
          #  S1 <- S1 + S2
          #  rm(S2)
          #  s2n.matrix <- M1/S1
          #  ##----Reformatting----##
          #  s2n.df <- as.data.frame(s2n.matrix)
          #  s2n.df$GeneID <- rownames(s2n.df)
          #  rownames(s2n.df) <- NULL
          #  data <- dplyr::select(s2n.df, GeneID, V1)
          #  data.gsea <- data$V1
          #  names(data.gsea) <- as.character(data$GeneID)
          #  s2n.matrix.s <- sort(data.gsea, decreasing = T)
          #  ##----GSEA----##
          #  gmt.i <- GStable.ubg()
          #  gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
          #  gsea.df <- as.data.frame(gsea.res@result)
          #  write_delim(gsea.df, file, delim = '\t')
          #}
          #else if (input$tables == 1){
          #  if (input$GenerateEST == TRUE) {
          #    gsea.df <- GeneratedMSigDBEST()
          #    gsea.df <- as.data.frame(gsea.df)
          #    write_delim(gsea.df, file, delim = '\t')
          #  }
          #}
        }
      )
      
      
      
      ####----Text----####
      
      
      #NES and Pval output
      output$NESandPval <- renderText({
        req(datasetInput())
        
        if (length(input$GeneSetTable_rows_selected) > 0){
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          pasteOut <- c()
          for (row in seq(length(input$GeneSetTable_rows_selected))) {
            GS = gsea.df[row,1]
            NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
            Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
            NES.o <- paste0("NES: ", NES)
            Pval.o <- paste0("Pvalue: ", Pval)
            if (input$volcanoCompChoice3 == "Two groups") {
              groupA <- input$comparisonA
              groupB <- input$comparisonB
            } else if (input$volcanoCompChoice3 == "One group") {
              groupA <- input$comparisonA_one
              groupB <- paste0("Not_",groupA)
            }
            if (NES > 0){
              UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", groupA , "group.")
            }
            else {
              UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", groupA , "group.")
            }
            out <- paste(GS,NES.o, Pval.o, UpOrDown, sep = '\n')
            pasteOut <- c(pasteOut,out)
          }
          #GS = unique(gmt_react()[,1])
          ##GS = as.character(msigdb.gsea2()[input$msigdbTable_rows_selected,3])
          #NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
          #Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
          #NES.o <- paste0("NES: ", NES)
          #Pval.o <- paste0("Pvalue: ", Pval)
          #if (input$volcanoCompChoice3 == "Limma: Two groups") {
          #  groupA <- input$comparisonA
          #  groupB <- input$comparisonB
          #} else if (input$volcanoCompChoice3 == "Limma: One group") {
          #  groupA <- input$comparisonA_one
          #  groupB <- paste0("Not_",groupA)
          #}
          #if (NES > 0){
          #  UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", groupA , "group.")
          #}
          #else {
          #  UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", groupA , "group.")
          #}
          #paste(NES.o, Pval.o, UpOrDown, sep = '\n')
          paste(pasteOut,collapse = "\n\n")
        }
        else if (length(input$GeneSetTable_rows_selected) == 0){
          paste("Please select gene set from side panel table to begin.", sep = '')
        }
        
        
        
        #if (input$tables == 1){
        #  if (length(input$msigdbTable_rows_selected) > 0){
        #    res <- datasetInput()
        #    gsea.df <- as.data.frame(res@result)
        #    GS = as.character(msigdb.gsea2()[input$msigdbTable_rows_selected,3])
        #    NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
        #    Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
        #    NES.o <- paste0("NES: ", NES)
        #    Pval.o <- paste0("Pvalue: ", Pval)
        #    if (NES > 0){
        #      UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", input$comparisonA , "group.")
        #    }
        #    else {
        #      UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", input$comparisonA , "group.")
        #    }
        #    paste(NES.o, Pval.o, UpOrDown, sep = '\n')
        #  }
        #  else if (length(input$msigdbTable_rows_selected) == 0){
        #    paste("Please select gene set from side panel table to begin.", sep = '')
        #  }
        #}
        #else if (input$tables == 3){
        #  if (length(input$tab2table_rows_selected) > 0){
        #    res <- datasetInput()
        #    gsea.df <- as.data.frame(res@result)
        #    GS = as.character(GeneSet2()[input$tab2table_rows_selected,1])
        #    NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
        #    Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
        #    NES.o <- paste0("NES: ", NES)
        #    Pval.o <- paste0("Pvalue: ", Pval)
        #    if (NES > 0){
        #      UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", input$comparisonA , "group.")
        #    }
        #    else {
        #      UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", input$comparisonA , "group.")
        #    }
        #    paste(NES.o, Pval.o, UpOrDown, sep = '\n')
        #  }
        #  else if (length(input$tab2table_rows_selected) == 0){
        #    paste("Please select gene set from side panel table to begin.", sep = '')
        #  }
        #}
        #else if (input$tables == 5){
        #  if (length(input$GStable.u_rows_selected) > 0){
        #    res <- datasetInput()
        #    gsea.df <- as.data.frame(res@result)
        #    GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        #    NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
        #    Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
        #    NES.o <- paste0("NES: ", NES)
        #    Pval.o <- paste0("Pvalue: ", Pval)
        #    if (NES > 0){
        #      UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", input$comparisonA , "group.")
        #    }
        #    else {
        #      UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", input$comparisonA , "group.")
        #    }
        #    paste(NES.o, Pval.o, UpOrDown, sep = '\n')
        #  }
        #  else if (length(input$GStable.u_rows_selected) == 0){
        #    paste("Please select gene set from side panel table to begin.", sep = '')
        #  }
        #}
      })
      
      output$GenesAboveCutoff1 <- renderText({
        
        req(topgenereact2())
        req(input$DEGcolHeat)
        FC_cutoff <- input$fc_cutoff_h              #FC cutoff for top gene selection
        P_cutoff <- input$p_cutoff_h                #P-value cutoff for top gene selections
        
        top1 <- topgenereact2()
        if (input$volcanoCompChoice4 == "DESeq2") {
          top_above_cutoff <- top1[which(top1$log2FoldChange > abs(FC_cutoff) & top1$pvalue < P_cutoff),]
        } else {
          top_above_cutoff <- top1[which(top1$logFC > abs(FC_cutoff) & top1$P.Value < P_cutoff),]
        }
        #top_above_cutoff <- top1[which(abs(top1$logFC) >= abs(FC_cutoff) & top1$P.Value <= P_cutoff),]
        
        Num_Genes_Above_Cutoff1 <- length(rownames(top_above_cutoff))
        
        paste("Number of Genes Above FC and P.Value Cutoffs: ",Num_Genes_Above_Cutoff1,sep = "")
        
      })
      
      output$VolGroupsText <- renderText({
        req(meta_react())
        req(metacol_reactVal())
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        covars <- input$DEGCovarSelect
        
        #print(meta)
        #print(c(colnames(meta)[1],metacol,covars))
        
        metaSub <- meta[,grepl(paste(c(colnames(meta)[1],metacol,covars),collapse = "|"),colnames(meta)), drop = F]
        metaSub_noNA <- metaSub[complete.cases(metaSub),]
        
        
        if (all(c(colnames(meta)[1],metacol) %in% colnames(metaSub_noNA))) {
          if (input$volcanoCompChoice == "Limma: Two groups") {
            groupA <- input$comparisonA2
            groupB <- input$comparisonB2
            A <- metaSub_noNA[which(metaSub_noNA[,metacol] == groupA),1]
            B <- metaSub_noNA[which(metaSub_noNA[,metacol] == groupB),1]
            form <- paste0("~0 + ",paste(c(metacol,covars),collapse = " + "))
          } else if (input$volcanoCompChoice == "Limma: One group") {
            groupA <- input$comparisonA2_one
            groupB <- paste0(unique(metaSub_noNA[which(metaSub_noNA[,metacol] != groupA),metacol]), collapse = " - ")
            #groupB <- paste0("Not_",groupA)
            metaSub_noNA[,paste0(metacol,"_Dichot")] <- NA
            metaSub_noNA[which(metaSub_noNA[,metacol] == groupA),paste0(metacol,"_Dichot")] <- groupA
            #metaSub_noNA[which(metaSub_noNA[,metacol] != groupA),paste0(metacol,"_Dichot")] <- paste0("Not_",groupA)
            metaSub_noNA[which(metaSub_noNA[,metacol] != groupA),paste0(metacol,"_Dichot")] <- groupB
            A <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == groupA),1]
            #B <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == paste0("Not_",groupA)),1]
            B <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == groupB),1]
            metacol2 <- paste0(metacol,"_Dichot")
            metaSub_noNA <- metaSub_noNA[which(metaSub_noNA[,metacol2] %in% c(groupA,groupB)),c(colnames(metaSub_noNA)[1],metacol2,covars)]
            metaSub_noNA <- metaSub_noNA[complete.cases(metaSub_noNA),]
            form <- paste0("~0 + ",paste(c(metacol,covars),collapse = " + "))
          } else if (input$volcanoCompChoice == "DESeq2") {
            metacol <- metacol_reactVal()
            metaSub <- meta[,c(colnames(meta)[1],metacol,covars),drop = F]
            metaSub_noNA <- metaSub[complete.cases(metaSub),]
            groupB <- input$DESeqDesignColRef_vol
            groupA <- paste0(unique(metaSub_noNA[which(metaSub_noNA[,metacol] != groupB),metacol]), collapse = " - ")
            #groupA <- paste0("Not_",groupB)
            metaSub_noNA[,paste0(metacol,"_Dichot")] <- NA
            metaSub_noNA[which(metaSub_noNA[,metacol] == groupB),paste0(metacol,"_Dichot")] <- groupB
            #metaSub_noNA[which(metaSub_noNA[,metacol] != groupB),paste0(metacol,"_Dichot")] <- paste0("Not_",groupB)
            metaSub_noNA[which(metaSub_noNA[,metacol] != groupB),paste0(metacol,"_Dichot")] <- groupA
            A <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == groupB),1]
            #B <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == paste0("Not_",groupB)),1]
            B <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == groupA),1]
            metacol2 <- paste0(metacol,"_Dichot")
            metaSub_noNA <- metaSub_noNA[which(metaSub_noNA[,metacol2] %in% c(groupA,groupB)),c(colnames(metaSub_noNA)[1],metacol2,covars)]
            metaSub_noNA <- metaSub_noNA[complete.cases(metaSub_noNA),]
            if (length(covars) > 0) {
              form <- paste0("~ ",paste(covars,collapse = " + "), " + ", metacol)
            } else {
              form <- paste0("~ ",metacol)
            }
          }
          
          
          
          if (nrow(metaSub_noNA) < nrow(metaSub)) {
            RowsTaken <- nrow(metaSub)-nrow(metaSub_noNA)
            paste(RowsTaken," samples removed due to NA values in covariates.",
                  "\nThis volcano plot is comparing group A: ",groupA," (N=",length(A),")", " and group B: ",groupB," (N=",length(B),")",
                  ".\nGenes with a positive log fold change are upregulated in the ",groupA,
                  " group.\nGenes with a negative log fold change are upregulated in the ",groupB, " group.",
                  "\nFormula: ",form, sep = "")
          } else {
            paste("This volcano plot is comparing group A: ",groupA," (N=",length(A),")", " and group B: ",groupB," (N=",length(B),")",
                  ".\nGenes with a positive log fold change are upregulated in the ",groupA,
                  " group.\nGenes with a negative log fold change are upregulated in the ",groupB, " group.",
                  "\nFormula: ",form, sep = "")
          }
        }
        
      })
      
      output$MAGroupsText <- renderText({
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        covars <- input$DEGCovarSelect
        metaSub <- meta[,c(colnames(meta)[1],metacol,covars)]
        metaSub_noNA <- metaSub[complete.cases(metaSub),]
        
        if (input$volcanoCompChoice == "Limma: Two groups") {
          groupA <- input$comparisonA2
          groupB <- input$comparisonB2
          A <- metaSub_noNA[which(metaSub_noNA[,metacol] == groupA),1]
          B <- metaSub_noNA[which(metaSub_noNA[,metacol] == groupB),1]
          form <- paste0("~0 + ",paste(c(metacol,covars),collapse = " + "))
        } else if (input$volcanoCompChoice == "Limma: One group") {
          groupA <- input$comparisonA2_one
          groupB <- paste0(unique(metaSub_noNA[which(metaSub_noNA[,metacol] != groupA),metacol]), collapse = " - ")
          #groupB <- paste0("Not_",groupA)
          metaSub_noNA[,paste0(metacol,"_Dichot")] <- NA
          metaSub_noNA[which(metaSub_noNA[,metacol] == groupA),paste0(metacol,"_Dichot")] <- groupA
          #metaSub_noNA[which(metaSub_noNA[,metacol] != groupA),paste0(metacol,"_Dichot")] <- paste0("Not_",groupA)
          metaSub_noNA[which(metaSub_noNA[,metacol] != groupA),paste0(metacol,"_Dichot")] <- groupB
          A <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == groupA),1]
          #B <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == paste0("Not_",groupA)),1]
          B <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == groupB),1]
          metacol2 <- paste0(metacol,"_Dichot")
          metaSub_noNA <- metaSub_noNA[which(metaSub_noNA[,metacol2] %in% c(groupA,groupB)),c(colnames(metaSub_noNA)[1],metacol2,covars)]
          metaSub_noNA <- metaSub_noNA[complete.cases(metaSub_noNA),]
          form <- paste0("~0 + ",paste(c(metacol,covars),collapse = " + "))
        } else if (input$volcanoCompChoice == "DESeq2") {
          metacol <- metacol_reactVal()
          metaSub <- meta[,c(colnames(meta)[1],metacol,covars), drop = F]
          metaSub_noNA <- metaSub[complete.cases(metaSub),]
          groupB <- input$DESeqDesignColRef_vol
          groupA <- paste0(unique(metaSub_noNA[which(metaSub_noNA[,metacol] != groupB),metacol]), collapse = " - ")
          #groupA <- paste0("Not_",groupB)
          metaSub_noNA[,paste0(metacol,"_Dichot")] <- NA
          metaSub_noNA[which(metaSub_noNA[,metacol] == groupB),paste0(metacol,"_Dichot")] <- groupB
          #metaSub_noNA[which(metaSub_noNA[,metacol] != groupB),paste0(metacol,"_Dichot")] <- paste0("Not_",groupB)
          metaSub_noNA[which(metaSub_noNA[,metacol] != groupB),paste0(metacol,"_Dichot")] <- groupA
          A <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == groupB),1]
          #B <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == paste0("Not_",groupB)),1]
          B <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == groupA),1]
          metacol2 <- paste0(metacol,"_Dichot")
          metaSub_noNA <- metaSub_noNA[which(metaSub_noNA[,metacol2] %in% c(groupA,groupB)),c(colnames(metaSub_noNA)[1],metacol2,covars)]
          metaSub_noNA <- metaSub_noNA[complete.cases(metaSub_noNA),]
          if (length(covars) > 0) {
            form <- paste0("~ ",paste(covars,collapse = " + "), " + ",metacol)
          } else {
            form <- paste0("~ ",metacol)
          }
        }
        
        if (nrow(metaSub_noNA) < nrow(metaSub)) {
          RowsTaken <- nrow(metaSub)-nrow(metaSub_noNA)
          paste(RowsTaken," samples removed due to NA values in covariates.",
                "\nThis MA plot is comparing group A: ",groupA," (N=",length(A),")", " and group B: ",groupB," (N=",length(B),")",
                ".\nGenes with a positive log fold change are upregulated in the ",groupA,
                " group.\nGenes with a negative log fold change are upregulated in the ",groupB, " group.",
                "\nFormula: ",form, sep = "")
        } else {
          paste("This MA plot is comparing group A: ",groupA," (N=",length(A),")", " and group B: ",groupB," (N=",length(B),")",
                ".\nGenes with a positive log fold change are upregulated in the ",groupA,
                " group.\nGenes with a negative log fold change are upregulated in the ",groupB, " group.",
                "\nFormula: ",form, sep = "")
        }
        
      })
      
      output$upregpath_text <- renderText({
        paste("Genes in these enriched terms are upregulated in group A: ", input$comparisonA2.path, " group.", sep = "")
      })
      
      output$downregpath_text <- renderText({
        paste("Genes in these enriched terms are upregulated in group B: ", input$comparisonB2.path," group", sep = "")
      })
      
      output$degtext <- renderText({
        
        meta <- meta_react()
        metacol <- metacol_reactVal()
        if (length(colnames(meta)) == 2) {
          metacol <- colnames(meta)[2]
        }
        covars <- input$DEGCovarSelectTab
        
        metaSub <- meta[,c(colnames(meta)[1],metacol,covars)]
        metaSub_noNA <- metaSub[complete.cases(metaSub),]
        
        if (input$volcanoCompChoice2 == "Limma: Two groups") {
          groupA <- input$comparisonA2.DEG
          groupB <- input$comparisonB2.DEG
          A <- metaSub_noNA[which(metaSub_noNA[,metacol] == groupA),1]
          B <- metaSub_noNA[which(metaSub_noNA[,metacol] == groupB),1]
          form <- paste0("~0 + ",paste(c(metacol,covars),collapse = " + "))
        } else if (input$volcanoCompChoice2 == "Limma: One group") {
          groupA <- input$comparisonA2.DEG_one
          groupB <- paste0(unique(metaSub_noNA[which(metaSub_noNA[,metacol] != groupA),metacol]), collapse = " - ")
          #groupB <- paste0("Not_",groupA)
          metaSub_noNA[,paste0(metacol,"_Dichot")] <- NA
          metaSub_noNA[which(metaSub_noNA[,metacol] == groupA),paste0(metacol,"_Dichot")] <- groupA
          #metaSub_noNA[which(metaSub_noNA[,metacol] != groupA),paste0(metacol,"_Dichot")] <- paste0("Not_",groupA)
          metaSub_noNA[which(metaSub_noNA[,metacol] != groupA),paste0(metacol,"_Dichot")] <- groupB
          A <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == groupA),1]
          #B <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == paste0("Not_",groupA)),1]
          B <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == groupB),1]
          metacol2 <- paste0(metacol,"_Dichot")
          metaSub_noNA <- metaSub_noNA[which(metaSub_noNA[,metacol2] %in% c(groupA,groupB)),c(colnames(metaSub_noNA)[1],metacol2,covars)]
          metaSub_noNA <- metaSub_noNA[complete.cases(metaSub_noNA),]
          form <- paste0("~0 + ",paste(c(metacol,covars),collapse = " + "))
        } else if (input$volcanoCompChoice2 == "DESeq2") {
          metacol <- metacol_reactVal()
          metaSub <- meta[,c(colnames(meta)[1],metacol,covars), drop = F]
          metaSub_noNA <- metaSub[complete.cases(metaSub),]
          groupB <- input$DESeqDesignColRef_tab
          groupA <- paste0(unique(metaSub_noNA[which(metaSub_noNA[,metacol] != groupB),metacol]), collapse = " - ")
          #groupA <- paste0("Not_",groupB)
          metaSub_noNA[,paste0(metacol,"_Dichot")] <- NA
          metaSub_noNA[which(metaSub_noNA[,metacol] == groupB),paste0(metacol,"_Dichot")] <- groupB
          #metaSub_noNA[which(metaSub_noNA[,metacol] != groupB),paste0(metacol,"_Dichot")] <- paste0("Not_",groupB)
          metaSub_noNA[which(metaSub_noNA[,metacol] != groupB),paste0(metacol,"_Dichot")] <- groupA
          A <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == groupB),1]
          #B <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == paste0("Not_",groupB)),1]
          B <- metaSub_noNA[which(metaSub_noNA[,paste0(metacol,"_Dichot")] == groupA),1]
          metacol2 <- paste0(metacol,"_Dichot")
          metaSub_noNA <- metaSub_noNA[which(metaSub_noNA[,metacol2] %in% c(groupA,groupB)),c(colnames(metaSub_noNA)[1],metacol2,covars)]
          metaSub_noNA <- metaSub_noNA[complete.cases(metaSub_noNA),]
          if (length(covars) > 0) {
            form <- paste0("~ ",paste(covars,collapse = " + ")," + ",metacol)
          } else {
            form <- paste0("~ ",metacol)
          }
        }
        
        
        
        if (nrow(metaSub_noNA) < nrow(metaSub)) {
          RowsTaken <- nrow(metaSub)-nrow(metaSub_noNA)
          paste(RowsTaken," samples removed due to NA values in covariates.",
                "This table represents differentially expressed genes when comparing group A: ",groupA," (N=",length(A),")"," and group B: ",input$comparisonB2.DEG," (N=",length(B),")",
                ".\nGenes with a positive logFC, indicate an upregulation in group A: ", groupA,
                ".\nGenes with a negative logFC, indicate an upregulation in group B: ", groupB,
                ".\nThe 'AveExpr' column represents the log transformed average expression between group A: ",groupA," and group B: ",input$comparisonB2.DEG,".",
                "\nFormula: ",form, sep = "")
        } else {
          paste("This table represents differentially expressed genes when comparing group A: ",groupA," (N=",length(A),")"," and group B: ",groupB," (N=",length(B),")",
                ".\nGenes with a positive logFC, indicate an upregulation in group A: ", groupA,
                ".\nGenes with a negative logFC, indicate an upregulation in group B: ", groupB,
                ".\nThe 'AveExpr' column represents the log transformed average expression between group A: ",groupA," and group B: ",groupB,".",
                "\nFormula: ",form, sep = "")
        }
      })
      
      output$UpRegPathLabel <- renderUI({
        
        FC <- input$pathFC
        h3(paste("(Up-regulated pathway (> ",FC," logFC)",sep = ""))
        
      })
      
      output$DnRegPathLabel <- renderUI({
        
        FC <- input$pathFC
        h3(paste("Down-regulated pathway (> -",FC," logFC)",sep = ""))
        
      })
      
    }
  })
  
  
}

shinyApp(ui = ui, server = server)


