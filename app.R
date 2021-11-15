
library(shiny)
library(shinythemes)
library(shinyjqui)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(dplyr)
library(DT)
library(enrichplot)
library(shinycssloaders)
library(ggplot2)
library(ggpubr)
library(msigdbr)
library(reshape2)
library(tibble)
library(plotly)
library(readr)
library(limma)
library(enrichR)
library(ggrepel)
library(tidyr)
library(GSVA)
library(tools)



####----User Data Input----####

#Input desired project name for webpage - will be followed by 'Expression Analysis'
ProjectName <- "USP7 Human Demo"


##--User Input File Names--##

#expression data
expr_file <- "~/R/GeneralExpressionAnalysisApp-main/ExampleData/htseq_gene_level_fpkm_T_geneName_max_1cutoff_v2.txt"

#meta data
meta_file <- "~/R/GeneralExpressionAnalysisApp-main/ExampleData/USP7_meta.tsv"
#Is there a header?
header <- TRUE

#Enriched Signatures data table
ES_tables <- c("~/R/GeneralExpressionAnalysisApp-main/ExampleData/USP7_Enrich_Sig.tsv")

#If human: set TRUE
#If mouse: set FALSE
human <- TRUE

##--User Gene Set Input--##

#write in the name of your gene set list for shiny UI
userGSlist_name <- 'CellMarker Gene Sets'

#path to your gene set file .gmt or .txt/.tsv
userGS_file <- '~/R/GeneralExpressionAnalysisApp-main/GeneSets/CellMarker_gsNsym_HS.tsv'
#Does gene set file have header?
header.gs <- TRUE

#path to your R data list object for ssGSEA
userRData_file <- '~/R/GeneralExpressionAnalysisApp-main/GeneSets/CellMarker_GS_HS.RData'





####----Backend Data Input----####


if (human == TRUE) {
    #MSigDB gene set
    msigdb <- '~/R/GeneralExpressionAnalysisApp-main/GeneSets/msigdb_gsNsym_HS.zip'
    #MSigDB gene set FOR UI
    msigdb2 <- '~/R/GeneralExpressionAnalysisApp-main/GeneSets/msigdb_gsNcat_HS.tsv'
    #gene set list for ssGSEA
    load('~/R/GeneralExpressionAnalysisApp-main/GeneSets/gs_list_HS.RData')
    #Cytokine genes for human
    CTKgenes <- c("IL2","IL12A","IL12B","IL17A","IFNA1","IFNB1","IFNG","IFNGR","CD11b",
                  "ITGAM","CD33","ENTPD1","ICOSLG","CD275","CD278","TNFSF9","TNFRSF9",
                  "CD40","CD40LG","CD70","CD27","TNFSF18","TNFRSF18","TNFSF14","TNFRSF14",
                  "TNFSF4","TNFRSF4","HLA-A","CD3","CEACAM1","CD80","CD86","CTLA4","CD276",
                  "VTCN1","PVR","CD226","TIGIT","CD96","LGALS3","LGALS3BP","LGALS9","LGALS9C",
                  "HAVCR2","HHLA2","TMIGD2","CD274","PDCD1LG2","PDCD1","VSIR")
}
if (human == FALSE) {
    #MSigDB gene set
    msigdb <- '~/R/GeneralExpressionAnalysisApp-main/GeneSets/msigdb_gsNsym_MM.zip'
    #MSigDB gene set FOR UI
    msigdb2 <- '~/R/GeneralExpressionAnalysisApp-main/GeneSets/msigdb_gsNcat_MM.tsv'
    #gene set list for ssGSEA 
    load('~/R/GeneralExpressionAnalysisApp-main/GeneSets/gs_list_MM.RData')
    #Cytokiny genes for mouse
    CTKgenes <- c("Il2","Il12a","Il12b","Il17a","Ifna13","Ifnb1","Ifng","Ifngr1","Cd11b","Itgam",
                  "Cd33","Entpd1","Icosl","Icos","Tnfsf9","Tnfrsf9","Cd40","Cd40lg","Cd70","Cd27",
                  "Tnfsf18","Tnfrsf18","Tnfsf14","Tnfrsf14","Tnfsf4","Tnfrsf4","H2-K1","CD3G",
                  "Ceacam1","Cd80","Cd86","Ctla4","Cd276","Vtcn1","Pvr","Cd226","Tigit","Cd96","Lgals3",
                  "Lgals3bp","Lgals9","Lgals9c","Havcr2","Hhla2","Cd274","Pdcd1lg2","Pdcd1","Vsir")
}




####----Read and Manipulate Files----####


##read files

#expression data
expr <- read.delim(expr_file, sep = '\t', header = T, strip.white = T)
#expr <- as.data.frame(read_delim(expr_file, delim = '\t', trim_ws = T))
colnames(expr)[1] <- "Gene"
expr <- expr %>%
    drop_na(Gene)
row.names(expr) <- expr[,1]
expr <- expr[,-1]
colnames(expr) <- gsub("[_.-]", "_", colnames(expr))
A <- as.matrix(expr)



#meta
meta <- read.delim(meta_file, sep = '\t', header = header, strip.white = T)
meta[,1] <- gsub("[_.-]", "_", meta[,1])
colnames(meta) <- c("SampleName","Group")
metagroups <- as.vector(levels(factor(meta[,2])))


#boxplot choices based on meta groups
if (length(metagroups) == 2) {
    boxopt <- c("wilcox.test", "t.test", "none")
}
if (length(metagroups) >= 3) {
    boxopt <- c("kruskal.test", "anova", "none")
}



#for heatmap sample selection
sampsames <- intersect(colnames(expr),meta[,1])


#Enriched Signatures
if (!is.null(ES_tables)) {
    ldf <- list()
    SigNames <- c()
    for (k in 1:length(ES_tables)){
        file <- basename(ES_tables[k])
        file2 <- gsub("\\..*$","",file)
        SigNames <- c(SigNames, file2)
        #ldf[[k]] <- read.delim(ES_tables[k], header = T, sep = '\t')
        ldf[[k]] <- as.data.frame(read_delim(ES_tables[k], delim = '\t'))
    }
    j <- 1
    for (i in 1:length(ldf)){
        names(ldf)[i] <- paste("ES_table",j,sep = "")
        j=j+1
    }
    list2env(ldf,globalenv())
}
#load blank enriched signatures table if none given by user
if (is.null(ES_tables)) {
    SigNames <- "None Loaded"
    ES_table1 <- data.frame("ID" = "No Enriched Signature Tables loaded","Description" = "NA","setSize" = "NA",
                            "enrichmentScore" = "NA","NES" = "NA","pvalue" = "NA","p.adjust" = "NA","qvalues" = "NA",
                            "rank" = "NA","leading_edge" = "NA","core_enrichment" = "NA")
    
}


#MSigDB gene sets
msigdb.gsea <- read_delim(msigdb, delim = '\t')
gmt <- msigdb.gsea


#MSigDB gene sets FOR UI
msigdb.gsea2 <- read_delim(msigdb2, delim = '\t')


#tab2 User gene set
if (file_ext(userGS_file) == "gmt") {
    tab2 <- read.gmt(userGS_file)
}
if (file_ext(userGS_file) == "tsv" || file_ext(userGS_file) == "txt") {
    tab2 <- read.delim(userGS_file, header = header.gs, sep = '\t')
}

#tab2 back end
GeneSet2 <- as.data.frame(unique(tab2[,1]))
rownames(GeneSet2) <- 1:nrow(GeneSet2)
colnames(GeneSet2)[1] <- "Gene_Set"

#tab2 R Data list
loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}
gs2 <- loadRData(userRData_file)


##gene list file from expr data
Gene <- rownames(expr)
geneList <- as.data.frame(Gene)


#CV function for variance
cv <- function(x){
    (sd(x)/mean(x))*100
}




####----Shiny UI----####

shinytheme("sandstone")

ui <-
    
navbarPage(paste("{",ProjectName,"Expression Analysis }", sep=" "),
    
    ####----Intro Tab----####
    
    tabPanel("Intro/Methodology",
        fluidPage(
            mainPanel(
                h2("GSEA and RNAseq Analysis Method"),
                p("Mapped reads were derived from BBSR and likely normalized to CPM. Pathway enrichment analysis was performed by GSEA [PMID: 16199517] and enrichR[PMID: 23586463]. Single-sample gene set enrichment analysis (ssGSEA) [PMID: 19847166, PMID: 30595505] was used to quantify the expression signatures. LIMMA was used to define differentially expressed genes [PMID: PMID: 25605792]")
                )
            )
        ),
    
    ####----GSEA Tab----####
    
    tabPanel("GSEA Analysis",
        fluidPage(
            title = "GSEA Analysis",
            sidebarLayout(
                sidebarPanel(
                    selectInput("comparisonA", "Comparison: GroupA",
                                choices = metagroups, selected = metagroups[1]),
                    selectInput("comparisonB", "Comparison: GroupB",
                                choices = metagroups, selected = metagroups[2]),
                    numericInput("userPval", "Pvalue Cutoff", value = 1.0, width = '100%'),
                    selectInput("boxplotcompare", "Boxplot Stat Compare Method:",
                                choices = boxopt),
                    tabsetPanel(
                        id = "tables",
                        tabPanel("MSigDB Gene Sets",
                                 div(DT::dataTableOutput("msigdbTable"), style = "font-size:10px; height:500px; overflow-X: scroll"),
                                 value = 1),
                        tabPanel(paste(userGSlist_name),
                                 div(DT::dataTableOutput("tab2table"), style = "font-size:10px; height:500px; overflow-X: scroll"),
                                 value = 3),
                        tabPanel("Use your own gene set",
                                 p(),
                                 uiOutput("user.gmt"),
                                 uiOutput("user.GStable"),
                                 uiOutput("user.gs.tsv"),
                                 uiOutput("user.gs.tsv.table"),
                                 value = 5)
                    )
                ),
                mainPanel(
                    tabsetPanel(
                        id = "dataset1",
                        tabPanel("Enrichment Plot",
                                 h3("GSEA Enrichment Plot"),
                                 verbatimTextOutput("NESandPval"),
                                 withSpinner(plotOutput("enrichplot0", width = "500px", height = "450px"), type = 6),
                                 h3("Leading Edge Genes (~Signal2Noise Ranking)"),
                                 downloadButton("LEGdownload", "Download .tsv"),
                                 div(DT::dataTableOutput("LeadingEdgeGenes"), style = "font-size:12px; height:500px; overflow-y: scroll"),
                                 value = 2),
                        tabPanel("Gene Expression Heatmap",
                                 withSpinner(jqui_resizable(plotOutput("heatmap0", width = "100%", height = "1000px")), type = 6),
                                 value = 4),
                        tabPanel("GSEA Summary Table",
                                 p(),
                                 selectInput("SigTableChoice", "Select Enriched Signatures Table:",
                                             choices = SigNames),
                                 h3("Enriched Signatures Table"),
                                 DT::dataTableOutput("enrich_sig_table", width = "100%"),
                                 downloadButton("enrich_sig_download","Download .tsv"),
                                 value = 6),
                        tabPanel("Generate Summary Table",
                                 p("Please note this may take several minutes depending on size and quantity of gene sets in GMT file."),
                                 withSpinner(DT::dataTableOutput("enrich_sig_table_gen"), type = 6),
                                 downloadButton("enrich_sig_download.u","Download .tsv"),
                                 value = 8),
                        tabPanel("ssGSEA Boxplots",
                                 p(),
                                 withSpinner(jqui_resizable(plotOutput('boxplot2', width = "700px", height = "500px")), type = 6),
                                 DT::dataTableOutput("ssGSEAtable"),
                                 downloadButton("ssGSEAdownload", "Download .tsv"))
                    )
                )
            )
        )
     ),
    
    ####----RNAseq Tab----####
    
    tabPanel("RNAseq Analysis",
        fluidPage(
            title = "RNAseq Analysis",
            sidebarLayout(
                sidebarPanel(
                    width = 3,
                    tabsetPanel(
                        id = "customs",
                        tabPanel("RNAseq Parameters",
                                 p(),
                                 selectInput("comparisonA2", "Comparison: GroupA",
                                             choices = metagroups, selected = metagroups[1]),
                                 selectInput("comparisonB2", "Comparison: GroupB",
                                             choices = metagroups, selected = metagroups[2]),
                                 numericInput("NumFeatures", step = 1, label = "Number of Genes", value = 100),
                                 selectInput("VarianceMeasure", "Selecte Variance Measure",
                                             choices = c("MAD","CV","VAR")),
                                 selectInput("ClusteringMethod",
                                             "Select clustering Method",
                                             choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")),
                                 numericInput("NumClusters", step = 1, label = "Number of Clusters (Cut Tree with ~k)", value = 2),
                                 h3("Cluster Result:"),
                                 div(DT::dataTableOutput("Clusters"), style = "font-size:10px"),
                                 h3("Most Variable Genes:"),
                                 div(DT::dataTableOutput("MostVariableGenesList"), style = "font-size:10px; height:500px; overflow-y: scroll"),
                                 p(),
                                 downloadButton("MVGdownload", "Download MVG .tsv"),
                                 downloadButton("MVGdownloadgmt", "Download MVG .gmt"),
                                 h3("DEG .gmt Download Parameters"),
                                 textInput("DEGfileName", "Input File Name:",value = "DEGgeneSet"),
                                 numericInput("fc_cutoff2", "LogFC Threshold (Absolute Value)",
                                              min = 0, max = 5, step = 0.1, value = 1),
                                 selectInput("UpDnChoice","Up-regulated or Down-regulated:",
                                             choices = c("UpAndDown_Regulated","Up_Regulated","Down_Regulated")),
                                 numericInput("p_cutoff2", "Adj.P.Value Cutoff:",
                                              min = 0, max = 10, step = 0.1, value = 0.05),
                                 numericInput("top_x2", "Number of Top Hits:", value = 100),
                                 downloadButton("DEGgmtDownload", "Download DEG .gmt")
                        ),
                        tabPanel("Custom Heatmap Genes",
                                 p(),
                                 selectizeInput("heatmapGeneSelec","Gene Selection:",
                                                choices = sort(as.vector(geneList[,1])),
                                                multiple = T, selected = CTKgenes),
                                 textInput("userheatgenes", "Text Input of Gene List (space or tab delimited):", value = ""),
                                 #textInput("userheatsamp", "Samples Selection (space delimited):", value = sampsames),
                                 selectizeInput("userheatsamp2", "Samples Selection:",
                                                choices = sampsames, multiple = T, selected = sampsames),
                                 value = 444
                        ),
                        tabPanel("Volcano & MA Plot Parameters",
                                 p(),
                                 numericInput("fc_cutoff", "LogFC Threshold (Absolute Value)",
                                              min = 0, max = 5, step = 0.1, value = 1),
                                 numericInput("p_cutoff", "Significance Threshold (-log10(P.Value)):",
                                              min = 0, max = 10, step = 0.1, value = 0.05),
                                 numericInput("top_x", "Number of Top Hits:", value = 10),
                                 selectizeInput("userGeneSelec", "User Selected Hits:",
                                                choices = sort(as.vector(geneList[,1])), multiple = T, selected = "-")
                        )
                    )
                ),
                mainPanel(
                    tabsetPanel(
                        id = "dataset2",
                        tabPanel("Heatmap",
                                 withSpinner(jqui_resizable(plotOutput("heatmap1", width = "100%", height = "1000px")), type = 6),
                                 value = 22),
                        tabPanel("Volcano Plot",
                                 verbatimTextOutput("VolGroupsText"),
                                 jqui_resizable(plotOutput('Volcano3', width = "800px", height = "550px",
                                                           hover = hoverOpts("plot_hover", delay = 10, delayType = "debounce"))),
                                 uiOutput("hover_info"),
                                 value = 44),
                        tabPanel("MA Plot",
                                 verbatimTextOutput("MAGroupsText"),
                                 jqui_resizable(plotOutput('MAPlot1', width = "800px", height = "550px",
                                                           hover = hoverOpts("plot_hover2", delay = 10, delayType = "debounce"))),
                                 uiOutput("hover_info2"),
                                 value = 66),
                        tabPanel("Box Plots",
                                 div(DT::dataTableOutput("GeneListTable"), style = "font-size:10px"),
                                 withSpinner(jqui_resizable(plotOutput('boxplot1', width = "700px", height = "400px")), type = 6),
                                 value = 88),
                        tabPanel("Gene Scatter Plot",
                                 fluidRow(
                                     column(4,
                                            p(),
                                            selectInput("scatterG1","Select Gene 1",
                                                        choices = Gene)),
                                     column(4,
                                            p(),
                                            selectInput("scatterG2","Select Gene 2",
                                                        choices = Gene, selected = Gene[2])),
                                     column(4,
                                            p(),
                                            checkboxInput("logask","Log2 Transform Expression Data"))
                                 ),
                                 withSpinner(jqui_resizable(plotlyOutput("geneScatter0", height = "600px")), type = 6),
                                 DT::dataTableOutput("geneScatterTable"),
                                 downloadButton("geneScatterDownload", "Download Non-log2 Transformed .tsv")),
                        tabPanel("Pathway Analysis",
                                 selectInput("SelectedPathway", "Select Pathway",
                                             choices = c("ChEA_2016", "MSigDB_Hallmark_2020","KEGG_2019_Human",
                                                         "KEGG_2019_Mouse","GO_Biological_Process_2021","GO_Cellular_Component_2021",
                                                         "GO_Molecular_Function_2021","TRANSFAC_and_JASPAR_PWMs", "TRRUST_Transcription_Factors_2019",
                                                         "DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019", "CellMarker_Augmented_2021")),
                                 h3("Up-regulated pathway (> 1 logFC)"),
                                 verbatimTextOutput("upregpath_text"),
                                 withSpinner(plotOutput('UpRegPathway1'), type = 6),
                                 div(DT::dataTableOutput("UpRegPathwayTable1"), style = "font-size:10px"),
                                 downloadButton("UpRegPathDownload", "Download Table .tsv"),
                                 downloadButton("UpRegPathDownloadgmt", "Download .gmt"),
                                 h3("Down-regulated pathway (< -1 logFC)"),
                                 verbatimTextOutput("downregpath_text"),
                                 withSpinner(plotOutput('DnRegPathway1'), type = 6),
                                 div(DT::dataTableOutput("DnRegPathwayTable1"), style = "font-size:10px"),
                                 downloadButton("DnRegPathDownload", "Download Table .tsv"),
                                 downloadButton("DnRegPathDownloadgmt", "Download .gmt"),
                                 value = 00),
                        tabPanel("DEG Table",
                                 p(),
                                 verbatimTextOutput("degtext"),
                                 div(DT::dataTableOutput("DEGtable1"), style = "font-size:12px"),
                                 downloadButton("DEGtableDownload", "Download DEG table .tsv"),
                                 value = 24)
                    )
                )
            )
        )
    )
)



####----Server----####


server <- function(input, output, session) {
    
    
    ####----Render UI----####
    
    
    #render user gmt data upload if indicated
    output$user.gmt <- renderUI({
        fileInput("user.gmt.file", "Gene Set File (.gmt, .tsv, or .txt)", accept = c(".gmt",".tsv",".txt"))
    })
    
    #render gene set table based off gmt file given
    output$user.GStable <- renderUI({
        req(input$user.gmt.file)
        div(DT::dataTableOutput("GStable.u"), style = "font-size:10px; height:500px; overflow-X: scroll")
    })
    
    #render UI for hover text in volcano plot
    output$hover_info <- renderUI({
        top2 <- topgenereact()
        df <- top2 %>%
            select(GeneName,logFC,P.Value,adj.P.Val)
        colnames(df)[3] <- "-log10(P.Value)"
        df$P.Value <- df$'-log10(P.Value)'
        df$`-log10(P.Value)` <- -log10(df$`-log10(P.Value)`)
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
        top2 <- topgenereact()
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
    
    
    ####----Reactives----####
    
    
    #reactive for ssGSEA function
    ssGSEAfunc <- reactive({
        if (input$tables == 1) {
            GS <- gs[(msigdb.gsea2[input$msigdbTable_rows_selected,3])]
        }
        if (input$tables == 3) {
            GS <- gs2[(GeneSet2[input$tab2table_rows_selected,1])]
        }
        gsva(A, GS, method = "ssgsea", verbose = F)
    })
    
    #create background GMT from user input gene set table
    GStable.ubg <- reactive({
        gmt.u <- input$user.gmt.file
        ext <- tools::file_ext(gmt.u$datapath)
        req(gmt.u)
        validate(need(ext == c("gmt","tsv","txt"), "Please upload .gmt, .tsv, or .txt file"))
        if (ext == "gmt") {
            read.gmt(gmt.u$datapath)
        }
        else {
            read_delim(gmt.u$datapath, delim = '\t')
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
        groupA <- meta[,1][meta[,2] == input$comparisonA]
        groupB <- meta[,1][meta[,2] == input$comparisonB]
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
        s2n.df <- as.data.frame(s2n.matrix)
        s2n.df$GeneID <- rownames(s2n.df)
        rownames(s2n.df) <- NULL
        data <- dplyr::select(s2n.df, GeneID, V1)
        data.gsea <- data$V1
        names(data.gsea) <- as.character(data$GeneID)
        s2n.matrix.s <- sort(data.gsea, decreasing = T)
        ##----GSEA----##
        if (input$tables == 1){
            GSEA(s2n.matrix.s, TERM2GENE = gmt[which(gmt$gs_name == msigdb.gsea2[input$msigdbTable_rows_selected,3]),],
                 verbose = F, pvalueCutoff = input$userPval)
        }
        else if (input$tables == 3){
            GSEA(s2n.matrix.s, TERM2GENE = tab2[which(tab2$term == as.character(GeneSet2[input$tab2table_rows_selected,1])),],
                 verbose = F, pvalueCutoff = input$userPval)
        }
        else if (input$tables == 5){
            GSEA(s2n.matrix.s, TERM2GENE = GStable.ubg()[which(GStable.ubg()$term == as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])),],
                 verbose = F, pvalueCutoff = input$userPval)
        }
    })
    
    #top genes data frame reactive
    topgenereact <- reactive({
        top_probes <- input$NumFeatures
        col_labels <- colnames(expr)
        isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
        exp <- expr[isexpr,]
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
        colnames(zdataset) <- names(dataset)
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
        bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
        hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
        results2 <- hclust(dist(t(dataset)), method = input$ClusteringMethod)
        m <- sort(cutree(results2, k = input$NumClusters))
        output <- cbind(colnames(m), as.matrix(m))
        colnames(output) <- c("Cluster")
        
        A <- meta[,1][meta[,2] == input$comparisonA2]
        B <- meta[,1][meta[,2] == input$comparisonB2]
        
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
        top2 <- top1
        top2["GeneName"] <- rownames(top2)
        Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
        top2['group'] <- "NotSignificant"
        top2[which(top2$P.Value < input$p_cutoff & abs(top2$logFC) < abs(input$fc_cutoff)), "group"] <- "Significant"
        top2[which(top2$P.Value > input$p_cutoff & abs(top2$logFC) > abs(input$fc_cutoff)), "group"] <- "FoldChange"
        top2[which(top2$P.Value < input$p_cutoff & abs(top2$logFC) > abs(input$fc_cutoff)), "group"] <- "Significant&FoldChange"
        top2['FCgroup'] <- "NotSignificant"
        top2[which(abs(top2$logFC) > abs(input$fc_cutoff)), "group2"] <- "FoldChange"
        top2
    })
    
    
    ####----Data Tables----####
    
    
    #render MSigDB gene set table
    output$msigdbTable <- DT::renderDataTable({
        DT::datatable(msigdb.gsea2,
                      selection = 'single',
                      options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100")))
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
            gmt.us <- read_delim(gmt.u$datapath, delim = '\t')
        }
        GeneSet <- as.data.frame(unique(gmt.us[1]))
        rownames(GeneSet) <- 1:nrow(GeneSet)
        colnames(GeneSet)[1] <- "Gene_Set"
        DT::datatable(GeneSet,
                      selection = 'single',
                      options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100")))
    })
    
    #render tab2 gene set table
    output$tab2table <- DT::renderDataTable({
        DT::datatable(GeneSet2,
                      selection = 'single',
                      options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100")))
        
    })
    
    #render leading edge genes list
    output$LeadingEdgeGenes <- DT::renderDataTable({
        if (input$tables == 1){
            if (length(input$msigdbTable_rows_selected) > 0){
                res <- datasetInput()
                gsea.df <- as.data.frame(res@result)
                GS <- msigdb.gsea2[input$msigdbTable_rows_selected,3]
                ## Subset core enriched genes
                genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
                genes2 <- strsplit(genes1,"/")
                GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
                GeneSymbol$Rank <- rownames(GeneSymbol)
                GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
                DT::datatable(GeneSymbol, options = list(paging = F), rownames = F)
            }
        }
        else if (input$tables == 3){
            if (length(input$tab2table_rows_selected) > 0){
                res <- datasetInput()
                gsea.df <- as.data.frame(res@result)
                GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
                ## Subset core enriched genes
                genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
                genes2 <- strsplit(genes1,"/")
                GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
                GeneSymbol$Rank <- rownames(GeneSymbol)
                GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
                DT::datatable(GeneSymbol, options = list(paging = F), rownames = F)
            }
        }
        else if (input$tables == 5){
            if (length(input$GStable.u_rows_selected) > 0){
                res <- datasetInput()
                gsea.df <- as.data.frame(res@result)
                GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
                ## Subset core enriched genes
                genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
                genes2 <- strsplit(genes1,"/")
                GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
                GeneSymbol$Rank <- rownames(GeneSymbol)
                GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
                DT::datatable(GeneSymbol, options = list(paging = F), rownames = F)
            }
        }
    })
    
    #render pre-loaded enriched signatures table
    output$enrich_sig_table <- DT::renderDataTable({
        if (input$SigTableChoice == SigNames[1]) {
            gsea.df <- as_tibble(ES_table1)
            DT::datatable(gsea.df,
                          extensions = c("KeyTable", "FixedHeader"),
                          caption = "Enriched Signatures",
                          options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"),scrollX = T)) %>%
                formatRound(columns = c(2:10), digits = 2)
        }
        else if (input$SigTableChoice == SigNames[2]) {
            gsea.df <- as_tibble(ES_table2)
            DT::datatable(gsea.df,
                          extensions = c("KeyTable", "FixedHeader"),
                          caption = "Enriched Signatures",
                          options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"),scrollX = T)) %>%
                formatRound(columns = c(2:10), digits = 2)
        }
        else if (input$SigTableChoice == SigNames[3]) {
            gsea.df <- as_tibble(ES_table3)
            DT::datatable(gsea.df,
                          extensions = c("KeyTable", "FixedHeader"),
                          caption = "Enriched Signatures",
                          options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"),scrollX = T)) %>%
                formatRound(columns = c(2:10), digits = 2)
        }
        
    })
    
    #render user generated enriched signature table based off other gmt
    output$enrich_sig_table_gen <- DT::renderDataTable({
        if (input$tables == 3) {
            groupA <- meta[,1][meta[,2] == input$comparisonA]
            groupB <- meta[,1][meta[,2] == input$comparisonB]
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
            s2n.df <- as.data.frame(s2n.matrix)
            s2n.df$GeneID <- rownames(s2n.df)
            rownames(s2n.df) <- NULL
            data <- dplyr::select(s2n.df, GeneID, V1)
            data.gsea <- data$V1
            names(data.gsea) <- as.character(data$GeneID)
            s2n.matrix.s <- sort(data.gsea, decreasing = T)
            ##----GSEA----##
            gmt.i <- tab2
            gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
            gsea.df <- as_tibble(gsea.res@result)
            ## displaying the GSEA results as interactive data table
            DT::datatable(gsea.df,
                          extensions = c("KeyTable", "FixedHeader"),
                          caption = "Enriched Signatures",
                          options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T)) %>%
                formatRound(columns = c(2:10), digits = 2)
        }
        else if (input$tables == 5) {
            groupA <- meta[,1][meta[,2] == input$comparisonA]
            groupB <- meta[,1][meta[,2] == input$comparisonB]
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
            s2n.df <- as.data.frame(s2n.matrix)
            s2n.df$GeneID <- rownames(s2n.df)
            rownames(s2n.df) <- NULL
            data <- dplyr::select(s2n.df, GeneID, V1)
            data.gsea <- data$V1
            names(data.gsea) <- as.character(data$GeneID)
            s2n.matrix.s <- sort(data.gsea, decreasing = T)
            ##----GSEA----##
            gmt.i <- GStable.ubg()
            gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
            gsea.df <- as_tibble(gsea.res@result)
            ## displaying the GSEA results as interactive data table
            DT::datatable(gsea.df,
                          extensions = c("KeyTable", "FixedHeader"),
                          caption = "Enriched Signatures",
                          options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T)) %>%
                formatRound(columns = c(2:10), digits = 2)
        }
    })
    
    #render variable genes list in sidebar from heatmap
    output$MostVariableGenesList <- DT::renderDataTable({
        top_probes <- input$NumFeatures
        col_labels <- colnames(expr)
        isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
        exp <- expr[isexpr,]
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
        DT::datatable(variable_gene_list, options = list(paging = F), rownames = F)
    })
    
    #render DEG table
    output$DEGtable1 <- DT::renderDataTable({
        top_probes <- input$NumFeatures
        col_labels <- colnames(expr)
        isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
        exp <- expr[isexpr,]
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
        colnames(zdataset) <- names(dataset)
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
        A <- meta[,1][meta[,2] == input$comparisonA2]
        B <- meta[,1][meta[,2] == input$comparisonB2]
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
        DT::datatable(top1, options = list(lengthMenu = c(50,100,1000, 5000, 10000), pageLength = 100, scrollX = TRUE),
                      selection=list(mode = "multiple"))
    })
    
    #render up regulated pathway enrichment data table
    output$UpRegPathwayTable1 <- DT::renderDataTable({
        top_probes <- input$NumFeatures
        col_labels <- colnames(expr)
        isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
        exp <- expr[isexpr,]
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
        colnames(zdataset) <- names(dataset)
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
        A <- meta[,1][meta[,2] == input$comparisonA2]
        B <- meta[,1][meta[,2] == input$comparisonB2]
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
        genes <- rownames(top1)[which(top1$adj.P.Val < 0.05 & top1$logFC > 1)]
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
        top_probes <- input$NumFeatures
        col_labels <- colnames(expr)
        isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
        exp <- expr[isexpr,]
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
        colnames(zdataset) <- names(dataset)
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
        A <- meta[,1][meta[,2] == input$comparisonA2]
        B <- meta[,1][meta[,2] == input$comparisonB2]
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
        genes <- rownames(top1)[which(top1$adj.P.Val < 0.05 & top1$logFC < -1)]
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
    
    #render cluster assignment data table
    output$Clusters <- DT::renderDataTable({
        top_probes <- input$NumFeatures
        col_labels <- colnames(expr)
        isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
        exp <- expr[isexpr,]
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
        colnames(zdataset) <- names(dataset)
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
        DT::datatable(output, options = list(lengthMenu = c(5,10,20,50,100,1000), pageLength = 5, scrollX = TRUE),
                      selection=list(mode = "multiple", selected = c(1)))
    })
    
    #render gene list table for boxplot selection
    output$GeneListTable <- DT::renderDataTable({
        DT::datatable(geneList,
                      selection = 'single',
                      options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100")))
    })
    
    #render gene scatter plot data table
    output$geneScatterTable <- DT::renderDataTable({
        
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
                rownames(expr_t2) == meta[,1] ~ meta[,2],
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
        
        #make table
        Sample <- rownames(expr_t3)
        Type <- expr_t3[,'type']
        gene1col <- expr_t3[,gene1.u]
        gene2col <- expr_t3[,gene2.u]
        scatterTab <- data.frame(Sample, Type, gene1col, gene2col)
        colnames(scatterTab)[c(3,4)] <- c(gene1.n, gene2.n)
        
        #table output
        DT::datatable(scatterTab,
                      options = list(keys = TRUE,
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = c("10", "25", "50", "100")))
        
    })
    
    #render ssGSEA table
    output$ssGSEAtable <- DT::renderDataTable({
        if (input$tables == 3) {
            if (length(input$tab2table_rows_selected) > 0){
                ssgsea <- ssGSEAfunc()
                ssgsea2 <- as.data.frame(t(ssgsea))
                samporder <- meta[,1]
                ssgsea3 <- as.data.frame(ssgsea2[samporder,])
                colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
                rownames(ssgsea3) <- samporder
                ssgsea4 <- ssgsea3 %>% 
                    mutate(type = case_when(
                        rownames(ssgsea3) == meta[,1] ~ meta[,2],
                    ))
                ssgsea4 <- ssgsea4 %>%
                    relocate(type)
                #table output
                DT::datatable(ssgsea4,
                              options = list(keys = TRUE,
                                             searchHighlight = TRUE,
                                             pageLength = 10,
                                             lengthMenu = c("10", "25", "50", "100")))
            }
        }
        else if (input$tables == 1) {
            if (length(input$msigdbTable_rows_selected) > 0){
                ssgsea <- ssGSEAfunc()
                ssgsea2 <- as.data.frame(t(ssgsea))
                samporder <- meta[,1]
                ssgsea3 <- as.data.frame(ssgsea2[samporder,])
                colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
                rownames(ssgsea3) <- samporder
                ssgsea4 <- ssgsea3 %>% 
                    mutate(type = case_when(
                        rownames(ssgsea3) == meta[,1] ~ meta[,2],
                    ))
                ssgsea4 <- ssgsea4 %>%
                    relocate(type)
                #table output
                DT::datatable(ssgsea4,
                              options = list(keys = TRUE,
                                             searchHighlight = TRUE,
                                             pageLength = 10,
                                             lengthMenu = c("10", "25", "50", "100")))
            }
        }
        
    })
    
    
    ####----Plots----####
    
    
    #render GSEA plot
    output$enrichplot0 <- renderPlot({
        if (input$tables == 1){
            if (length(input$msigdbTable_rows_selected) > 0){
                res <- datasetInput()
                gsea.df <- as.data.frame(res@result)
                geneset <- msigdb.gsea2[input$msigdbTable_rows_selected,3]
                title <- geneset
                gseaplot2(res,
                          geneset,
                          title,
                          pvalue_table = F)
            }
        }
        else if (input$tables == 3){
            if (length(input$tab2table_rows_selected) > 0){
                res <- datasetInput()
                gsea.df <- as.data.frame(res@result)
                geneset <- as.character(GeneSet2[input$tab2table_rows_selected,1])
                title <- geneset
                gseaplot2(res,
                          geneset,
                          title,
                          pvalue_table = F)
            }
        }
        else if (input$tables == 5){
            if (length(input$GStable.u_rows_selected) > 0){
                res <- datasetInput()
                gsea.df <- as.data.frame(res@result)
                geneset <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
                title <- geneset
                gseaplot2(res,
                          geneset,
                          title,
                          pvalue_table = F)
            }
        }
    })
    
    #render heatmap
    output$heatmap0 <- renderPlot({
        if (input$tables == 1) {
            if (length(input$msigdbTable_rows_selected) > 0){
                groupA <- meta[,1][meta[,2] == input$comparisonA]
                groupB <- meta[,1][meta[,2] == input$comparisonB]
                res <- datasetInput()
                gsea.df <- as.data.frame(res@result)
                GS <- msigdb.gsea2[input$msigdbTable_rows_selected,3]
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
                hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
                pheatmap(dataset,            #data
                         cluster_cols = F,    #cluster columns - NO
                         cluster_row = F,     #cluster rows - YES
                         fontsize_col = 12,   #column fontsize
                         show_rownames = T,  
                         show_colnames = T,
                         color=hmcols,
                         annotation_col = meta3)
            }
        }
        else if (input$tables == 3) {
            if (length(input$tab2table_rows_selected) > 0){
                groupA <- meta[,1][meta[,2] == input$comparisonA]
                groupB <- meta[,1][meta[,2] == input$comparisonB]
                res <- datasetInput()
                gsea.df <- as.data.frame(res@result)
                GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
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
                hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
                pheatmap(dataset,            #data
                         cluster_cols = F,    #cluster columns - NO
                         cluster_row = F,     #cluster rows - YES
                         fontsize_col = 12,   #column fontsize
                         show_rownames = T,  
                         show_colnames = T,
                         color=hmcols,
                         annotation_col = meta3)
            }
        }
        else if (input$tables == 5) {
            if (length(input$GStable.u_rows_selected) > 0){
                groupA <- meta[,1][meta[,2] == input$comparisonA]
                groupB <- meta[,1][meta[,2] == input$comparisonB]
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
                hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
                pheatmap(dataset,            #data
                         cluster_cols = F,    #cluster columns - NO
                         cluster_row = F,     #cluster rows - YES
                         fontsize_col = 12,   #column fontsize
                         show_rownames = T,  
                         show_colnames = T,
                         color=hmcols,
                         annotation_col = meta3)
            }
        }
    })
    
    #render RNAseq heatmap
    output$heatmap1 <- renderPlot({
        if (input$customs == 444){
            if (length(input$heatmapGeneSelec) >= 2 || length(input$userheatgenes) >= 1) {
                genelist.uih <- NULL
                genelist.ush <- NULL
                genelist.uih2 <- NULL
                genelist.ush <- input$heatmapGeneSelec
                genelist.uih <- unlist(strsplit(input$userheatgenes, " "))
                genelist.uih2 <- unlist(strsplit(input$userheatgenes, "\t"))
                #usersamps <- unlist(strsplit(input$userheatsamp, " "))
                usersamps <- input$userheatsamp2
                exp <- expr[,usersamps]
                meta <- meta[which(meta[,1] %in% usersamps),]
                mad <- NULL
                var <- NULL
                cv <- NULL
                var_type <- input$VarianceMeasure
                if (var_type == "MAD"){
                    mad <- apply(log2(exp + 1), 1, mad)
                    mad <- sort(mad, decreasing = T)
                    mad1 <- mad[which(names(mad) %in% genelist.ush)]
                    mad2 <- mad[which(names(mad) %in% genelist.uih)]
                    mad3 <- mad[which(names(mad) %in% genelist.uih2)]
                    mad <- c(mad1,mad2,mad3)
                    out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
                    colnames(out) <- c("Gene", "MAD", colnames(exp))
                    dataset <- exp[names(mad),]
                    variable_gene_list <- names(mad)
                }
                if (var_type == "VAR"){
                    var <- apply(log2(exp + 1), 1, var)
                    var <- sort(var, decreasing = T)
                    var1 <- var[which(names(var) %in% genelist.ush)]
                    var2 <- var[which(names(var) %in% genelist.uih)]
                    mad3 <- mad[which(names(mad) %in% genelist.uih2)]
                    mad <- c(mad1,mad2,mad3)
                    var <- c(var1,var2)
                    out <- cbind(names(var), var[names(var)], exp[names(var),])
                    colnames(out) <- c("Gene", "VAR", colnames(exp))
                    dataset <- exp[names(var),]
                    variable_gene_list <- names(var)
                }
                if (var_type == "CV"){
                    cv <- apply(log2(exp + 1), 1, cv)
                    cv <- sort(cv, decreasing = T)
                    cv1 <- cv[which(names(cv) %in% genelist.ush)]
                    cv2 <- cv[which(names(cv) %in% genelist.uih)]
                    mad3 <- mad[which(names(mad) %in% genelist.uih2)]
                    mad <- c(mad1,mad2,mad3)
                    cv <- c(cv1,cv2)
                    out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
                    colnames(out) <- c("Gene", "VAR", colnames(exp))
                    dataset <- exp[names(cv),]
                    variable_gene_list <- names(cv)
                }
                
                dataset <- log2(dataset + 1)
                zdataset <- apply(dataset, 1, scale)
                zdataset <- apply(zdataset, 1, rev)
                colnames(zdataset) <- names(dataset)
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
                meta <- meta[order(meta[,2]),]
                type <- meta[,2]
                meta2 <- as.data.frame(type)
                rownames(meta2) <- meta[,1]
                bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
                hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
                pheatmap(dataset,
                         cluster_col = F,
                         cluster_row = T,
                         fontsize_row = 9,
                         fontsize_col = 10,
                         show_rownames = T ,
                         show_colnames = T,
                         annotation_col = meta2,
                         clustering_method = input$ClusteringMethod,
                         color=hmcols,
                         border_color = NA)
            }
        }
        else {
            top_probes <- input$NumFeatures
            exp <- expr
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
            colnames(zdataset) <- names(dataset)
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
            type <- meta[,2]
            meta2 <- as.data.frame(type)
            rownames(meta2) <- meta[,1]
            bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
            hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
            pheatmap(dataset,
                     cluster_col = T,
                     cluster_row = T,
                     fontsize_row = 9,
                     fontsize_col = 10,
                     show_rownames = T ,
                     annotation_col = meta2,
                     show_colnames = T,
                     clustering_method = input$ClusteringMethod,
                     color=hmcols,
                     border_color = NA)
        }
        
    })
    
    #render MA plot
    output$MAPlot1 <- renderPlot({
        top2 <- topgenereact()
        top_hits_up <- top2[head(which(top2$logFC > abs(input$fc_cutoff)), n = input$top_x),]
        top_hits_dn <- top2[head(which(top2$logFC < -abs(input$fc_cutoff)), n = input$top_x),]
        Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
        x <- ggplot(data = top2, aes(x = AveExpr, y = logFC)) +
            geom_point(size = 2, shape = 16) +
            theme_light(base_size = 16)
        x <- x + aes(color = group2) +
            scale_color_manual(values = Okabe_Ito)
        x <- x + geom_hline(yintercept = c(-abs(input$fc_cutoff),abs(input$fc_cutoff)), linetype="dashed", color="gray20")
        x <- x + geom_text_repel(
            data =  top_hits_up,
            aes(label = rownames(top_hits_up)),
            color="gray20",
            size = 6,
            nudge_x = 0.2,
            nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = 50)
        x <- x + geom_text_repel(
            data =  top_hits_dn,
            aes(label = rownames(top_hits_dn)),
            color="gray20",
            size = 6,
            nudge_x = 0.2,
            nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = 50)
        top2_selec <- top2 %>%
            filter(GeneName %in% input$userGeneSelec)
        x <- x + geom_text_repel(
            data =  top2_selec,
            aes(label = rownames(top2_selec)),
            color="gray20",
            size = 6,
            nudge_x = 0.2,
            nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = 50)
        x <- x + theme(legend.position="none")
        x <- x + labs(x = "Average Expression")
        x <- x + theme(axis.text = element_text(size=18))
        x <- x + theme(axis.title = element_text(size=24))
        x
    })
    
    #render boxplot
    output$boxplot1 <- renderPlot({
        if (length(input$GeneListTable_rows_selected) > 0){
            top_probes <- input$NumFeatures
            col_labels <- colnames(expr)
            isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
            exp <- expr[isexpr,]
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
            colnames(zdataset) <- names(dataset)
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
            bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
            hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
            results2 <- hclust(dist(t(dataset)), method = input$ClusteringMethod)
            m <- sort(cutree(results2, k = input$NumClusters))
            output <- cbind(colnames(m), as.matrix(m))
            colnames(output) <- c("Cluster")
            
            rownames(output) <- rownames(as.matrix(m))
            gene <- geneList[input$GeneListTable_rows_selected, 1]
            min <- min(log2(expr[gene,] + 1.0))
            max <- max(log2(expr[gene,] + 1.0))
            meta_temp <- meta
            rownames(meta_temp) <- meta[,1]
            meta_temp <- meta_temp %>%
                select(Group)
            data = merge(t(expr[gene,]), meta_temp, by=0)
            colnames(data) = c("SampleName", "GeneExpr", "Cluster")
            ggplot(data, aes(x=Cluster, y=log2(GeneExpr + 1.0))) +
                geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=1) +
                geom_dotplot(binaxis='y', stackdir='center', dotsize=0.75) +
                stat_compare_means(label = "p.signif") + ylim(min * 0.9, max * 1.3) +
                theme_bw() +
                labs(title= paste(gene, "Expression (log2)")) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
        }
    })
    
    #render up regulated pathway enrichment plot
    output$UpRegPathway1 <- renderPlot({
        top1 <- topgenereact()
        genes <- rownames(top1)[which(top1$adj.P.Val < 0.05 & top1$logFC > 1)]
        dbs <- listEnrichrDbs() 
        enrichRLive <- TRUE 
        if (is.null(dbs)) { 
            enrichRLive <- FALSE 
        }
        dbs <- input$SelectedPathway
        enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value 
        plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value") 
    })
    
    #render down regulated pathway enrichment plot
    output$DnRegPathway1 <- renderPlot({
        top1 <-top2 <- topgenereact()
        genes <- rownames(top1)[which(top1$adj.P.Val < 0.05 & top1$logFC < -1)]
        dbs <- listEnrichrDbs() 
        enrichRLive <- TRUE 
        if (is.null(dbs)) { 
            enrichRLive <- FALSE 
        }
        dbs <- input$SelectedPathway
        enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value 
        plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value") 
    })
    
    #render volcano plot
    output$Volcano3 <- renderPlot({
        top2 <- topgenereact()
        top_hits_up <- top2[head(which(top2$logFC > abs(input$fc_cutoff) & top2$P.Value < input$p_cutoff), n = input$top_x),]
        top_hits_dn <- top2[head(which(top2$logFC < -abs(input$fc_cutoff) & top2$P.Value < input$p_cutoff), n = input$top_x),]
        Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
        x <- ggplot(data = top2, aes(x = logFC, y = -log10(P.Value))) +
            geom_point(size = 2, shape = 16) +
            theme_light(base_size = 16)
        x <- x + aes(color = group) +
            scale_color_manual(values = Okabe_Ito)
        x <- x + geom_vline(xintercept = c(-abs(input$fc_cutoff),abs(input$fc_cutoff)), linetype="dashed", color="gray20")
        x <- x + geom_hline(yintercept = -log10(input$p_cutoff), linetype="dashed", color="gray20")
        x <- x + geom_text_repel(
            data =  top_hits_up,
            aes(label = rownames(top_hits_up)),
            size = 6,
            color="gray20",
            nudge_x = 0.2,
            nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = 50)
        x <- x + geom_text_repel(
            data =  top_hits_dn,
            aes(label = rownames(top_hits_dn)),
            size = 6,
            color="gray20",
            nudge_x = 0.2,
            nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = 50)
        top2_selec <- top2 %>%
            filter(GeneName %in% input$userGeneSelec)
        x <- x + geom_text_repel(
            data =  top2_selec,
            aes(label = rownames(top2_selec)),
            size = 6,
            color="gray20",
            nudge_x = 0.2,
            nudge_y=0.2,
            box.padding = unit(0.9, "lines"),
            point.padding = unit(.3+4*0.1, "lines"),
            max.overlaps = 50)
        x <- x + theme(legend.position="none")
        x <- x + theme(axis.text = element_text(size=18))
        x <- x + theme(axis.title = element_text(size=24))
        x
    })
    
    #render ssGSEA boxplot
    output$boxplot2 <- renderPlot({
        if (input$tables == 3) {
            if (length(input$tab2table_rows_selected) > 0){
                ssgsea <- ssGSEAfunc()
                ssgsea2 <- as.data.frame(t(ssgsea))
                samporder <- meta[,1]
                ssgsea3 <- as.data.frame(ssgsea2[samporder,])
                colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
                rownames(ssgsea3) <- samporder
                ssgsea4 <- ssgsea3 %>% 
                    mutate(type = case_when(
                        rownames(ssgsea3) == meta[,1] ~ meta[,2],
                    ))
                ssgsea4 <- ssgsea4 %>%
                    relocate(type)
                if (input$boxplotcompare == "wilcox.test"){
                    ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
                        geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
                        geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 0.75) +
                        labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                             title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
                        theme(plot.title = element_text(size=12),
                              axis.text.x = element_text(size=12),
                              axis.text.y = element_text(size=12),
                              axis.title = element_text(size=12),
                              legend.text = element_text(size=12),
                              legend.title = element_text(size=12)) +
                        theme_bw() +
                        stat_compare_means() +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
                }
                else if (input$boxplotcompare == "t.test"){
                    ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
                        geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
                        geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 0.75) +
                        labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                             title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
                        theme(plot.title=element_text(size=12),
                              axis.text.x = element_text(size=12),
                              axis.text.y = element_text(size=12),
                              axis.title = element_text(size=12),
                              legend.text = element_text(size=12),
                              legend.title = element_text(size=12)) +
                        theme_bw() +
                        stat_compare_means(method = "t.test") +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
                }
                else if(input$boxplotcompare == "kruskal.test"){
                    ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
                        geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
                        geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 0.75) +
                        labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                             title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
                        theme(plot.title=element_text(size=12),
                              axis.text.x = element_text(size=12),
                              axis.text.y = element_text(size=12),
                              axis.title = element_text(size=12),
                              legend.text = element_text(size=12),
                              legend.title = element_text(size=12)) +
                        theme_bw() +
                        stat_compare_means(method = "kruskal.test") +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
                }
                else if(input$boxplotcompare == "anova"){
                    ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
                        geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
                        geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 0.75) +
                        labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                             title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
                        theme(plot.title=element_text(size=12),
                              axis.text.x = element_text(size=12),
                              axis.text.y = element_text(size=12),
                              axis.title = element_text(size=12),
                              legend.text = element_text(size=12),
                              legend.title = element_text(size=12)) +
                        theme_bw() +
                        stat_compare_means(method = "anova") +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
                }
                else if(input$boxplotcompare == "none"){
                    ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
                        geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
                        geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 0.75) +
                        labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                             title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
                        theme(plot.title=element_text(size=12),
                              axis.text.x = element_text(size=12),
                              axis.text.y = element_text(size=12),
                              axis.title = element_text(size=12),
                              legend.text = element_text(size=12),
                              legend.title = element_text(size=12)) +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
                }
            }
        }
        else if (input$tables == 1) {
            if (length(input$msigdbTable_rows_selected) > 0){
                ssgsea <- ssGSEAfunc()
                ssgsea2 <- as.data.frame(t(ssgsea))
                samporder <- meta[,1]
                ssgsea3 <- as.data.frame(ssgsea2[samporder,])
                colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
                rownames(ssgsea3) <- samporder
                ssgsea4 <- ssgsea3 %>% 
                    mutate(type = case_when(
                        rownames(ssgsea3) == meta[,1] ~ meta[,2],
                    ))
                ssgsea4 <- ssgsea4 %>%
                    relocate(type)
                if (input$boxplotcompare == "wilcox.test"){
                    ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
                        geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
                        geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 0.75) +
                        labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                             title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
                        theme(plot.title=element_text(size=12),
                              axis.text.x = element_text(size=12),
                              axis.text.y = element_text(size=12),
                              axis.title = element_text(size=12),
                              legend.text = element_text(size=12),
                              legend.title = element_text(size=12)) +
                        theme_bw() +
                        stat_compare_means() +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
                }
                else if (input$boxplotcompare == "t.test"){
                    ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
                        geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
                        geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 0.75) +
                        labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                             title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
                        theme(plot.title=element_text(size=12),
                              axis.text.x = element_text(size=12),
                              axis.text.y = element_text(size=12),
                              axis.title = element_text(size=12),
                              legend.text = element_text(size=12),
                              legend.title = element_text(size=12)) +
                        theme_bw() +
                        stat_compare_means(method = "t.test") +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
                }
                else if(input$boxplotcompare == "kruskal.test"){
                    ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
                        geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
                        geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 0.75) +
                        labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                             title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
                        theme(plot.title=element_text(size=12),
                              axis.text.x = element_text(size=12),
                              axis.text.y = element_text(size=12),
                              axis.title = element_text(size=12),
                              legend.text = element_text(size=12),
                              legend.title = element_text(size=12)) +
                        theme_bw() +
                        stat_compare_means(method = "kruskal.test") +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
                }
                else if(input$boxplotcompare == "anova"){
                    ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
                        geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
                        geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 0.75) +
                        labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                             title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
                        theme(plot.title=element_text(size=12),
                              axis.text.x = element_text(size=12),
                              axis.text.y = element_text(size=12),
                              axis.title = element_text(size=12),
                              legend.text = element_text(size=12),
                              legend.title = element_text(size=12)) +
                        theme_bw() +
                        stat_compare_means(method = "anova") +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
                }
                else if(input$boxplotcompare == "none"){
                    ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
                        geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
                        geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 0.75) +
                        labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                             title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
                        theme(plot.title=element_text(size=12),
                              axis.text.x = element_text(size=12),
                              axis.text.y = element_text(size=12),
                              axis.title = element_text(size=12),
                              legend.text = element_text(size=12),
                              legend.title = element_text(size=12)) +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
                }
            }
        }
    })
    
    #render gene expression comparison scatter plot
    output$geneScatter0 <- renderPlotly({
        
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
                rownames(expr_t2) == meta[,1] ~ meta[,2],
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
                                 color = type,
                                 text = paste("</br> Sample: ", rownames(expr_t3),
                                              "</br> Type: ", type))) +
            geom_point() +
            theme_minimal() +
            labs(x = gene1.n, y = gene2.n,
                 title = paste(gene1.n, " vs. ", gene2.n, sep = ''))
        
        ggplotly(p, tooltip = 'text')
        
    })
    
    
    ####----Download Handlers----####
    
    
    output$geneScatterDownload <- downloadHandler(
        filename = function() {
            #user gene input
            gene1.u <- input$scatterG1
            gene2.u <- input$scatterG2
            paste(gene1.u, "_vs_", gene2.u, "_Expression.tsv", sep = "")
        },
        content = function(file) {
            #transpose
            expr_t <- as.data.frame(t(expr))
            
            #reorder rowname to match meta for merging
            samporder <- meta[,1]
            expr_t2 <- as.data.frame(expr_t[samporder,])
            
            #add type
            expr_t3 <- expr_t2 %>% 
                mutate(type = case_when(
                    rownames(expr_t2) == meta[,1] ~ meta[,2],
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
    
    output$ssGSEAdownload <- downloadHandler(
        filename = function() {
            if (input$tables == 1) {
                paste('ssGSEAscore_',names(gs[(msigdb.gsea2[input$msigdbTable_rows_selected,3])]),'.tsv',sep = '')
            }
            else if (input$tables == 3) {
                paste('ssGSEAscore_',names(gs2[(GeneSet2[input$tab2table_rows_selected,1])]),'.tsv',sep = '')
            }
        },
        content = function(file) {
            if (input$tables == 1) {
                GS <- gs[(msigdb.gsea2[input$msigdbTable_rows_selected,3])]
            }
            else if (input$tables == 2) {
                GS <- gs2[(GeneSet2[input$tab2table_rows_selected,1])]
            }
            ssgsea <- ssGSEAfunc()
            ssgsea2 <- as.data.frame(t(ssgsea))
            samporder <- meta[,1]
            ssgsea3 <- as.data.frame(ssgsea2[samporder,])
            colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
            rownames(ssgsea3) <- samporder
            ssgsea4 <- ssgsea3 %>% 
                mutate(type = case_when(
                    rownames(ssgsea3) == meta[,1] ~ meta[,2],
                ))
            ssgsea4 <- ssgsea4 %>%
                relocate(type)
            ssgsea4$sample <- rownames(ssgsea4)
            ssgsea5 <- ssgsea4[,c(3,1,2)]
            rownames(ssgsea5) <- 1:nrow(ssgsea5)
            write_tsv(ssgsea5, file)
        }
    )
    
    #render download button for DEG GMT
    output$DEGgmtDownload <- downloadHandler(
        filename = function() {
            paste(input$DEGfileName,".gmt",sep = "")
        },
        content = function(file) {
            top_probes <- input$NumFeatures
            col_labels <- colnames(expr)
            isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
            exp <- expr[isexpr,]
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
            colnames(zdataset) <- names(dataset)
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
            A <- meta[,1][meta[,2] == input$comparisonA2]
            B <- meta[,1][meta[,2] == input$comparisonB2]
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
            if (input$UpDnChoice == "UpAndDown_Regulated"){
                genes <- rownames(top1)[which(abs(top1$logFC) > abs(input$fc_cutoff2) & top1$adj.P.Val < input$p_cutoff2)]
                genes.h <- t(as.data.frame(head(genes, n=input$top_x2)))
                genes.h.gmt <- data.frame(paste(input$top_x2,input$UpDnChoice,"DEgenes", sep = ""),
                                          paste(input$UpDnChoice,"DEgenes",sep = ""),
                                          genes.h)
            }
            else if (input$UpDnChoice == "Up_Regulated"){
                genes <- rownames(top1)[which(top1$logFC > input$fc_cutoff2 & top1$adj.P.Val < input$p_cutoff2)]
                genes.h <- t(as.data.frame(head(genes, n=input$top_x2)))
                genes.h.gmt <- data.frame(paste(input$top_x2,input$UpDnChoice,"DEgenes", sep = ""),
                                          paste(input$UpDnChoice,"DEgenes",sep = ""),
                                          genes.h)
            }
            else if (input$UpDnChoice == "Down_Regulated"){
                genes <- rownames(top1)[which(top1$logFC < -abs(input$fc_cutoff2) & top1$adj.P.Val < input$p_cutoff2)]
                genes.h <- t(as.data.frame(head(genes, n=input$top_x2)))
                genes.h.gmt <- data.frame(paste(input$top_x2,input$UpDnChoice,"DEgenes", sep = ""),
                                          paste(input$UpDnChoice,"DEgenes",sep = ""),
                                          genes.h)
            }
            write_delim(genes.h.gmt, file, delim = '\t', col_names = F)
        }
    )
    
    #render download button for upreg gene pathways GMT
    output$UpRegPathDownloadgmt <- downloadHandler(
        filename = function() {
            paste("UpReg_",input$SelectedPathway,"_pathway.gmt",sep = "")
        },
        content = function(file) {
            top1 <- topgenereact()
            genes <- rownames(top1)[which(top1$adj.P.Val < 0.05 & top1$logFC > 1)]
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
            genes <- rownames(top1)[which(top1$adj.P.Val < 0.05 & top1$logFC > 1)]
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
            genes <- rownames(top1)[which(top1$adj.P.Val < 0.05 & top1$logFC > 1)]
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
            genes <- rownames(top1)[which(top1$adj.P.Val < 0.05 & top1$logFC < -1)]
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
            col_labels <- colnames(expr)
            isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
            exp <- expr[isexpr,]
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
    
    #render Most Variable Gene download button
    output$MVGdownloadgmt <- downloadHandler(
        filename = function() {
            top_probes <- input$NumFeatures
            paste(top_probes,"_Most_Variable_Genes", ".gmt", sep = "")
        },
        content = function(file){
            top_probes <- input$NumFeatures
            col_labels <- colnames(expr)
            isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
            exp <- expr[isexpr,]
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
                                 "Most_Variable_Genes",genes.n)
            rownames(gmt.df) <- NULL
            write_delim(gmt.df, file, delim = '\t', col_names = F)
        }
    )
    
    #download button for leading edge genes
    output$LEGdownload <- downloadHandler(
        filename = function() {
            if (input$tables == 1){
                if (length(input$msigdbTable_rows_selected) > 0){
                    GS <- msigdb.gsea2[input$msigdbTable_rows_selected,3]
                }
            }
            else if (input$tables == 3){
                if (length(input$tab2table_rows_selected) > 0){
                    GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
                }
            }
            else if (input$tables == 5){
                if (length(input$GStable.u_rows_selected) > 0){
                    GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
                }
            }
            paste(GS,"_leading_edge_genes", ".tsv", sep = "")
        },
        content = function(file){
            if (input$tables == 1){
                if (length(input$msigdbTable_rows_selected) > 0){
                    GS <- msigdb.gsea2[input$msigdbTable_rows_selected,3]
                }
            }
            else if (input$tables == 3){
                if (length(input$tab2table_rows_selected) > 0){
                    GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
                }
            }
            else if (input$tables == 5){
                if (length(input$GStable.u_rows_selected) > 0){
                    GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
                }
            }
            res <- datasetInput()
            gsea.df <- as.data.frame(res@result)
            ## Subset core enriched genes
            genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
            genes2 <- strsplit(genes1,"/")
            GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
            GeneSymbol$Rank <- rownames(GeneSymbol)
            GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
            write_delim(GeneSymbol, file, delim = '\t')
        }
    )
    
    #render DEG table download button
    output$DEGtableDownload <- downloadHandler(
        filename = function() {
            paste("DEG_Table_",Sys.Date(), ".tsv", sep = "")
        },
        content = function(file){
            top_probes <- input$NumFeatures
            col_labels <- colnames(expr)
            isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
            exp <- expr[isexpr,]
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
            colnames(zdataset) <- names(dataset)
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
            A <- meta[,1][meta[,2] == input$comparisonA2]
            B <- meta[,1][meta[,2] == input$comparisonB2]
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
            top1$Gene <- rownames(top1)
            top1 <- top1 %>%
                relocate(Gene)
            write_delim(top1, file, delim = '\t')
        }
    )
    
    #download button for leading edge genes
    output$enrich_sig_download <- downloadHandler(
        filename = function() {
            paste("Enriched_Signatures_Table_",Sys.Date(), ".tsv", sep = "")
        },
        content = function(file){
            if (input$SigTableChoice == SigNames[1]) {
                gsea.df <- as_tibble(ES_table1)
                DT::datatable(gsea.df,
                              extensions = c("KeyTable", "FixedHeader"),
                              caption = "Enriched Signatures",
                              options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"),scrollX = T)) %>%
                    formatRound(columns = c(2:10), digits = 2)
            }
            else if (input$SigTableChoice == SigNames[2]) {
                gsea.df <- as_tibble(ES_table2)
                DT::datatable(gsea.df,
                              extensions = c("KeyTable", "FixedHeader"),
                              caption = "Enriched Signatures",
                              options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"),scrollX = T)) %>%
                    formatRound(columns = c(2:10), digits = 2)
            }
            else if (input$SigTableChoice == SigNames[3]) {
                gsea.df <- as_tibble(ES_table3)
                DT::datatable(gsea.df,
                              extensions = c("KeyTable", "FixedHeader"),
                              caption = "Enriched Signatures",
                              options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"),scrollX = T)) %>%
                    formatRound(columns = c(2:10), digits = 2)
            }
            
            enrich_sig_table <- ES.df
            write_delim(enrich_sig_table, file, delim = '\t')
        }
    )
    
    #download button for user enriched signature table
    output$enrich_sig_download.u <- downloadHandler(
        filename = function() {
            paste("Enriched_Signatures_Table_",Sys.Date(), ".tsv", sep = "")
        },
        content = function(file) {
            if (input$tables == 3) {
                groupA <- meta[,1][meta[,2] == input$comparisonA]
                groupB <- meta[,1][meta[,2] == input$comparisonB]
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
                s2n.df <- as.data.frame(s2n.matrix)
                s2n.df$GeneID <- rownames(s2n.df)
                rownames(s2n.df) <- NULL
                data <- dplyr::select(s2n.df, GeneID, V1)
                data.gsea <- data$V1
                names(data.gsea) <- as.character(data$GeneID)
                s2n.matrix.s <- sort(data.gsea, decreasing = T)
                ##----GSEA----##
                gmt.i <- tab2
                gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
                gsea.df <- as.data.frame(gsea.res@result)
                write_delim(gsea.df, file, delim = '\t')
            }
            else if (input$tables == 5) {
                groupA <- meta[,1][meta[,2] == input$comparisonA]
                groupB <- meta[,1][meta[,2] == input$comparisonB]
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
                s2n.df <- as.data.frame(s2n.matrix)
                s2n.df$GeneID <- rownames(s2n.df)
                rownames(s2n.df) <- NULL
                data <- dplyr::select(s2n.df, GeneID, V1)
                data.gsea <- data$V1
                names(data.gsea) <- as.character(data$GeneID)
                s2n.matrix.s <- sort(data.gsea, decreasing = T)
                ##----GSEA----##
                gmt.i <- GStable.ubg()
                gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
                gsea.df <- as.data.frame(gsea.res@result)
                write_delim(gsea.df, file, delim = '\t')
            }
        }
    )
    
    
    ####----Text----####
    
    
    #NES and Pval output
    output$NESandPval <- renderText({
        if (input$tables == 1){
            if (length(input$msigdbTable_rows_selected) > 0){
                res <- datasetInput()
                gsea.df <- as.data.frame(res@result)
                GS = msigdb.gsea2[input$msigdbTable_rows_selected,3]
                NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
                Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
                NES.o <- paste0("NES: ", NES)
                Pval.o <- paste0("Pvalue: ", Pval)
                if (NES > 0){
                    UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", input$comparisonA , "group.")
                }
                else {
                    UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", input$comparisonA , "group.")
                }
                paste(NES.o, Pval.o, UpOrDown, sep = '\n')
            }
        }
        else if (input$tables == 3){
            if (length(input$tab2table_rows_selected) > 0){
                res <- datasetInput()
                gsea.df <- as.data.frame(res@result)
                GS = as.character(GeneSet2[input$tab2table_rows_selected,1])
                NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
                Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
                NES.o <- paste0("NES: ", NES)
                Pval.o <- paste0("Pvalue: ", Pval)
                if (NES > 0){
                    UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", input$comparisonA , "group.")
                }
                else {
                    UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", input$comparisonA , "group.")
                }
                paste(NES.o, Pval.o, UpOrDown, sep = '\n')
            }
        }
        else if (input$tables == 5){
            if (length(input$GStable.u_rows_selected) > 0){
                res <- datasetInput()
                gsea.df <- as.data.frame(res@result)
                GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
                NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
                Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
                NES.o <- paste0("NES: ", NES)
                Pval.o <- paste0("Pvalue: ", Pval)
                if (NES > 0){
                    UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", input$comparisonA , "group.")
                }
                else {
                    UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", input$comparisonA , "group.")
                }
                paste(NES.o, Pval.o, UpOrDown, sep = '\n')
            }
        }
    })
    
    output$VolGroupsText <- renderText({
        paste("This volcano plot is comparing group A: ",input$comparisonA2, " and group B: ",input$comparisonB2, ".\nGenes with a positive log fold change are upregulated in the ",
              input$comparisonA2, " group.\nGenes with a negative log fold change are upregulated in the ",input$comparisonB2, " group.", sep = "")
    })
    
    output$MAGroupsText <- renderText({
        paste("This MA plot is comparing group A: ",input$comparisonA2, " and group B: ",input$comparisonB2, ".\nGenes with a positive log fold change are upregulated in the ",
              input$comparisonA2, " group.\nGenes with a negative log fold change are upregulated in the ",input$comparisonB2, " group.", sep = "")
    })
    
    output$upregpath_text <- renderText({
        paste("Genes in these enriched terms are upregulated in group A: ", input$comparisonA2, " group.", sep = "")
    })
    
    output$downregpath_text <- renderText({
        paste("Genes in these enriched terms are upregulated in group B: ", input$comparisonB2," group", sep = "")
    })
    
    output$degtext <- renderText({
        paste("This table represents differentially expressed genes when comparing group A: ",input$comparisonA2," and group B: ",input$comparisonB2,
              ".\nGenes with a positive logFC, indicate an upregulation in group A: ", input$comparisonA2,
              ".\nGenes with a negative logFC, indicate an upregulation in group B: ", input$comparisonB2,
              ".\nThe 'AveExpr' column represents the log transformed average expression between group A: ",input$comparisonA2," and group B: ",input$comparisonB2,".", sep = "")
    })
}




# Run the application 
shinyApp(ui = ui, server = server)
