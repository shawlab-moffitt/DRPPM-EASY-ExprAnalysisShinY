

####----Install and load packages----####

packages <- c("clusterProfiler","dplyr","DT","msigdbr","GSVA","readr","ggpubr","reshape2")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


####-------------------------------------------------------------####
#                                                                   #
#  This is a script to generate Enriched Signature Table(s) for the #
#  GSEA portion of the General Expression Analysis R Shiny App      #
#                                                                   #
####-------------------------------------------------------------####


####----Files----####

##User Input

#Assign expression data path
expr_file <- "~/R/GeneralExpressionAnalysisApp-main/ExampleData/htseq_gene_level_fpkm_T_geneName_max_1cutoff_v2.txt"

#Assign meta data path
meta_file <- "~/R/GeneralExpressionAnalysisApp-main/ExampleData/USP7_meta.tsv"
#Is there a header?
header <- TRUE

#Gene Set File
GeneSet_file <- "~/R/GeneralExpressionAnalysisApp-main/GeneSets/msigdb_gsNsym_HS.zip"

#assign outfile path - where you want the enriched signatures table
OutPath <- "~/R/GeneralExpressionAnalysisApp-main/ExampleData/"





#Code to generate gene set for MSigDB based on species if not using human or mouse
#msigdbr_species() #to see species available
#msigdb_gs <- msigdbr(species = "Homo sapiens") %>%
#  dplyr::select(gs_name, gene_symbol)
#gmt <- msigdb_gs



####----Reading Files----####

##reading expression data
expr <- read.delim(expr_file, sep = '\t', header = T, strip.white = T)
colnames(expr)[1] <- "Gene"
expr <- expr %>%
  drop_na(Gene)
row.names(expr) <- expr[,1]
expr <- expr[,-1]
colnames(expr) <- gsub("[_.-]", "_", colnames(expr))
A <- as.matrix(expr)

##reading meta data
#turn header on or off if in file
meta <- read.delim(meta_file, sep = '\t', header = header, strip.white = T)
colnames(meta) <- c("SampleName","Group")
meta[,1] <- gsub("[_.-]", "_", meta[,1])



##reading GeneSet data

#generate gmt with pre-loaded MSigDB gmt file - Human or Mouse
if (file_ext(GeneSet_file) == "gmt") {
  gmt <- read.gmt(GeneSet_file)
}
if (file_ext(GeneSet_file) == "tsv" || file_ext(GeneSet_file) == "txt") {
  gmt <- read_delim(GeneSet_file, delim = '\t')
  colnames(gmt) <- c("Term","Gene")
}






####----Create groups based on meta----####
levlst <- list()
for (i in 1:length(levels(factor(meta[,2])))){
  group.i <- meta[,1][meta[,2] == levels(factor(meta[,2]))[i]]
  levlst <- append(levlst, list(group.i))
}
rm(group.i)
j <- 1
for (i in 1:length(levlst)){
  names(levlst)[i] <- paste("group",LETTERS[j],sep = "")
  j=j+1
}
list2env(levlst,globalenv())




####----Generate Enrich Sig Tables----####
##loop through combinations of the groups
##This will likely take several minutes due to the GSEA function

g_combos <- combn(names(levlst), 2)
for (i in 1:ncol(g_combos)){
  ####----Create groups----####
  #for calculations
  group1 <- unlist(levlst[which(names(levlst) == g_combos[1,i])])
  group2 <- unlist(levlst[which(names(levlst) == g_combos[2,i])])
  #for file names
  group1.n <- meta[,2][which(meta[,1] == group1[1])]
  group1.n <- gsub(" ","",group1.n)
  group2.n <- meta[,2][which(meta[,1] == group2[1])]
  group2.n <- gsub(" ","",group2.n)
  
  ####----Signal-to-Noise Calculation----####
  A <- A + 0.00000001
  P = as.matrix(as.numeric(colnames(A) %in% group1))
  n1 <- sum(P[,1])
  M1 <- A %*% P
  M1 <- M1/n1
  A2 <- A*A
  S1 <- A2 %*% P
  S1 <- S1/n1 - M1*M1 
  S1 <- sqrt(abs((n1/(n1-1)) * S1))
  P = as.matrix(as.numeric(colnames(A) %in% group2))
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
  
  ####----Reformatting----####
  test1 <- as.data.frame(s2n.matrix)
  test1$GeneID <- rownames(test1)
  rownames(test1) <- NULL
  data <- dplyr::select(test1, GeneID, V1)
  data.gsea <- data$V1
  names(data.gsea) <- as.character(data$GeneID)
  s2n.matrix.s <- sort(data.gsea, decreasing = T)
  
  ####----GSEA----####
  #perform GSEA
  gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, pvalueCutoff = 1)
  #extract results and convert to tibble
  gsea.df <- as_tibble(gsea.res@result)
  
  ####----Write to file----####
  write_delim(gsea.df, paste(OutPath,"Enrich_Sig_Table_",group1.n,"v",group2.n,".tsv",sep = ""), delim = "\t")
}

