
library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(DT)
library(enrichplot)
library(ggplot2)
library(GSVA)
library(reshape2)
library(pheatmap)
library(ggpubr)
library(readr)



####----Files----####

##User Input

#Assign expression data
expr_file <- "~/R/PauloRodriguez/read_counts_normalized_geneSymbol_max_1FPKM.tsv"

#Assign meta data
meta_file <- "~/R/PauloRodriguez/metafile.txt"
#Is there a header?
header <- FALSE

#Assign .gmt file
#GeneSet_file <- "..."

#OR

#2 column file with GeneSetName and GeneSymbol as columns
GeneSet_file <- "~/R/USP7appsNstuff/data/MSigDB/msigdb_gsNsym_MS.tsv"


#Code to generate gene set for MSigDB based on species
#msigdbr_species() #to see species available
#msigdb_gs <- msigdbr(species = "Homo sapiens") %>%
#  dplyr::select(gs_name, gene_symbol)
#gmt <- msigdb_gs






####----Reading Files----####

##reading expression data
expr <- read.delim(expr_file, sep = '\t', header = T, row.names = 1, strip.white = T)
A <- as.matrix(expr)

##reading meta data
#turn header on or off if in file
meta <- read.delim(meta_file, sep = '\t', header = header, strip.white = T)

##reading GeneSet data
GeneSet <- read.delim(GeneSet_file, header = T, sep = '\t')
gmt <- GeneSet

##generate genelist data
#Gene <- rownames(expr)
#GeneList_MS <- as.data.frame(Gene)


##Create groups based off meta
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




####----Signal-to-Noise Calculation----####

#assign group from environment
groupA.r <- groupA
groupB.r <- groupB
#assign outfile name and path
OutPathAndName <- "~/R/PauloRodriguez/Enrich_Sig_Table"


A <- A + 0.00000001
P = as.matrix(as.numeric(colnames(A) %in% groupA.r))
n1 <- sum(P[,1])
M1 <- A %*% P
M1 <- M1/n1
A2 <- A*A
S1 <- A2 %*% P
S1 <- S1/n1 - M1*M1 
S1 <- sqrt(abs((n1/(n1-1)) * S1))
P = as.matrix(as.numeric(colnames(A) %in% groupB.r))
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

#Write table to file
write_delim(gsea.df, paste(OutPathAndName,".tsv",sep = ""), delim = "\t")




