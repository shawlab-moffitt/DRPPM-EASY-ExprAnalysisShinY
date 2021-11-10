

####-------------------------------------------------------------####
#                                                                   #
#  This is a script to easily generate an RData list for the        #
#  ssGSEA portion of the General Expression Analysis R Shiny App    #
#                                                                   #
####-------------------------------------------------------------####



####----User Input----####

#will accept .gmt, .tsv, or .txt
GeneSet_file <- '~/path/to/gene/set/file'
header.gs <- TRUE

#must use .RData extension
OutFile_PathAndName <- '~/path/to/desired/list.RData'




####----Read File----####

if (file_ext(GeneSet_file) == "gmt") {
  GeneSet <- read.gmt(GeneSet_file)
}
if (file_ext(GeneSet_file) == "tsv" || file_ext(GeneSet_file) == "txt") {
  GeneSet <- read.delim(GeneSet_file, header = header.gs, sep = '\t')
}

colnames(GeneSet) <- c("term","gene")



####----Generate RData List----####

gsDataList <- list()
for (i in unique(GeneSet[,1])){
  gsDataList[[i]] <- GeneSet[GeneSet[,1] == i,]$gene
}

save(gsDataList, file = OutFile_PathAndName)











