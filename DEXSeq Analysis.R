#The most important thing about this program is having the right input files... These must be generated using python
#in the terminal. Install pysam for use in python and the terminal.
#First you need to determine (in R) the address of the python modules that you need:
pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
list.files(pythonScriptsDir)

#From here, refer to DEXSeq counts preparation file for the Command Line inputs to generate the appropriate exon counts.

setwd("~/Desktop/Bam_files")
#set the working directory to whichever folder contains all your aligned .txt files and their reference .gff
WD <- getwd()
countsFiles = list.files(WD, pattern="-scounts_1.txt$", full.names=TRUE)
#countsFiles
flattenedFile = list.files(WD, pattern="gff$", full.names=TRUE)
#flattenedFile
sampleTable = data.frame(row.names = c("WT1","WT2","SOD1","SOD2"), condition = c(rep("WT",2), rep("SOD", 2)), libType=c(rep("single-end",4)))
#change the above line according to the number of samples and the sample groups
suppressPackageStartupMessages(library("DEXSeq"))
dxd = DEXSeqDataSetFromHTSeq(countsFiles, sampleData=sampleTable, design=~sample + exon + condition:exon, flattened=flattenedFile)
colData(dxd)
#returns column info?
head(counts(dxd))
split (seq_len(ncol(dxd)), colData(dxd) $exon)
head(rowData(dxd))
 
sampleAnnotation(dxd)
dxd = estimateSizeFactors(dxd) #accounts for differences in library size
dxd = estimateDispersions(dxd) #calls DESeq2 to figure out dispersion across conditions for each exon. This one takes foreveeeeerrrrrrrr, maybe 45 min. or so. Faster if you can free up some memory
plotDispEsts(dxd) #Only works if you're running the code in a new session that has never seen this command before
 
dxd = testForDEU(dxd) 
dxd = estimateExonFoldChanges(dxd, fitExpToVar="sample") #This one also takes some time, maybe 20 min
dxr1 = DEXSeqResults(dxd)

table(dxr1$padj < 0.05) #tells you how many exons change usage with adjusted p-value less than 0.05.
 
DEXSeqHTML(dxr1, FDR=0.05, color.samples=c("#FF0000", "#00FF00", "#0000FF", "#000000"))
 +write.table(dxr1, file='DEXSeq_File.txt',sep="\t",row.names=FALSE,quote=FALSE)
plotDEXSeq(genename, "encode id", fitExpToVar="condition",splicing=TRUE,displayTranscripts=FALSE,legend=TRUE, names=TRUE, color.samples=c("#FF0000", "#00FF00", "#0000FF", "#000000"), cex.axis=1.2, cex=1.3, lwd=2)
