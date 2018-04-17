# This is a script for running edgeR.
# For help, type edgeRUsersGuide() after loading the library

library(edgeR)
library(Hmisc)
library(ggplot2)
library(reshape2)
library(Heatplus)
library(RColorBrewer)
library(pheatmap)

setwd("WD")

data <- read.delim("Counts.txt", header=TRUE)
data.matrix <- as.matrix(data)
rownames(data.matrix) <- data.matrix[,1]
data.matrix <- data.matrix[, 2:29] # This line and the one above move the gene name out of the first column and into the matrix data 
data.num <- apply(data.matrix[, 1:28], 2, as.numeric) # Makes the values numeric instead of string
rownames(data.num) <- rownames(data.matrix) # Changes the rownames again (which get reverted to integers after changing to numeric)

data.num.adj <- subset(data.num, select=-c(list_of_replicates_to_skip) #IF you want to remove any samples. Make sure to adjust 'group' accordingly

group <- factor(c(rep("Cell type 1", 3), rep("Cell type 2", 3), rep("Cell type 3", 3)))# Outlines the locations of each replicate in the individual matrices
group <- factor(c(rep('WT Pre', 3), rep('SOD Pre', 2), rep('WT Post', 2), rep('SOD Post', 3))) #If instead of cell types, you're comparing disease groups

## The following steps are for performing negative binomial analysis on the data:
y <- DGEList(counts=data.num.adj, group=group)
y <- calcNormFactors(y, method="TMM", doWeighting=TRUE)
##

## Clustering:
y$samples
mds <- plotMDS(y) # Plots the differential dispersions of the different conditions to see how closely they cluster to one another...

# These commands below will turn that clustering plot into something actually legible!:
clustering_data <- data.frame(mds$x, mds$y) #Takes the coordinates from the mds object created above and creates a plottable frame from it
ggplot(clustering_data, aes(mds.x, mds.y, color=group)) + geom_point(size=3) + xlab('Leading logDE dimension 1') + ylab('Leading logDE dimension 2')
##

## Continues with Differential expression preparations:
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
plotBCV(y)
##

## Here we perform replicate correlation analyses using Hmisc:
logcpm <- cpm(y$pseudo.counts, prior.count = 2, log=TRUE) #Creates a matrix of normalized log2 counts per million values from the DGElist object

pcor <- rcorr(logcpm)
melted_pcor <- melt(pcor$r)
ggplot(data=melted_pcor, aes(x=Var1, y=Var2, fill=value))+theme(axis.text.x=element_text(angle=90,hjust=1))+title('Pearson Correlation')+geom_tile()+scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
      midpoint = 0.985, limit = c(0.97,1), space = "Lab", 
      name="Pearson\nCorrelation")
##

## This is for plotting a hierarchical cluster dendrogram of the samples
sampleDists <- dist(t(logcpm))
sampleDistMatrix <- as.matrix(sampleDists)
celltypes <- factor(c("Cell type 1", "Cell type 2", "Cell type 3"))
celltypes <- factor(c('WT Pre', 'WT Post', 'SOD Pre', 'SOD Post'))
rownames(sampleDistMatrix) <- group
colnames(sampleDistMatrix) <- group
annotation <- data.frame(Var1 = factor(1:4, labels = celltypes))
rownames(annotation) <- celltypes
pheatmap(sampleDistMatrix, annotation = annotation)
#df <- as.data.frame(sampleDistMatrix)
#gdf <- grouped_df(df, colnames(sampleDistMatrix))
#colors = colorRampPalette(rev(brewer.pal(9, "Reds")) )(255)
#ha = HeatmapAnnotation(df=data.frame(type=group))
#heatmap.2(sampleDistMatrix, Rowv=sampleDists, Colv=sampleDists, col=colors, legend=3)
##

## This is where we test Differential expression using Fischer's Exact Test
et <- exactTest(y, dispersion="tagwise", pair=c("Cell type 1", "Cell type 2")) # Diserpsion can be changed to "auto", "tagwise", or "trended". The first name in pair is the control
top <- topTags(et)
##

## STOP!!!! The following code is for running GLM linear regression analysis for multiple experimental groups:
data.num <- data.num[1:32600, 1:17] +1 # Add a 1 to all values to remove 0's
y <- DGEList(counts=data.num, group=group.whole)
design <- model.matrix(~0+group.whole, data=y$samples)
colnames(design) <- levels(y$samples$group)
y <- calcNormFactors(y, method="TMM", doWeigthing=TRUE)
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotMDS(y, main="TMM Normalization on whole GLM group")
fit <- glmFit(y, design)
lrt <- glmLRT(fit, contrast=c(0,0,-1,0, 0, 1))
topTags(lrtPre)
##

## Continuing with differential expression analyses and plotting of differentially expressed genes:
summary(de <- decideTestsDGE(et, adjust.method="none", p.value=0.05, lfc=1)) # As written, tells you total number of up or down regulated genes with p-value less than 0.05, lfc greater than 1, regardless of fdr
detags <-rownames(lrt)[as.logical(de)]

tt <- topTags(et, n=32578, adjust.method="fdr")
ttdata<-data.frame(tt)
de <- row.names(ttdata[ttdata$FDR<0.05,]) #This makes a lis(t of the gene names that meet the condition of PValue<0.05. It's meant to bypass the "topTags" feature, which returns FDRs < 0.05
plotSmear(et, de.tags = de, main="Cell Type 2 vs. Cell Type 1", col=ifelse(ttdata[order(ttdata$PValue),]$FDR<0.05, "light coral", "grey55"))
plotSmear(et, main='Cell Type 2 vs. Cell Type 1', col='grey55', pch=16, cex=0.6)
points(subset(et$table$logCPM, rownames(et$table) %in% de), subset(et$table$logFC, rownames(et$table) %in% de), col="lightcoral", pch=16, cex=0.6)
points(subset(et$table$logCPM, rownames(et$table) %in% genesOfInterest), subset(et$table$logFC, rownames(et$table) %in% genesOfInterest), col='green', pch=16, cex=0.8)
abline(0, 0, col="black")

## Ploting a Scatterplot of logcpm values, and labeling specific gene classes:
genesOfInterest  <- c("Gene1", "Gene2", "Gene3")
cellType1.values <- rowMeans(logcpm[, 9:10])
cellType2.values <- rowMeans(logcpm[, 11:12])
plot(cellType1.values, cellType2.values, ylim=c(-4, 15), xlim=c(-4, 15), col="grey55", main="Gene expression changes", pch=16, cex=0.6, ylab= "Expression in Cell Type 2, log10", xlab="Expression in Cell Type 1, log10", legend=TRUE)
points(subset(cellType1.values, rownames(logcpm) %in% de), subset(cellType2.values, rownames(logcpm) %in% de), col="lightcoral", pch=16, cex=0.6)
points(subset(cellType1.values, rownames(logcpm) %in% genesOfInterest), subset(cellType2.values, rownames(logcpm) %in% genesOfInterest), col='black', pch=16, cex=0.8)
abline(0, 1, col='blue')
## Save the results:
write.table(topTags(et, n=32578), file="DE_file.txt", sep="\t", row.names=TRUE, col.names=TRUE)
