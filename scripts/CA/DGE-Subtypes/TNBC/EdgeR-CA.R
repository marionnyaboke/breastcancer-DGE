library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(HTSFilter)
library(ggplot2) #Best plots
library(ggrepel)
library(dplyr)

# set the directory from which files are imported
directory <- "C:/Users/Marion/Desktop/DGE/CA/all-counts/TNBC/"
dir(directory)

# read files
rawCountTable <- read.csv(paste0(directory,"reads_count.csv"), header=TRUE,
                            row.names=1)
sampleInfo <- read.csv(paste0(directory,"metadata-TNBC.csv"), header=TRUE,
                         row.names=1)

head(rawCountTable)
nrow(rawCountTable)

rawCountTable <- rawCountTable[,match(rownames(sampleInfo),
                                      colnames(rawCountTable))]
colnames(rawCountTable)
rownames(sampleInfo)

# create the ‘condition’ column
condition = c("Primary Tumor", "Primary Tumor", "Solid Tissue Normal", "Solid Tissue Normal", 
              "Primary Tumor", "Primary Tumor", "Primary Tumor", "Primary Tumor", "Primary Tumor", 
              "Solid Tissue Normal", "Solid Tissue Normal", "Solid Tissue Normal", "Solid Tissue Normal",
              "Solid Tissue Normal", "Solid Tissue Normal","Solid Tissue Normal", "Solid Tissue Normal",
              "Solid Tissue Normal", "Primary Tumor", "Solid Tissue Normal", "Primary Tumor", "Primary Tumor",
              "Primary Tumor", "Primary Tumor")

sampleInfo = data.frame(sampleInfo, condition)
#sampleInfo = subset(sampleInfo, select= -condition)

# create a DGEList data object
dgeFull <- DGEList(rawCountTable, group=sampleInfo$condition)


# add the sample information object in the DGEList data
dgeFull$sampleInfo <- sampleInfo

# check the number of genes with no expression in all samples
table(rowSums(dgeFull$counts==0)==15)

#FALSE  TRUE 
#59618   995  

# Filtering non-expressed and lowly-expressed genes.
keep.exprs <- filterByExpr(dgeFull, group=sampleInfo$condition)
filtered.counts <- dgeFull[keep.exprs,, keep.lib.sizes=FALSE]

# preparing the data object for the analysis 
# select the subset paired-end samples from dgeFull
dge <- DGEList(filtered.counts$counts)
dge$sampleInfo <- dgeFull$sampleInfo[dgeFull$sampleInfo$type=="paired-end",]

table(rowSums(filtered.counts$counts==0)==15)
#FALSE 
#20677  

# preparing the data object for the analysis 
# select the subset paired-end samples from degFull
#dge <- DGEList(dgeFull$counts)
#dge$sampleInfo <- dgeFull$sampleInfo[dgeFull$sampleInfo$type=="paired-end",]

# data exploration and quality assessment
# extract pseudo-counts (ie \(\log_2(K+1)\))
pseudoCounts <- log2(dge$counts+1)
head(pseudoCounts)

# histogram for pseudo-counts
hist(pseudoCounts[,"SRR10729843"])

# boxplot for pseudo-counts
boxplot(pseudoCounts, col="cyan")

# MDS for pseudo-counts (using limma package)
plotMDS(pseudoCounts)

# PCA plot
plotMDS(pseudoCounts, gene.selection="common")

# Differential expression analysis
# remove genes with zero counts for all samples

dge <- DGEList(dge$counts[apply(dge$counts, 1, sum) != 0, ],
               group=sampleInfo$condition)
dge$sampleInfo <- dge$sampleInfo

# estimate the normalization factors
dge <- calcNormFactors(dge, method="TMM")
dge$samples$group <- relevel(dge$samples$group, ref="Solid Tissue Normal") #sets non-african as the control group

# estimate common and tagwise dispersion
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

# perform an exact test for the difference in expression between the conditions “treated” and “control”
dgeTest <- exactTest(dge)

# Independant filtering
#  remove low expressed genes

filtData <- HTSFilter(dge)$filteredData

dgeTestFilt <- exactTest(filtData)

# Diagnostic plot for multiple testing
# plot a histogram of unadjusted p-values
hist(dgeTest$table[,"PValue"], breaks=50)

# plot a histogram of unadjusted p-values after filtering
hist(dgeTestFilt$table[,"PValue"], breaks=50)

# Inspecting the results
# extract a summary of the differential expression statistics
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))
head(resNoFilt)

resFilt <- topTags(dgeTestFilt, n=nrow(dgeTest$table))
head(resFilt)

#resFilt$genelabels <- factor(resFilt$gene, levels = c("H2BC5", "CCNF", "H2BC21", "H3C4"))
#resFilt$
# compare the number of differentially expressed genes with and without filtering
# before independent filtering
sum(resNoFilt$table$FDR < 0.05)

# after independent filtering
sum(resFilt$table$FDR < 0.05)

# extract and sort differentially expressed genes
sigDownReg <- resFilt$table[resFilt$table$FDR<0.05,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
head(sigDownReg)

sigUpReg <- sigDownReg[order(sigDownReg$logFC, decreasing=TRUE),]
head(sigUpReg)

# write the results in csv files
write.csv(sigDownReg, file="C:/Users/Marion/Desktop/DGE/CA/all-counts/TNBC/sigDownReg-filtered.csv")
write.csv(sigUpReg, file="C:/Users/Marion/Desktop/DGE/CA/all-counts/TNBC/sigUpReg-filtered.csv")

# Interpreting the DE analysis results
# create a MA plot with 5% differentially expressed genes

plotSmear(dgeTestFilt,
          de.tags = rownames(resFilt$table)[which(resFilt$table$FDR<0.05)])

# create a Volcano plot

# create a column with thresholds
significance = -log10(resFilt$table$FDR) < 1
resFilt$table$significance <- significance
dim(resFilt$table)

## setting the values
resFilt$table$diffexpressed <- "No"

# if logFC > 0.1 and FDR < 0.05, set as "UP"
resFilt$table$diffexpressed[resFilt$table$logFC > 1 & resFilt$table$FDR < 0.05] <- "Up"

# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
resFilt$table$diffexpressed[resFilt$table$logFC < -1 & resFilt$table$FDR < 0.05] <- "Down"

# set different colors
mycolors <- c("cornflowerblue", "black", "firebrick")
names(mycolors) <- c("Down", "No", "Up")

dev.off()


input <- cbind(gene=rownames(resFilt$table), resFilt$table) #convert the rownames to a column
volc = ggplot(input, aes(logFC, -log10(PValue))) + #volcanoplot with logFC versus pvalue
  geom_point() +
  theme_classic()+
  geom_vline(xintercept=c(-1.0, 0.8), col="orange", linetype="dashed") +
  geom_hline(yintercept=-log10(0.01), col="orange" , linetype="dashed")+
  geom_point(aes(col=diffexpressed)) + #add points colored by significance
  scale_color_manual(values= mycolors) + 
  ggtitle("Tumor vs Normal")
volc <- volc+geom_text_repel(data=head(input, 500), aes(label=gene), size = 3) #adding text for the top 20 genes
volc
#traceback()

# transform the normalized counts in log-counts-per-million
y <- cpm(dge, log=TRUE, prior.count = 1)

# select 1% differentially expressed genes and produce a heatmap
selY <- y[rownames(resFilt$table)[resFilt$table$FDR<0.01 & 
                                    abs(resFilt$table$logFC)>1.5],]

cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]

dev.off() 

cim(selY, margins = c(10, 10))

