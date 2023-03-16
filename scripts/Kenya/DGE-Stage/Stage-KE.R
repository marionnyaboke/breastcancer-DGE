library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(HTSFilter)
library(ggplot2) #Best plots
library(ggrepel)
library(dplyr)

# set the directory from which files are imported
directory <- "C:/Users/Marion/Desktop/DGE/DGE-Kenya/"
dir(directory)

# read files
rawCountTable <- read.csv(paste0(directory,"reads_count.csv"), header=TRUE,
                          row.names=1)
sampleInfo <- read.csv(paste0(directory,"metadata.csv"), header=TRUE,
                       row.names=1)

head(rawCountTable)
nrow(rawCountTable)

rawCountTable <- rawCountTable[,match(rownames(sampleInfo),
                                      colnames(rawCountTable))]
colnames(rawCountTable)
rownames(sampleInfo)

# create the ‘condition’ column
condition = sampleInfo$sample_type

sampleInfo = data.frame(sampleInfo, condition)
#sampleInfo = subset(sampleInfo, select= -condition)

# create a DGEList data object
dgeFull <- DGEList(rawCountTable, group=sampleInfo$condition)


# add the sample information object in the DGEList data
dgeFull$sampleInfo <- sampleInfo

# check the number of genes with no expression in all samples
table(rowSums(dgeFull$counts==0)==15)

#FALSE  TRUE 
#39218   371 

# Filtering non-expressed and lowly-expressed genes.
keep.exprs <- filterByExpr(dgeFull, group=sampleInfo$condition)
filtered.counts <- dgeFull[keep.exprs,, keep.lib.sizes=FALSE]

# preparing the data object for the analysis 
# select the subset paired-end samples from dgeFull
filtered.counts <- DGEList(filtered.counts$counts)
filtered.counts$sampleInfo <- dgeFull$sampleInfo[dgeFull$sampleInfo$type=="paired-end",]

table(rowSums(filtered.counts$counts==0)==15)
#FALSE  TRUE 
#19882    92 

# Apply sample grouping based on age range from which the sample was derived
design <- model.matrix(~0+sampleInfo$ajcc_pathologic_stage+sampleInfo$condition)
#colnames(design) <- levels(samples.metadata.clean$Tissue)
colnames(design) <- c("Missing", "StageI", "StageII", "StageIII", "sample_type")
ncol(y)
nrow(design)

# Estimate dispersions for tags
filtered.counts.dge <- estimateDisp(filtered.counts, design, robust = TRUE)

ncol(filtered.counts)
# Fit a generalized likelihood model to the DGELIST using sample grouping
fit <- glmFit(filtered.counts.dge,design)

#################################################################
# code in this section adapted from https://github.com/iscb-dc-rsg/2016-summer-workshop
# generate a list of all possible pairwise contrasts
sampleInfo$ajcc_pathologic_stage <- as.factor(sampleInfo$ajcc_pathologic_stage)
condition_pairs <- t(combn(levels(sampleInfo$ajcc_pathologic_stage), 2))

comparisons <- list()
for (i in 1:nrow(condition_pairs)) {
  comparisons[[i]] <- as.character(condition_pairs[i,])
}

# remove MG to SG comparison
comparisons[[4]] <- NULL

# vector to store deferentially expressed genes
sig_genes <- c()

# iterate over the contrasts, and perform a differential expression test for each pair
for (conds in comparisons) {
  # generate string contrast formula
  contrast_formula <- paste(conds, collapse=' - ')
  contrast_mat <- makeContrasts(contrasts=contrast_formula, levels=design)
  contrast_lrt <- glmLRT(fit, contrast=contrast_mat)
  topGenes <- topTags(contrast_lrt, n=Inf, p.value=0.05, adjust.method = "fdr")
  
  # Grab highly ranked genes
  sig_genes <- union(sig_genes, rownames(topGenes$table))
}

# Filter out genes which were not differentially expressed for any contrast
de.genes <- filtered.counts.dge[rownames(filtered.counts.dge) %in% sig_genes,]
dim(de.genes$counts)
#576  47

################################################################
# Obtain the counts of genes expressed for each contrast individually
# This aims to obtain the number of genes differentially expressed between 
# the 3 stages of development i.e. MG -> PV, PV -> SG
# Likelihood ratio test to identify DEGs

my.contrasts <- makeContrasts(IvsII=StageI-StageII, IvsIII=StageI-StageIII, 
                                IIvsIII=StageII-StageIII, levels=design)

IvsII_lrt <- glmLRT(fit, contrast=my.contrasts[,"IvsII"])

IvsIII_lrt <- glmLRT(fit, contrast=my.contrasts[,"IvsIII"])

IIvsIII_lrt <- glmLRT(fit, contrast=my.contrasts[,"IIvsIII"])


# Genes with most significant differences (using topTags)
# IvsII
topGenes_IvsII <- topTags(IvsII_lrt, adjust.method = "fdr", p.value = 0.05, n=Inf)
dim(topGenes_IvsII)
#5214    5

# IvsIII
topGenes_IvsIII <- topTags(IvsIII_lrt, adjust.method = "fdr", p.value = 0.05, n=Inf)
dim(topGenes_IvsIII)
#6838    5

# IIvsIII
topGenes_IIvsIII <- topTags(IIvsIII_lrt, adjust.method = "fdr", p.value = 0.05, n=Inf)
dim(topGenes_IIvsIII)
#236   5


#Total number of genes: 776

#######################################################################################
# DE genes at 5% FDR (using decideTestsDGE function)
#
# IvsII
IvsII_de.genes <- decideTestsDGE(IvsII_lrt, adjust.method = "fdr", p.value = 0.05)
# get summary
summary(IvsII_de.genes)
#       1*StageI -1*StageII
#Down                  4480
#NotSig               14760
#Up                     734

# IvsIII
IvsIII_de.genes <- decideTestsDGE(IvsIII_lrt, adjust.method = "fdr", p.value = 0.05)
# get summary
summary(IvsIII_de.genes)
#       1*StageI -1*StageIII
#Down                   6292
#NotSig                13136
#Up                      546

# IIvsIII
IIvsIII_de.genes <- decideTestsDGE(IIvsIII_lrt, adjust.method = "fdr", p.value = 0.05)
# get summary
summary(IIvsIII_de.genes)
#       1*StageII -1*StageIII
#Down                     230
#NotSig                 19738
#Up                         6

# DE genes in the Young that are common in all comparisons
de.common <- which(IIIvsIV_de.genes[,1]!=0 & IIvsIV_de.genes[,1]!=0 & IIvsIII_de.genes[,1]!=0)
length(de.common)
#1

de.common.df <- as.data.frame(de.common)
de.common.df <- tibble::rownames_to_column(de.common.df, var = "gene_id")

# create a dataframe with data on PV and SG differential gene expression
IvsII_data <- topGenes_IvsII$table
IvsIII_data <- topGenes_IvsIII$table
IIvsIII_data <- topGenes_IIvsIII$table

IvsII_data <- tibble::rownames_to_column(IvsII_data, var = "gene_id")
IvsIII_data <- tibble::rownames_to_column(IvsIII_data, var = "gene_id")
IIvsIII_data <- tibble::rownames_to_column(IIvsIII_data, var = "gene_id")

# obtain the common genes for each comparison
Young_data_common_de_genes <- Young_data %>% filter(gene_id %in% de.common.df$gene_id)
Mid_data_common_de_genes <- Mid_data %>% filter(gene_id %in% de.common.df$gene_id)
Old_data_common_de_genes <- Old_data %>% filter(gene_id %in% de.common.df$gene_id)

#PV_data_common_de_genes <- PV_data_common_de_genes[order(PV_data_common_de_genes$logFC, 
#                                                         decreasing = TRUE),]
#SG_data_common_de_genes <- SG_data_common_de_genes[order(SG_data_common_de_genes$logFC, 
#                                                         decreasing = TRUE),]

# write out to excel upregulated and downregulated genes and the commmon genes between contrasts
#PV_vs_MG <- PV_data[order(PV_data$logFC, decreasing = TRUE),]

IvsII <- IvsII_data[order(IvsII_data$logFC, decreasing = TRUE),]
write.csv(IvsII, file ="C:/Users/Marion/Desktop/DGE/DGE-Kenya/Stage-KE/IvsII-DGE.csv", row.names = FALSE)

IvsIII <- IvsIII_data[order(IvsIII_data$logFC, decreasing = TRUE),]
write.csv(IvsIII, file ="C:/Users/Marion/Desktop/DGE/DGE-Kenya/Stage-KE/IvsIII-DGE.csv", row.names = FALSE)

IIvsIII <- IIvsIII_data[order(IIvsIII_data$logFC, decreasing = TRUE),]
write.csv(IIvsIII, file ="C:/Users/Marion/Desktop/DGE/DGE-Kenya/Stage-KE/IIvsIII-DGE.csv", row.names = FALSE)

# write.xlsx(PV_vs_MG, file = "../results/differentially_expressed_genes.xlsx",
#            sheetName = "MG vs PV", row.names = FALSE)
# write.xlsx(SG_vs_PV, file = "../results/differentially_expressed_genes.xlsx",
#            sheetName = "PV vs SG", append = TRUE, row.names = FALSE)
# 
# write.xlsx(PV_data_common_de_genes, 
#            file = "../results/differentially_expressed_genes_common_in_contrasts.xlsx",
#            sheetName = "MG vs PV Common genes", append = TRUE, row.names = FALSE)
# write.xlsx(SG_data_common_de_genes, 
#            file = "../results/differentially_expressed_genes_common_in_contrasts.xlsx",
#            sheetName = "PV vs SG Common genes", append = TRUE, row.names = FALSE)
# assign "select" function to dplyr in the evironment
#assign('select', dplyr::select, envir=.GlobalEnv)
#PV_contrast_common_de_genes_logFC <- PV_data_common_de_genes %>% select(gene_id, logFC_MGvsPV=logFC)
#SG_contrast_common_de_genes_logFC <- SG_data_common_de_genes %>% select(gene_id, logFC_PVvsSG=logFC)
#contrast_common_de_genes_logFC <- merge(PV_contrast_common_de_genes_logFC, 
#                                        SG_contrast_common_de_genes_logFC, by="gene_id")
# write.xlsx(contrast_common_de_genes_logFC, 
#            file = "../results/differentially_expressed_genes_common_in_contrasts_logFC_comparison.xlsx",
#            sheetName = "Common DEGs LogFC Comparisons", row.names = FALSE)
#```


#```{r, eval=FALSE, echo=FALSE}

# Plotting to visually inspect differential gene expression results.
# Differential expression analysis - plots
#
# Interpreting the DE analysis results

# create a Volcano plot- Young

# create a column with thresholds
significance = -log10(topGenes_IIvsIII$table$FDR) < 1
topGenes_IIvsIII$table$significance <- significance
dim(topGenes_IIvsIII$table)

## setting the values
topGenes_IIvsIII$table$diffexpressed <- "No"

# if logFC > 0.1 and FDR < 0.05, set as "UP"
topGenes_IIvsIII$table$diffexpressed[topGenes_IIvsIII$table$logFC > 1 & topGenes_IIvsIII$table$FDR < 0.05] <- "Up"

# if logFC < -1 and pvalue < 0.05, set as "DOWN"
topGenes_IIvsIII$table$diffexpressed[topGenes_IIvsIII$table$logFC < -1 & topGenes_IIvsIII$table$FDR < 0.05] <- "Down"

# set different colors
mycolors <- c("cornflowerblue", "black", "firebrick")
names(mycolors) <- c("Down", "No", "Up")

dev.off()

IvsII 
IvsIII
IvsIV
IIvsIII
IIvsIV
IIIvsIV

par(mfrow=c(3,3))
plot(volc)
graphics.off()
dev.list()

input <- cbind(gene=rownames(topGenes_IIvsIII$table), topGenes_IIvsIII$table) #convert the rownames to a column
volc = ggplot(input, aes(logFC, -log10(PValue))) + #volcanoplot with logFC versus pvalue
  geom_point() +
  theme_classic()+
  geom_vline(xintercept=c(-1.0, 0.8), col="orange", linetype="dashed") +
  geom_hline(yintercept=-log10(0.003), col="orange" , linetype="dashed")+
  geom_point(aes(col=diffexpressed)) + #add points colored by significance
  scale_color_manual(values= mycolors) + 
  ggtitle("Stage II vs Stage III KE patients")
volc+geom_text_repel(data=head(input, 1000), aes(label=gene), size = 3) #adding text for the top 20 genes

#traceback()


# transform the normalized counts in log-counts-per-million
y <- cpm(filtered.counts, log=TRUE, prior.count = 1)

# select 10% differentially expressed genes and produce a heatmap
selY <- y[rownames(topGenes_IvsIII$table)[topGenes_IvsIII$table$FDR<0.01 & 
                                             abs(topGenes_IvsIII$table$logFC)>1.5],]

cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]

dev.off() 

cim(selY, margins = c(10, 10))
