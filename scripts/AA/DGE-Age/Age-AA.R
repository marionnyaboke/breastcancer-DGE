library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(HTSFilter)
library(ggplot2) #Best plots
library(ggrepel)
library(dplyr)

# set the directory from which files are imported
directory <- "C:/Users/Marion/Desktop/DGE/AA/"
dir(directory)

# read files
rawCountTable <- read.csv(paste0(directory,"reads_count_AA.csv"), header=TRUE,
                          row.names=1)
sampleInfo <- read.csv(paste0(directory,"Metadata-final-AA.csv"), header=TRUE,
                       row.names=1)

head(rawCountTable)
nrow(rawCountTable)

rawCountTable <- rawCountTable[,match(rownames(sampleInfo),
                                      colnames(rawCountTable))]
colnames(rawCountTable)
rownames(sampleInfo)

# create the ‘condition’ column
condition = sampleInfo$cases.0.samples.0.sample_type

sampleInfo = data.frame(sampleInfo, condition)
#sampleInfo = subset(sampleInfo, select= -condition)

# create a DGEList data object
dgeFull <- DGEList(rawCountTable, group=sampleInfo$condition)


# add the sample information object in the DGEList data
dgeFull$sampleInfo <- sampleInfo

# check the number of genes with no expression in all samples
table(rowSums(dgeFull$counts==0)==15)

#FALSE  TRUE 
#59246   431  

# Filtering non-expressed and lowly-expressed genes.
keep.exprs <- filterByExpr(dgeFull, group=sampleInfo$condition)
filtered.counts <- dgeFull[keep.exprs,, keep.lib.sizes=FALSE]

# preparing the data object for the analysis 
# select the subset paired-end samples from dgeFull
filtered.counts <- DGEList(filtered.counts$counts)
filtered.counts$sampleInfo <- dgeFull$sampleInfo[dgeFull$sampleInfo$type=="paired-end",]

table(rowSums(filtered.counts$counts==0)==15)
#FALSE  TRUE 
#27168   137

# Apply sample grouping based on age range from which the sample was derived
design <- model.matrix(~0+sampleInfo$age_range+sampleInfo$condition)
#colnames(design) <- levels(samples.metadata.clean$Tissue)
colnames(design) <- c("Mid", "Old", "Young", "sample_type")

# Estimate dispersions for tags
filtered.counts.dge <- estimateDisp(filtered.counts, design, robust = TRUE)

# Fit a generalized likelihood model to the DGELIST using sample grouping
fit <- glmFit(filtered.counts.dge,design)

#################################################################
# code in this section adapted from https://github.com/iscb-dc-rsg/2016-summer-workshop
# generate a list of all possible pairwise contrasts
sampleInfo$age_range <- as.factor(sampleInfo$age_range)
condition_pairs <- t(combn(levels(sampleInfo$age_range), 2))

comparisons <- list()
for (i in 1:nrow(condition_pairs)) {
  comparisons[[i]] <- as.character(condition_pairs[i,])
}

# remove MG to SG comparison
#comparisons[[2]] <- NULL

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
#1566   47

################################################################
# Obtain the counts of genes expressed for each contrast individually
# This aims to obtain the number of genes differentially expressed between 
# the 3 stages of development i.e. MG -> PV, PV -> SG
# Likelihood ratio test to identify DEGs
# Under50 compared to 50- 75
YoungvsMid <- makeContrasts(Young-Mid, levels=design)
Young_vs_Mid_lrt <- glmLRT(fit, contrast=YoungvsMid)

# Mid compared to Old
MidvsOld <- makeContrasts(Mid-Old, levels=design)
Mid_vs_Old_lrt <- glmLRT(fit, contrast=MidvsOld)

# Old compared to Young
OldvsYoung <- makeContrasts(Old-Young, levels=design)
Old_vs_Young_lrt <- glmLRT(fit, contrast=OldvsYoung)

# Genes with most significant differences (using topTags)
# Young compared to Mid
topGenes_Young <- topTags(Young_vs_Mid_lrt, adjust.method = "fdr", p.value = 0.05, n=Inf)
dim(topGenes_Young)
#1483    5

# Mid compared to Old
topGenes_Mid <- topTags(Mid_vs_Old_lrt, adjust.method = "fdr", p.value = 0.05, n=Inf)
dim(topGenes_Mid)
#60  5

# Old compared to Young
topGenes_Old <- topTags(Old_vs_Young_lrt, adjust.method = "fdr", p.value = 0.05, n=Inf)
dim(topGenes_Old)
#69  5

#Total number of genes: 1612

#######################################################################################
# DE genes at 5% FDR (using decideTestsDGE function)
#
# Young compared to Mid
Young_vs_Mid_de.genes <- decideTestsDGE(Young_vs_Mid_lrt, adjust.method = "fdr", p.value = 0.05)
# get summary
summary(Young_vs_Mid_de.genes)
#       -1*Mid 1*Young
#Down              427
#NotSig          25822
#Up               1056

# Mid compared to Old
Mid_vs_Old_de.genes <- decideTestsDGE(Mid_vs_Old_lrt, adjust.method = "fdr", p.value = 0.05)
# get summary
summary(Mid_vs_Old_de.genes)
#       1*Mid -1*Old
#Down             60
#NotSig        27245
#Up                0

# Old compared to Young
Old_vs_Young_de.genes <- decideTestsDGE(Old_vs_Young_lrt, adjust.method = "fdr", p.value = 0.05)
# get summary
summary(Old_vs_Young_de.genes)
#       1*Old -1*Young
#Down                2
#NotSig          27236
#Up                 67

# DE genes in the Young that are common in all comparisons
de.common <- which(Young_vs_Mid_de.genes[,1]!=0 & Mid_vs_Old_de.genes[,1]!=0 & Old_vs_Young_de.genes[,1]!=0)
length(de.common)
#1

de.common.df <- as.data.frame(de.common)
de.common.df <- tibble::rownames_to_column(de.common.df, var = "gene_id")

# create a dataframe with data on PV and SG differential gene expression
Young_data <- topGenes_Young$table
Mid_data <- topGenes_Mid$table
Old_data <- topGenes_Old$table

Young_data <- tibble::rownames_to_column(Young_data, var = "gene_id")
Mid_data <- tibble::rownames_to_column(Mid_data, var = "gene_id")
Old_data <- tibble::rownames_to_column(Old_data, var = "gene_id")

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

Young_vs_Mid <- Young_data[order(Young_data$logFC, decreasing = TRUE),]
write.csv(Young_vs_Mid, file ="C:/Users/Marion/Desktop/DGE/AA/Age-AA/YvsM-DGE.csv", row.names = FALSE)

Mid_vs_Old <- Mid_data[order(Young_data$logFC, decreasing = TRUE),]
write.csv(Mid_vs_Old, file ="C:/Users/Marion/Desktop/DGE/AA/Age-AA/MvsO-DGE.csv", row.names = FALSE)

Old_vs_Young <- Old_data[order(Young_data$logFC, decreasing = TRUE),]
write.csv(Old_vs_Young, file ="C:/Users/Marion/Desktop/DGE/AA/Age-AA/OvsY-DGE.csv", row.names = FALSE)

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
significance = -log10(topGenes_Young$table$FDR) < 1
topGenes_Young$table$significance <- significance
dim(topGenes_Young$table)

## setting the values
topGenes_Young$table$diffexpressed <- "No"

# if logFC > 0.1 and FDR < 0.05, set as "UP"
topGenes_Young$table$diffexpressed[topGenes_Young$table$logFC > 1 & topGenes_Young$table$FDR < 0.05] <- "Up"

# if logFC < -1 and pvalue < 0.05, set as "DOWN"
topGenes_Young$table$diffexpressed[topGenes_Young$table$logFC < -1 & topGenes_Young$table$FDR < 0.05] <- "Down"

# set different colors
mycolors <- c("cornflowerblue", "black", "firebrick")
names(mycolors) <- c("Down", "No", "Up")

dev.off()


input <- cbind(gene=rownames(topGenes_Young$table), topGenes_Young$table) #convert the rownames to a column
volc = ggplot(input, aes(logFC, -log10(PValue))) + #volcanoplot with logFC versus pvalue
  geom_point() +
  theme_classic()+
  geom_vline(xintercept=c(-1.0, 0.8), col="orange", linetype="dashed") +
  geom_hline(yintercept=-log10(0.003), col="orange" , linetype="dashed")+
  geom_point(aes(col=diffexpressed)) + #add points colored by significance
  scale_color_manual(values= mycolors) + 
  ggtitle("25-50yo AA patients vs 50-75yo AA patients")
volc <- volc+geom_text_repel(data=head(input, 1000), aes(label=gene), size = 3) #adding text for the top 20 genes
volc
#traceback()

# create a Volcano plot- Mid

# create a column with thresholds
significance = -log10(topGenes_Mid$table$FDR) < 1
topGenes_Mid$table$significance <- significance
dim(topGenes_Mid$table)

## setting the values
topGenes_Mid$table$diffexpressed <- "No"

# if logFC > 0.1 and FDR < 0.05, set as "UP"
topGenes_Mid$table$diffexpressed[topGenes_Mid$table$logFC > 1 & topGenes_Mid$table$FDR < 0.05] <- "Up"

# if logFC < -1 and pvalue < 0.05, set as "DOWN"
topGenes_Mid$table$diffexpressed[topGenes_Mid$table$logFC < -1 & topGenes_Mid$table$FDR < 0.05] <- "Down"

# set different colors
mycolors <- c("cornflowerblue", "black", "firebrick")
names(mycolors) <- c("Down", "No", "Up")

dev.off()


input <- cbind(gene=rownames(topGenes_Mid$table), topGenes_Mid$table) #convert the rownames to a column
volc = ggplot(input, aes(logFC, -log10(PValue))) + #volcanoplot with logFC versus pvalue
  geom_point() +
  theme_classic()+
  geom_vline(xintercept=c(-1.0, 0.8), col="orange", linetype="dashed") +
  geom_hline(yintercept=-log10(0.003), col="orange" , linetype="dashed")+
  geom_point(aes(col=diffexpressed)) + #add points colored by significance
  scale_color_manual(values= mycolors) + 
  ggtitle("50-75yo AA patients vs 75+yo AA patients")
volc <- volc+geom_text_repel(data=head(input, 1000), aes(label=gene), size = 3) #adding text for the top 20 genes
volc
#traceback()


# create a Volcano plot- Old

# create a column with thresholds
significance = -log10(topGenes_Old$table$FDR) < 1
topGenes_Old$table$significance <- significance
dim(topGenes_Old$table)

## setting the values
topGenes_Old$table$diffexpressed <- "No"

# if logFC > 0.1 and FDR < 0.05, set as "UP"
topGenes_Old$table$diffexpressed[topGenes_Old$table$logFC > 1 & topGenes_Old$table$FDR < 0.05] <- "Up"

# if logFC < -1 and pvalue < 0.05, set as "DOWN"
topGenes_Old$table$diffexpressed[topGenes_Old$table$logFC < -1 & topGenes_Old$table$FDR < 0.05] <- "Down"

# set different colors
mycolors <- c("cornflowerblue", "black", "firebrick")
names(mycolors) <- c("Down", "No", "Up")

dev.off()


input <- cbind(gene=rownames(topGenes_Old$table), topGenes_Old$table) #convert the rownames to a column
volc = ggplot(input, aes(logFC, -log10(PValue))) + #volcanoplot with logFC versus pvalue
  geom_point() +
  theme_classic()+
  geom_vline(xintercept=c(-1.0, 0.8), col="orange", linetype="dashed") +
  geom_hline(yintercept=-log10(0.003), col="orange" , linetype="dashed")+
  geom_point(aes(col=diffexpressed)) + #add points colored by significance
  scale_color_manual(values= mycolors) + 
  ggtitle("75+yo AA patients vs 25-50yo AA patients")
volc <- volc+geom_text_repel(data=head(input, 1000), aes(label=gene), size = 3) #adding text for the top 20 genes
volc
#traceback()

# transform the normalized counts in log-counts-per-million
y <- cpm(filtered.counts, log=TRUE, prior.count = 1)

# select 10% differentially expressed genes and produce a heatmap
selY <- y[rownames(topGenes_Old$table)[topGenes_Old$table$FDR<0.01 & 
                                             abs(topGenes_Old$table$logFC)>1.5],]

cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]

dev.off() 

cim(selY, margins = c(10, 10))
