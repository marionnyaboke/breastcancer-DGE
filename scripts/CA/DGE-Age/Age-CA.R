library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(HTSFilter)
library(ggplot2) #Best plots
library(ggrepel)
library(dplyr)

# set the directory from which files are imported
directory <- "C:/Users/Marion/Desktop/DGE/CA/"
dir(directory)

# read files
rawCountTable <- read.csv(paste0(directory,"reads_count_CA.csv"), header=TRUE,
                          row.names=1)
sampleInfo <- read.csv(paste0(directory,"Metadata-final-CA.csv"), header=TRUE,
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
#60270   343  

# Filtering non-expressed and lowly-expressed genes.
keep.exprs <- filterByExpr(dgeFull, group=sampleInfo$condition)
filtered.counts <- dgeFull[keep.exprs,, keep.lib.sizes=FALSE]

# preparing the data object for the analysis 
# select the subset paired-end samples from dgeFull
filtered.counts <- DGEList(filtered.counts$counts)
filtered.counts$sampleInfo <- dgeFull$sampleInfo[dgeFull$sampleInfo$type=="paired-end",]

table(rowSums(filtered.counts$counts==0)==15)
#FALSE  TRUE 
#22610    39  

# Apply sample grouping based on age range from which the sample was derived
design <- model.matrix(~0+sampleInfo$age_range+sampleInfo$condition)
#colnames(design) <- levels(samples.metadata.clean$Tissue)
colnames(design) <- c("Young","Mid", "sample_type")
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
#118  66

################################################################
# Obtain the counts of genes expressed for each contrast individually
# This aims to obtain the number of genes differentially expressed between 
# the 3 stages of development i.e. MG -> PV, PV -> SG
# Likelihood ratio test to identify DEGs
# Under50 compared to 50- 75
under50vsMid <- makeContrasts(Young-Mid, levels=design)
Under50_vs_Mid_lrt <- glmLRT(fit, contrast=under50vsMid)

# PV compared to MG
#PV_vs_MG_lrt <- glmLRT(fit, contrast = c(-1,1,0,0))

# Genes with most significant differences (using topTags)
# Young compared to Mid
topGenes_Under50 <- topTags(Under50_vs_Mid_lrt, adjust.method = "fdr", p.value = 0.05, n=Inf)
dim(topGenes_Under50)
#118   5

# PV compared to MG
#topGenes_PV <- topTags(PV_vs_MG_lrt, adjust.method = "BH", p.value = 0.05, n=Inf)
#dim(topGenes_PV)
#1908    5
#Total number of genes: 5074

#######################################################################################
# DE genes at 5% FDR (using decideTestsDGE function)
#
# Under50 compared to 50- 75
Under50_vs_Mid_de.genes <- decideTestsDGE(Under50_vs_Mid_lrt, adjust.method = "fdr", p.value = 0.05)
# get summary
summary(Under50_vs_Mid_de.genes)
#       -1*Under50 1*Mid
#Down                 94
#NotSig            22531
#Up                   24

# PV compared to MG
#PV_vs_MG_de.genes <- decideTestsDGE(PV_vs_MG_lrt, adjust.method = "BH", p.value = 0.05)
# summary
#summary(PV_vs_MG_de.genes)
#       -1*MG 1*PV
#Down         987
#NotSig       5482
#Up           921

# DE genes in the PV that are common in both comparisons
#de.common <- which(PV_vs_MG_de.genes[,1]!=0 & SG_vs_PV_de.genes[,1]!=0)
#length(de.common)
#1140 
#de.common.df <- as.data.frame(de.common)
#de.common.df <- tibble::rownames_to_column(de.common.df, var = "gene_id")

# create a dataframe with data on PV and SG differential gene expression
Under50_data <- topGenes_Under50$table
#SG_data <- topGenes_SG$table
Under50_data <- tibble::rownames_to_column(Under50_data, var = "gene_id")
#SG_data <- tibble::rownames_to_column(SG_data, var = "gene_id")
# obtain the common genes for each comparison
#PV_data_common_de_genes <- PV_data %>% filter(gene_id %in% de.common.df$gene_id)
#SG_data_common_de_genes <- SG_data %>% filter(gene_id %in% de.common.df$gene_id)
#PV_data_common_de_genes <- PV_data_common_de_genes[order(PV_data_common_de_genes$logFC, 
#                                                         decreasing = TRUE),]
#SG_data_common_de_genes <- SG_data_common_de_genes[order(SG_data_common_de_genes$logFC, 
#                                                         decreasing = TRUE),]
# write out to excel upregulated and downregulated genes and the commmon genes between contrasts
#PV_vs_MG <- PV_data[order(PV_data$logFC, decreasing = TRUE),]
Under50_vs_Mid <- Under50_data[order(Under50_data$logFC, decreasing = TRUE),]
write.csv(Under50_vs_Mid, file ="C:/Users/Marion/Desktop/DGE/CA/Age/differentially_expressed_genes.csv", row.names = FALSE)

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
                                        SG_contrast_common_de_genes_logFC, by="gene_id")
# write.xlsx(contrast_common_de_genes_logFC, 
#            file = "../results/differentially_expressed_genes_common_in_contrasts_logFC_comparison.xlsx",
#            sheetName = "Common DEGs LogFC Comparisons", row.names = FALSE)
```


```{r, eval=FALSE, echo=FALSE}
# Plotting to visually inspect differential gene expression results.
# Differential expression analysis - plots
#
# Interpreting the DE analysis results

# create a Volcano plot

# create a column with thresholds
significance = -log10(topGenes_Under50$table$FDR) < 1
topGenes_Under50$table$significance <- significance
dim(topGenes_Under50$table)

## setting the values
topGenes_Under50$table$diffexpressed <- "No"

# if logFC > 0.1 and FDR < 0.05, set as "UP"
topGenes_Under50$table$diffexpressed[topGenes_Under50$table$logFC > 1 & topGenes_Under50$table$FDR < 0.05] <- "Up"

# if logFC < -1 and pvalue < 0.05, set as "DOWN"
topGenes_Under50$table$diffexpressed[topGenes_Under50$table$logFC < -1 & topGenes_Under50$table$FDR < 0.05] <- "Down"

# set different colors
mycolors <- c("cornflowerblue", "black", "firebrick")
names(mycolors) <- c("Down", "No", "Up")

dev.off()


input <- cbind(gene=rownames(topGenes_Under50$table), topGenes_Under50$table) #convert the rownames to a column
volc = ggplot(input, aes(logFC, -log10(PValue))) + #volcanoplot with logFC versus pvalue
  geom_point() +
  theme_classic()+
  geom_vline(xintercept=c(-1.0, 0.8), col="orange", linetype="dashed") +
  geom_hline(yintercept=-log10(0.01), col="orange" , linetype="dashed")+
  geom_point(aes(col=diffexpressed)) + #add points colored by significance
  scale_color_manual(values= mycolors) + 
  ggtitle("25-50yo CA patients vs 50-75yo CA patients")
volc <- volc+geom_text_repel(data=head(input, 500), aes(label=gene), size = 3) #adding text for the top 20 genes
volc
#traceback()

# transform the normalized counts in log-counts-per-million
y <- cpm(filtered.counts, log=TRUE, prior.count = 1)

# select 10% differentially expressed genes and produce a heatmap
selY <- y[rownames(topGenes_Under50$table)[topGenes_Under50$table$FDR<0.10 & 
                                    abs(topGenes_Under50$table$logFC)>1.5],]

cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]

dev.off() 

cim(selY, margins = c(10, 10))


