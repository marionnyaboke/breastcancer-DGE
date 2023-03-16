library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(pathview)
library(GOSemSim)
library(ggnewscale)
library(europepmc)
library(fgsea)
BiocManager::install("ggplot2", force = TRUE)
install.packages("fgseaMultilevelCpp")
install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_1.0.5.tar.gz", repos = NULL, type="source")


# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# Prepare Input

{
# reading in data from EdgeR
df = read.csv("all-apoptotic-genes.csv", header=TRUE)

# we want the log fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- df$Gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
}

# Gene Set Enrichment

{
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             minGSSize = 5, 
             maxGSSize = 100, 
             pvalueCutoff = 1, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "fdr")
}

#Dotplot
{
require(DOSE)
dotplot(gse, showCategory=8, split=".sign") + facet_grid(.~.sign)
}

#Enrichment plot
{
  d <- godata('org.Hs.eg.db', ont="BP")
  ego2 <- pairwise_termsim(gse, method="Wang", semData = d)
  emapplot(ego2)
}

#Category netplot
#Depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network (helpful to see which genes are involved in enriched pathways and genes that may belong to multiple annotation categories)
{
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory =7)
}

#PubMed trend of enriched terms
{
  terms <- gse$Description[1:3]
  pmcplot(terms, 2010:2018, proportion=FALSE)
}

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
#keytypes(org.Hs.eg.db)

{

ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "SYMBOL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$Gene %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$logFC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
}

# Create gsee object
{
  kegg_organism = "hsa"
  kk2 <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = kegg_organism,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr",
                 keyType       = "ncbi-geneid")
}

#KEGG dotplot
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

