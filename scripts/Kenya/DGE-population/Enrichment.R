#install.packages("enrichR")

# By default human genes are selected otherwise select your organism of choice.
library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE

#  find the list of all available databases from Enrichr.
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

library(knitr)
if (websiteLive) kable(head(dbs[c(1:6),-4]))

# query enrichr for downregulated genes
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")
if (websiteLive) {
  enriched <- enrichr(c("IGF1R", "TNFSF11", "PSMD3", "PSMB8", "ACIN1", "KRT18", "PTHLH", "DAPK3", "MAPT", "AKT1", "ADCY1", "IKBKB", "NAIP", "H1-4", "MAPKAPK2", "PRKCZ", "SREBF1", "CHMP4B", "CREBBP", "HSPB1", "MAP2K2", "GNA14", "CDC37", "MAPKAPK3", "PLEC", "PIK3R2", "PPP2R3B", "ARHGDIB", "GNAS", "SQSTM1", "HLA-A", "BCL2L12", "UBC", "NUMA1", "H1-2", "FADD", "HSPA1A", "H1-0", "BBC3"), dbs)
  }

# view the results table (Mol. Funct.)
if (websiteLive) enriched[["GO_Molecular_Function_2021"]]

if (websiteLive) {
  x <- head(enriched[["GO_Molecular_Function_2021"]])
  x[,1] <- gsub("GO:", "GO_", x[,1])
  kable(x)
}

# Plot Enrichr GO-MF output
if (websiteLive) plotEnrich(enriched[[1]], showTerms =15, numChar = 30, 
                            y = "Count", orderBy = "P.value", 
                            title = "GO Molecular Function"
)

# view the results table (Cell comp..)
if (websiteLive) enriched[["GO_Cellular_Component_2021"]]

if (websiteLive) {
  x <- head(enriched[["GO_Cellular_Component_2021"]])
  x[,1] <- gsub("GO:", "GO_", x[,1])
  kable(x)
}

# Plot Enrichr GO-CC output
if (websiteLive) plotEnrich(enriched[[2]], showTerms = 15, numChar = 30, 
                            y = "Count", orderBy = "P.value", 
                            title = "GO Cellular Component"
)

# view the results table (Bio process)
if (websiteLive) enriched[["GO_Biological_Process_2021"]]

if (websiteLive) {
  x <- head(enriched[["GO_Biological_Process_2021"]])
  x[,1] <- gsub("GO:", "GO_", x[,1])
  kable(x)
}

# Plot Enrichr GO-BP output
if (websiteLive) plotEnrich(enriched[[3]], showTerms = 15, numChar = 30, 
                            y = "Count", orderBy = "P.value", 
                            title = "GO Biological Process"
)

# GSEA with Cluster profiler
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(pathview)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# reading in data from EdgeR
setwd("C:/Users/Marion/Desktop/DGE/DGE-Kenya")

df = read.csv("all-apoptotic-genes.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- df$Gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ALIAS", 
             nPerm = 5, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

