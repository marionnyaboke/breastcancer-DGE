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
  enriched <- enrichr(c("GSN", "CFL2", "GNG11", "ELANE", "PPM1F", "TP53AIP1", "FOS", "BMX", "BCL2L2", "PRKAR2B", "AIFM2", "ARHGAP10", "TP63", "PIK3R1", "GAS2", "BCL6", "JUN", "NGFR", "PPP2R2B", "TNFRSF10D", "PPP2R1B", "GNG2", "TLR4", "IGF1", "GNG7", "CFLAR", "TFDP2", "VIM", "BOK", "HIPK2", "TNFSF12", "DAPK2", "MAP2K5", "AXIN2", "MAPK10", "UACA", "GCK", "RPS6KA2", "RPS6KA5", "TJP2", "Y_RNA", "GNG12", "TNFSF11", "ADD1", "PRKD1", "BNIP3L", "CARD8", "IRS1", "TRAF6", "MYC", "PPM1L", "Y_RNA", "OGT", "CRADD", "PTPN13", "SOS2", "PTPN11", "MAP3K8", "GNB5", "RPS6KA3", "ATM", "Y_RNA", "CHMP3", "DSG1", "PPTC7", "TNFRSF8", "PPP3CC", "SOS1", "EGF", "SATB1", "PPP3CB"), dbs)
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

