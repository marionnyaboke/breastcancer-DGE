
library("sf")
library(ggplot2)
library(ggVennDiagram)


# use list as an input
x <-list('A.American'=c("BMX", "GNG11", "PIK3R1", "PRKAR2B", "GNG2", "TLR4", "IGF1", "NGFR", "ITGA6"),
         'Caucasian'=c("ELANE", "TP63", "TP53AIP1", "FOS", "GSN", "Y_RNA", "PRKAR2B", "BMX", "NGFR", "GNG11", "TNFSF11", "IGF1", "AIFM2", "TNFRSF10D", "PPP2R2B", "CFL2", "GAS2", "MAPK10", "GNG2", "JUN", "DSG1", "GNG7", "ARHGAP10", "PIK3R1", "BCL6", "GCK", "TLR4", "BOK", "PPP2R1B", "IRS1", "AXIN2", "EGF", "VIM", "BCL2L2"),
         'Kenyan'=c("IGF1R", "TNFSF11", "PTHLH", "MAPT", "PSMB8", "KRT18", "HSPA1A", "PSMD3", "HSPB1", "SQSTM1", "GNA14", "PPP2R3B", "ADCY1", "DAPK3", "HLA-A"))
         
#UP-REG
#'Caucasian'=c("BIRC5", "GNG4", "CDC25C", "CDK1", "CDKN2A", "CLSPN", "LMNB1", "CDC25A", "BCL2A1", "E2F1", "PPP2R2C", "UNC5A", "GZMB", "TNFRSF18", "PSMD3", "TNFSF4", "TNFSF18", "PSMB3", "BCL2L14", "UNC5B", "ACTA1", "PMAIP1", "CDC25B", "CD3D", "HMGB2", "LMNB2", "RIPK2", "CDK5R1", "CD27", "TNFRSF4", "CD3E", "TNFRSF13B", "STK26", "HLA-A", "NFKBIE", "MAPK15", "PSME2", "IL18", "KRT18", "HSPB1", "LIMK1", "PPP1CA", "CDK5", "FASLG", "PARP1", "BIK", "HRK", "PRDX1", "BAK1", "CD3G", "HLA-C", "BCL2L12", "HLA-DQA1", "TNFSF13B", "HLA-DOB", "MAGED1", "IL1A", "BCL2L15", "TNFRSF12A", "TNFSF15", "PSMB9", "SFN", "CYCS", "PSMB4", "HLA-B", "IKBKE", "BMF", "CD70", "YWHAZ", "PSMD11", "TNF", "HLA-F", "RPS6KB2", "PSMA5", "PSMD14", "GNAS", "PTRH2", "PSMC4", "BAX", "CASP6", "MAPKAPK2", "SMPD3", "PIK3R3", "IRF1", "PPM1G", "RELB", "PSMB8", "BBC3", "CASP3", "PSMB2", "DAPK1", "PSMA7", "PSME4", "PRDX2", "SLC25A5", "TRAF2", "HLA-DRB1", "PTPA", "HLA-DQA2"),
#'Kenyan'=c("RNF103-CHMP3", "DSG1", "CHMP3", "TNFRSF6B", "ACTA1", "DSG3", "STK26", "HRK", "CDK1", "PMAIP1", "DCC", "BIRC3", "IL1B", "GNB4"),
#'A.American'=c("E2F1", "BIRC5", "LMNB2", "CDK1", "TNFRSF12A", "BAX", "CDC25C", "DAPK3", "CDK5", "PSMA7", "BCL2L12", "PPP1CA", "STUB1", "PSMD4", "GNB2", "PPP2R2C", "BBC3", "TRAF2", "RPS6K")

#DOWN-REG
#'Caucasian'=c("ELANE", "TP63", "TP53AIP1", "FOS", "GSN", "Y_RNA", "PRKAR2B", "BMX", "NGFR", "GNG11", "TNFSF11", "IGF1", "AIFM2", "TNFRSF10D", "PPP2R2B", "CFL2", "GAS2", "MAPK10", "GNG2", "JUN", "DSG1", "GNG7", "ARHGAP10", "PIK3R1", "BCL6", "GCK", "TLR4", "BOK", "PPP2R1B", "IRS1", "AXIN2", "EGF", "VIM", "BCL2L2")
#'Kenyan'=c("IGF1R", "TNFSF11", "PTHLH", "MAPT", "PSMB8", "KRT18", "HSPA1A", "PSMD3", "HSPB1", "SQSTM1", "GNA14", "PPP2R3B", "ADCY1", "DAPK3", "HLA-A")
#'A.American'=c("BMX", "GNG11", "PIK3R1", "PRKAR2B", "GNG2", "TLR4", "IGF1", "NGFR", "ITGA6")
# creating Venn diagram and displaying 
# all sets

ggVennDiagram(x) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
  scale_color_manual(values = c("Kenyan" = "black", "Caucasian" ="black", 'A.American'="black")) +
  scale_fill_gradient(low = "#D1EEEE", high = "#528B8B")
