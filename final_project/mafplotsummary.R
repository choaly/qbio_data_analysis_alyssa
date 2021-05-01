if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

#------------Recreate the Top 10 mutated genes MAF Summary Plot (so we can pick the top 2 mutated genes to investigate)------------#
mutation <- GDCquery_Maf(tumor = "BRCA",save.csv=TRUE, pipeline="varscan2") # only need to query once
maf_dataframe = read.maf(mutation)
#See here: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/ for information on each column

#to save the pdf, execute all your plotting code after the pdf() line:
pdf("final_project/maf_summary.pdf")
plotmafSummary(maf = maf_dataframe, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

# top 2 mutated genes are PIK3CA and TP53, and TTN
