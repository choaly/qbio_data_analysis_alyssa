#Installs TCGAbiolinks if not already present
if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

setwd("Desktop/qbioresearch/qbio_data_analysis_alyssa")

#------------Recreate the Top 10 mutated genes MAF Summary Plot (so we can pick the top 3 mutated genes to investigate)------------#
mutation <- GDCquery_Maf(tumor = "BRCA",save.csv=TRUE, pipeline="varscan2") # only need to query once
maf_dataframe = read.maf(mutation)
#See here: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/ for information on each column

#Remember to use commented code when on the cluster and then use rsnyc to copy to local machine

#to save the pdf (execute all your plotting code after the pdf() line:
pdf("final_project/maf_summary.pdf")
plotmafSummary(maf = maf_dataframe, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

# top 3 mutated genes are PIK3CA, TP53, and TTN
#-----------------------------------------------------------------------------------------------------------------------------------#


#----------RNA–protein Spearman correlation analysis (working on--Alyssa)--------#



# #add barcodes argument to query if you want to run on your local machine for smaller files downloaded
# barcodes_rnaseq <- c("TCGA-BH-A0DG-01A-21R-A12P-07","TCGA-A2-A0YF-01A-21R-A109-07",
#                      "TCGA-AN-A04A-01A-21R-A034-07","TCGA-AR-A1AP-01A-11R-A12P-07",
#                      "TCGA-A2-A0T3-01A-21R-A115-07", "TCGA-E2-A154-01A-11R-A115-07" )
# barcodes_clinic <- c("TCGA-BH-A0DG","TCGA-A2-A0YF","TCGA-AN-A04A","TCGA-AR-A1AP", "TCGA-A2-A0T3",
#                      "TCGA-E2-A154", "TCGA-AO-A12F", "TCGA-A2-A0YM", "TCGA-BH-A0E0", "TCGA-AN-A0FL")

library(survival)
library(survminer)
library(arsenal)

clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type="xml")
setwd("data/experiments/2021-02-24_Week_4")
# GDCdownload( clin_query ) #should only need this command once. This downloads the files onto your system.
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"

#Add a new column to clinic called "age_category"
age_clinical = clinic$age_at_initial_pathologic_diagnosis
clinic$age_category = ifelse(age_clinical < 40, "Young", ifelse(age_clinical >= 60, "Old", "Mid"))



#use tableone to create summary of clinic data
# Add the age_category column in the same way as above
# How many patients are in subtypes vs clinic? Why?
# install.packages("tableone")
library(tableone)
clinic_summary <- CreateTableOne(data = clinic)
subtypes <- TCGAquery_subtype(tumor = "BRCA")

age_subs <- subtypes$age_at_initial_pathologic_diagnosis
subtypes$age_category = ifelse(age_subs < 40, "Young", ifelse(age_subs >= 60, "Old", "Mid"))



#name(subtypes) = make.names(names(subtypes))
names(subtypes)[names(subtypes) == "mRNA Clusters"] = "mRNA_Clusters"
table_arse <- tableby(age_category ~ (pathologic_stage) +(mRNA_Clusters) +(BRCA_Pathology),
                      data = subtypes, numeric.test="kwt", cat.test="chisq",
                      numeric.stats = c("Nmiss", "meansd"), total=FALSE)
df <- as.data.frame( summary(table_arse, text=TRUE, pfootnote=TRUE) )
write.csv(df, 'subtypes_summarize_by_age.csv', row.names=FALSE)


# Use rsync to copy figure onto local system and view
overall_survival <- as.integer( ifelse( is.na(clinic$days_to_death), clinic$days_to_last_follow_up, clinic$days_to_death) )
clinic$overall_survival <- overall_survival
clinic$death_event <- ifelse(clinic$vital_status == "Alive", 0,1)
# colnames(clinic) #use if want to visually see affect of ^^
cox_fit <- coxph(Surv(overall_survival, death_event)~age_at_initial_pathologic_diagnosis, data=clinic)
jpeg("cox_plot_age_continuous.jpg")
ggadjustedcurves(cox_fit, data=clinic)
dev.off()
