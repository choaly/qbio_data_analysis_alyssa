import cptac
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats 

# cptac.download(dataset="Brca")
br = cptac.Brca()

protein_data = br.get_proteomics()

#The dataframes are MultIndex pandas dataframes. 
#However, to teach the basics of pandas, we will remove the "multi" part of the dataframe.
protein_data = protein_data.droplevel(1, axis=1)

rna_data = br.get_transcriptomics()
clinical_data = br.get_clinical()

clinical_data["Age_in_years"] = clinical_data["Age.in.Month"]/12 #fill in with correct math


# clinical_data.to_csv("~/Desktop/qbioresearch/qbio_data_analysis_alyssa/data/clinical_data.csv")
# protein_data.to_csv("~/Desktop/qbioresearch/qbio_data_analysis_alyssa/data/protein_data.csv")
# rna_data.to_csv("~/Desktop/qbioresearch/qbio_data_analysis_alyssa/data/rna_data.csv")
# #Read back in your data
# #Try with different arguments. The index_col=0 creates the index (rownames) as the first column. 
# #Double check that the index and columns are what you want
# clinical_data_readin = pd.read_csv("~/Desktop/qbioresearch/qbio_data_analysis_alyssa/data/clinical_data.csv", index_col=1)
# clinical_data_readin


#Step 1: Check that the patients are in the same order.
#This is important because we are looking at the gene expression and protein information of EACH patient
#The data needs to be in pairs. This is easiest when the patients are in the same order in rna_data and protein_data

#Are the index or columns (rna_data.index or rna_data.columns) the patients?
#A: rna_data.index

#Write an assert statement that checks the patients rna_data are equal to the patients of protein_data 
#Fill in parantheses 
assert list(protein_data.index) == list(rna_data.index)



#Plot the data
#You can present the rho value in various ways. Table, scatterplot, etc. 
#Here you will make a scatterplot
def spearman_plots (gene_name, gene_rna_data, gene_protein_data, output_file):
    '''make scatterplots showing correlation between transcriptomic and proteomic data'''

    rho, spear_pvalue = stats.spearmanr( gene_rna_data, gene_protein_data ) #using scipy

    plt.figure( figsize=(10,10) ) #using matplotlib.pyplot
    plt.scatter( gene_rna_data, gene_protein_data )

    #trend line
    #https://stackoverflow.com/questions/41635448/how-can-i-draw-scatter-trend-line-on-matplot-python-pandas
    z = np.polyfit(gene_rna_data, gene_protein_data, 1)
    p = np.poly1d(z)
    plt.plot(gene_rna_data,p(gene_rna_data),"r--")

    title = "rho: {} for {}".format(rho, gene_name)
    plt.title(title)

    #Fill in informative x and y labels
    plt.xlabel("RNA") #how do we know RNA is x-axis, not y-axis?  plt.scatter( gene_rna_data, gene_protein_data )
    plt.ylabel("Protein")

    pngfile = '{}.png'.format( output_file)

    plt.savefig( pngfile, bbox_inches="tight" ) #Use this when saving figure in script (use plot.show() in jupyter ntbk)


#Access the transcriptomic and proteomic information of the specific gene (ESR1)

#Fill in the brackets []. We want ALL rows of the [specifc gene] column. 
#Remember [row,col] format and : refers to ALL.

rna_pik3ca = rna_data.loc[:, "PIK3CA"]
protein_pik3ca = protein_data.loc[:, "PIK3CA"]
spearman_plots("PIK3CA", rna_pik3ca, protein_pik3ca, "pik3ca_spear")

rna_tp53 = rna_data.loc[:, "TP53"]
protein_tp53 = protein_data.loc[:, "TP53"]
spearman_plots("TP53", rna_tp53, protein_tp53, "tp53_spear")



# #----------UNNECESSARY--DROPPED TTN. Test if TTN has null values-----------
# rna_ttn = rna_data.loc[:, "TTN"]
# protein_ttn = protein_data.loc[:, "TTN"]

# rna_data_readin = pd.read_csv("~/Desktop/qbioresearch/qbio_data_analysis_alyssa/data/cptac_data/rna_data.csv", index_col=1)
# protein_data_readin = pd.read_csv("~/Desktop/qbioresearch/qbio_data_analysis_alyssa/data/cptac_data/protein_data.csv", index_col=1)


# # # creating bool series True for NaN values
# # rna_bool_series = pd.isnull(rna_data_readin["TTN"])
 
# # # filtering data; displayind data only with team = NaN
# # rna_ttn_null = rna_data_readin[rna_bool_series]
# # rna_ttn_null.to_csv("rna_ttn_null.csv")

# rna_ttn_null = rna_ttn.isnull()
# protein_ttn_null = protein_ttn.isnull()

# # rna_ttn_null.to_csv("rna_ttn_null.csv")
# # protein_ttn_null.to_csv("protein_ttn_null.csv")
# #----------Test if TTN has null values-----------


# # print( "TTN RNA data length: ", len(rna_ttn), "\n TTN Protein data length: ", len(protein_ttn) )
# # #Write an assert statement that checks the patients rna_data are equal to the patients of protein_data 
# # assert len(protein_ttn) == len(rna_ttn)
# # spearman_plots("TTN", rna_ttn, protein_ttn, "ttn_spear")