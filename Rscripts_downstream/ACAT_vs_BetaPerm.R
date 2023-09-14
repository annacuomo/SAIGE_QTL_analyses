# this script compares performance of using the ACAT p-value 
# combination test as proposed in APEX vs the Beta approximation
# based on permutation proposed in the Fast QTL paper
# to obtain gene-level p-values in an eQTL study

# load R libraries
library(cowplot)
library(data.table)
library(ggplot2)

# Angli's TensorQTL results
# using a slightly updated version of the OneK1K data
# * removing two low QC samples (tot n: 980)
# * removing from analysis samples (donor x cell type) with <5cells

# files copied from the following directory (in Brenner)
# Anglis_dir = "/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/scale_gene_wise/CD4_NC_cell_5_mean_mx/"
Anglis_dir = "/share/ScratchGeneral/anncuo/OneK1K/anglis_files/tensorqtl/"

# load file
# chromosome: 1
# cell type: CD4 NC 
myfile = paste0(Anglis_dir,"OneK1K_CD4_NC.sig_cis_qtl_pairs.chr1.csv")
df = read.csv(myfile, sep="\t")
nrow(df)

