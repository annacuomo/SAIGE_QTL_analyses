"""
Takes as input the directory containing all input files
creates a data frame with a row per file, 
extracts gene names from file names as an additional column
"""

library(data.table)

args = commandArgs(trailingOnly=TRUE)
input_dir = args[1]
output_file = args[2]

files = list.files(input_dir)
df = data.frame(genes = gsub(".cis.geno.pheno.txt","",files), files = files)
fwrite(df, output_file, sep="\t", col.names = FALSE)
