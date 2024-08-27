library(data.table)

args = commandArgs(trailingOnly=TRUE)
input_dir = args[1]
output_file = args[2]

files = list.files(input_dir)
df = data.frame(genes = gsub(".cis.geno.pheno.txt","",files), files = files)
fwrite(df, output_file, sep="\t", col.names = FALSE)
