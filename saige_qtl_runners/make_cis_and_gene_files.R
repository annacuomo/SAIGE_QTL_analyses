library(data.table)

dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/input/Sept23/SCT/"

# contains info on genes, including chromosome, start and end
gene_loc_file = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/GeneLocations.tsv"
gene_loc_df = fread(gene_loc_file)

# non-gene columns from pheno cov file
drop_cols = c('individual','barcode','sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2')

# define window size, here 1Mb
w = 1000000

# inputs will be stored here
output_dir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/input/power_analysis_subsets/"

# pick one at random, as genes are the same for all, the subset is on the cells
pct=10

file = paste0(dir,"CD4_NC_sc_pheno_cov_",pct,"pct_subset.tsv")
df = fread(file)
genes = cols[!(cols %in% drop_cols)]

#### cis region files
cis_region_files_dir = paste0(output_dir,"cis_region_files/")

# loop over genes
for (gene in genes){
  if(nrow(gene_loc_df[gene_loc_df$gene_name == gene,])==0){next}
  region_df = data.frame(chrom = gene_loc_df[gene_loc_df$gene_name == gene,"seqid"]$seqid,
                     window_start = max(0,(gene_loc_df[gene_loc_df$gene_name == gene,"start"]$start-w)),
                     window_end = gene_loc_df[gene_loc_df$gene_name == gene,"end"]$end+w)
  colnames(region_df) <- c()
  cis_out_file = paste0(cis_region_files_dir, gene,"_cis_region_1Mb.txt")
  fwrite(region_df, cis_out_file, sep = "\t")
}

#### gene list files by chromosome
gene_list_files_dir = paste0(output_dir,"gene_lists/")

for (chrom in 1:22){
    chrom_locs = gene_loc_df[gene_loc_df$seqid == chrom,]
    chrom_genes = chrom_locs$gene_name[chrom_locs$gene_name %in% genes]
    genes_df = data.frame(gene = chrom_genes)
    colnames(genes_df) <- c()
    gene_out_file = paste0(gene_list_files_dir, "chrom", chrom ,"_genes.txt")
    fwrite(genes_df, gene_out_file, sep = "\t")
}
