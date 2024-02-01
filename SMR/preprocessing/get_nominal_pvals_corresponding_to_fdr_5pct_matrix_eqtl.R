library(data.table)

matrix_dir = "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/smr_analysis/matrix_eQTL/"

celltypes = c('BimmNaive','Bmem','CD4all','CD4effCM','CD4TGFbStim','CD8all','CD8eff',
              'CD8unknown','DC','MonoC','MonoNC','NKact','NKmat','Plasma')

acat_dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/genes/matrixeqtl_egenes/"

for (celltype in celltypes){
  ct_dir = paste0(matrix_dir, celltype, "/output_files/")
  acat_file = paste0(acat_dir, celltype,".csv")
  acat_df = fread(acat_file)
  # find genes with highest q-value still smaller than 0.05
  selected_df = acat_df[acat_df$qv==max(acat_df[acat_df$qv<0.05,]$qv),]
  if (nrow(selected_df) == 1){
    chrom = selected_df$chrom
    gene = selected_df$gene
    results_file = paste0(ct_dir,celltype,"_chr",chrom[2],"_cis_eqtls_210211.tsv")
    results_df = read.csv(gzfile(results_file), sep="\t")
    gene_df = results_df[results_df$gene == gene,]
    gene_df = gene_df[order(gene_df$p.value),]
    pval = head(gene_df,1)$p.value
    print(c(celltype, pval))
  } else {
      chroms = selected_df$chrom
      genes = selected_df$gene
      pvals = c()
      for (i in 1:length(genes)){
        results_file = paste0(ct_dir,celltype,"_chr",chroms[i],"_cis_eqtls_210211.tsv")
        gene_df = results_df[results_df$gene == genes[i],]
        gene_df = gene_df[order(gene_df$p.value),]
        pvals = c(pvals, head(gene_df,1)$p.value)
    }   
    pval = max(pvals)
    print(c(celltype, pval))
  }
}
