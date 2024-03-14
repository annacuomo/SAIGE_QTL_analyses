library(data.table)
library(qvalue)

# ACAT implementation functions (from Wei)
# Code adapted from the STAR package https://github.com/xihaoli/STAAR/blob/dc4f7e509f4fa2fb8594de48662bbd06a163108c/R/CCT.R wtih a modifitcaiton: when indiviudal p-value = 1, use minimum p-value 
#' An analytical p-value combination method using the Cauchy distribution
#'
#' The \code{CCT} function takes in a numeric vector of p-values, a numeric
#' vector of non-negative weights, and return the aggregated p-value using Cauchy method.
#' @param pvals a numeric vector of p-values, where each of the element is
#' between 0 to 1, to be combined.
#' @param weights a numeric vector of non-negative weights. If \code{NULL}, the
#' equal weights are assumed.
#' @return the aggregated p-value combining p-values from the vector \code{pvals}.
#' @examples pvalues <- c(2e-02,4e-04,0.2,0.1,0.8)
#' @examples CCT(pvals=pvalues)
#' @references Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test
#' with analytic p-value calculation under arbitrary dependency structures.
#' \emph{Journal of the American Statistical Association 115}(529), 393-402.
#' (\href{https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1554485}{pub})
#' @export

CCT <- function(pvals, weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }

  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }

  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  #if(is.zero && is.one){
  #  stop("Cannot have both 0 and 1 p-values!")
  #}
  if(is.zero){
    return(0)
  }
  if(is.one){
    #warning("There are p-values that are exactly 1!")
    return(min(1,(min(pvals))*(length(pvals))))
  }

  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }

  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}

get_CCT_pvalue = function(pvalue, weights=NULL){
   pvals = pvalue
   notna = which(!is.na(pvals))
   if(length(notna) > 0){
     pvals = pvals[!is.na(pvals)]
     cctpval = CCT(pvals, weights=weights)
   }else{
     cctpval = NA
   }
   return(cctpval)
}

saige_dir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/output/subsets/"
tensor_dir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/output/tensorqtl/"

sct_dir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K_from_ScratchGeneral/OneK1K/saige_eqtl/input/Sept23/SCT/"

subsets_analysis_dir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/input/power_analysis_subsets/"
gene_list_dir = paste0(subsets_analysis_dir,"gene_lists/")

out_dir = paste0(subsets_analysis_dir,"gene_level_results_atleast_1pct_expressed/")

# 1, 5, 10, 20, 50%
for (chrom in c(1:22)){
    genes_file = paste0(gene_list_dir,"chrom",chrom,"_genes.txt")
    genes_df = fread(genes_file, header=F)
    genes = genes_df$V1
    for (pct in c(1,5,10,20,50)){
        pheno_cov_filename_sc = paste0(sct_dir,"CD4_NC_sc_pheno_cov_",pct,"pct_subset.tsv")
        sc_df = fread(pheno_cov_filename_sc)
        sc_df = as.data.frame(sc_df)
        sc_mat = sc_df[,colnames(sc_df) %in% genes]
        n_cells = data.frame(gene = colnames(sc_mat), n_cell = colSums(sc_mat))
        n_cells$pct_expr = n_cells$n_cell / nrow(sc_mat)*100
        n_cells_1pct = n_cells[n_cells$pct_expr >= 1,]
        retained_genes = rownames(n_cells_1pct)
        tensor_file = paste0(tensor_dir,"CD4_NC_",pct,"pct_subset/CD4_NC_",pct,"pct_subset.cis_qtl_pairs.chr",chrom,".csv")
        tensor_df_all = as.data.frame(fread(tensor_file))
        df_gene = data.frame(gene = retained_genes)
        for (gene in retained_genes){
            # SAIGE-QTL
            saige_file = paste0(saige_dir,"subset_",pct,"pct/",gene,"_",pct,"pct_CD4_NC_poisson_test_cis")
            if (file.exists(saige_file) == FALSE){next}
            saige_df = as.data.frame(fread(saige_file))
            if (nrow(saige_df) == 0){next}
            saige_snps = unique(saige_df$MarkerID)
            # TensorQTL
            tensor_df = tensor_df_all[tensor_df_all$phenotype_id == gene,]
            tensor_snps = unique(tensor_df$variant_id)
            # get variants tested in both
            common_snps = saige_snps[saige_snps %in% tensor_snps]
            # get CCT values
            saige_acat_all = get_CCT_pvalue(saige_df$p.value)
            saige_acat_common = get_CCT_pvalue(saige_df[saige_df$MarkerID %in% common_snps,"p.value"])
            tensor_acat_all = get_CCT_pvalue(tensor_df$pval_nominal)
            tensor_acat_common = get_CCT_pvalue(tensor_df[tensor_df$variant_id %in% common_snps,"pval_nominal"])
            # add to data frame
            df_gene[df_gene$gene == gene,"saige_acat_all"] = saige_acat_all
            df_gene[df_gene$gene == gene,"saige_acat_common"] = saige_acat_common
            df_gene[df_gene$gene == gene,"tensor_acat_all"] = tensor_acat_all
            df_gene[df_gene$gene == gene,"tensor_acat_common"] = tensor_acat_common
      }
      df_gene = df_gene[rowSums(is.na(df_gene)) == 0, ] 
      df_gene$qv_saige_common = qvalue(df_gene$saige_acat_common, pi0=1)$qvalues
      df_gene$qv_tensor_common = qvalue(df_gene$tensor_acat_common, pi0=1)$qvalues
      df_gene$qv_saige_all = qvalue(df_gene$saige_acat_all, pi0=1)$qvalues
      df_gene$qv_tensor_all = qvalue(df_gene$tensor_acat_all, pi0=1)$qvalues
      out_file = paste0(out_dir, "chr",chrom,"_",pct,"pct.csv")
      fwrite(df_gene, out_file)
   }
}

# 100%
pct=100
for (chrom in c(1:22)){
    genes_file = paste0(gene_list_dir,"chrom",chrom,"_genes.txt")
    genes_df = fread(genes_file, header=F)
    genes = genes_df$V1
    # part1
    pheno_cov_filename_sc1 = paste0(sct_dir,"CD4_NC_sc_pheno_cov_part1.tsv")
    sc_df1 = fread(pheno_cov_filename_sc1)
    sc_df1 = as.data.frame(sc_df1)
    sc_mat1 = sc_df1[,colnames(sc_df1) %in% genes]
    n_cells1 = data.frame(gene = colnames(sc_mat1),
                        n_cell = colSums(sc_mat1))
    n_cells1$pct_expr = n_cells1$n_cell / nrow(sc_mat1)*100
    n_cells_1pct1 = n_cells1[n_cells1$pct_expr >= 1,]
    retained_genes1 = rownames(n_cells_1pct1)
    # part2
    pheno_cov_filename_sc2 = paste0(sct_dir,"CD4_NC_sc_pheno_cov_part2.tsv")
    sc_df2 = fread(pheno_cov_filename_sc2)
    sc_df2 = as.data.frame(sc_df2)
    sc_mat2 = sc_df2[,colnames(sc_df2) %in% genes]
    n_cells2 = data.frame(gene = colnames(sc_mat2),
                        n_cell = colSums(sc_mat2))
    n_cells2$pct_expr = n_cells2$n_cell / nrow(sc_mat2)*100
    n_cells_1pct2 = n_cells2[n_cells2$pct_expr >= 1,]
    retained_genes2 = rownames(n_cells_1pct2)
    # combine
    retained_genes = c(retained_genes1, retained_genes2)
    tensor_file = paste0(tensor_dir,"CD4_NC_",pct,"pct_subset/CD4_NC_",pct,"pct_subset.cis_qtl_pairs.chr",chrom,".csv")
    tensor_df_all = fread(tensor_file)
    df_gene = data.frame(gene = retained_genes)
    for (gene in retained_genes){
        # SAIGE-QTL
        saige_file = paste0(saige_dir,"subset_",pct,"pct/",gene,"_",pct,"pct_CD4_NC_poisson_test_cis")
        if (file.exists(saige_file) == FALSE){next}
        saige_df = fread(saige_file)
        saige_snps = unique(saige_df$MarkerID)
        # TensorQTL
        tensor_df = tensor_df_all[tensor_df_all$phenotype_id == gene,]
        tensor_snps = unique(tensor_df$variant_id)
        # get variants tested in both
        common_snps = saige_snps[saige_snps %in% tensor_snps]
        # get CCT values
        saige_acat_all = get_CCT_pvalue(saige_df$p.value)
        saige_acat_common = get_CCT_pvalue(saige_df[saige_df$MarkerID %in% common_snps,"p.value"]$p.value)
        tensor_acat_all = get_CCT_pvalue(tensor_df$pval_nominal)
        tensor_acat_common = get_CCT_pvalue(tensor_df[tensor_df$variant_id %in% common_snps,"pval_nominal"]$pval_nominal)
        # add to data frame
        df_gene[df_gene$gene == gene,"saige_acat_all"] = saige_acat_all
        df_gene[df_gene$gene == gene,"saige_acat_common"] = saige_acat_common
        df_gene[df_gene$gene == gene,"tensor_acat_all"] = tensor_acat_all
        df_gene[df_gene$gene == gene,"tensor_acat_common"] = tensor_acat_common
    }
    df_gene = df_gene[rowSums(is.na(df_gene)) == 0, ] 
    df_gene$qv_saige_common = qvalue(df_gene$saige_acat_common, pi0=1)$qvalues
    df_gene$qv_tensor_common = qvalue(df_gene$tensor_acat_common, pi0=1)$qvalues
    df_gene$qv_saige_all = qvalue(df_gene$saige_acat_all, pi0=1)$qvalues
    df_gene$qv_tensor_all = qvalue(df_gene$tensor_acat_all, pi0=1)$qvalues
    out_file = paste0(out_dir, "chr",chrom,"_",pct,"pct.csv")
    fwrite(df_gene, out_file)
}




