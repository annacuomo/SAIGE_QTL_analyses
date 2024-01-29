# SAIGE-QTL analyses

Analyses as part of the [SAIGE-QTL](https://github.com/weizhou0/qtl) paper writing.

* [SAIGE-QTL runners](saige_qtl_runners)
* [downstream analyses](Rscripts_downstream) (plotting, summarising) scripts
* (private) [input file making scripts](https://github.com/annacuomo/Notebooks_private/tree/main/scripts/saigeqtl_onek1k)

## Rare variant annotation

Using AnnoVar:

```bash
perl table_annovar.pl ../../Documents/Garvan/onek1k_data/imputed_genotypes_filter_vcf_r08/chr2.dose.filtered.R2_0.8.vcf.gz humandb/ -buildver hg19 -out ../../Documents/Garvan/onek1k_data/imputed_genotypes_filter_vcf_r08/chr2_annotated -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -vcfinput -polish
```

## Count number of variants

```bash
cd /directflow/SCCGGroupShare/projects/anncuo/OneK1K/genotypes/filter_vcf_r08/
bcftools stats chr1.dose.filtered.R2_0.8.vcf.gz
```

Number of SNPs: 
* chr1: 928,708
* chr2: 1,032,942
* chr3: 895,348
* chr4:
* chr5:
* chr6:
* chr7:
* chr8:
* chr9:
* chr10:
* chr11:
* chr12:
* chr13:
* chr14:
* chr15:
* chr16:
* chr17:
* chr18:
* chr19:
* chr20:
* chr21:
* chr22:
