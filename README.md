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
* chr4: 890,388
* chr5: 815,838
* chr6: 822,879
* chr7: 708,513
* chr8: 690,326
* chr9: 520,275
* chr10: 624,823
* chr11: 619,520
* chr12: 580,103
* chr13: 456,146
* chr14: 399,653
* chr15: 328,591
* chr16: 371,767
* chr17: 295,102
* chr18: 342,772
* chr19: 234,779
* chr20: 256,191
* chr21: 150,535
* chr22: 143,083
