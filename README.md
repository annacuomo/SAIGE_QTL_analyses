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
