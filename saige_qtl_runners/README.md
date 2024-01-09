# Runners for SAIGE-QTL on OneK1K data

For the SAIGE-QTL method paper, we will evaluate performance on simulated data as well as real data.
For the latter, we will use the OneK1K dataset, as the current largest population-scale single-cell RNA-seq dataset.

## remember to add .csi index to VCF

```bash
module use /share/ClusterShare/apps/brenner/Modules/modulefiles
module load bcftools
bcftools index -c chr2.dose.filtered.R2_0.8.vcf.gz
```
