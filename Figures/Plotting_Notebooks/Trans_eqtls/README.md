# Summary of notebooks here

The notebooks in this folder contain code to plot and count the number of trans eGenes and eQTLs identified when we run SAIGE-QTL genome-wide for three representative cell types, CD4 NC (naive and central memory T cells), B IN (immature and naïve B cells), and Plasma.

* trans eGenes & eQTLs counter (across the three cell types):
  * [counter notebook](count_trans_egenes.ipynb): all SNPs at MAF>=10%, p<5e-6.
  * [counter notebook genome-wide significance](count_trans_egenes_gw_significance.ipynb): filtered for SNPs at MAF>=10%, **p<5e-8**.
* trans eQTLs heatmap plotting (trans + cis as well as trans only):
  * [CD4 NC heatmap notebook](trans_eqtl_heatmap_CD4_NC.ipynb): trans results for CD4 NC.
  * [CD4 NC heatmap notebook genome-wide significance](trans_eqtl_heatmap_CD4_NC_p5e_8.ipynb): trans results for CD4 NC, filtered for SNPs at p<5e-8.
  * [B IN heatmap notebook](trans_eqtl_heatmap_B_IN.ipynb): trans results for B IN.
  * [B IN heatmap notebook genome-wide significance](trans_eqtl_heatmap_B_IN_p5e_8.ipynb): trans results for B IN, filtered for SNPs at p<5e-8.
  * [Plasma heatmap notebook](trans_eqtl_heatmap_Plasma.ipynb): trans results for Plasma.
  * [Plasma heatmap notebook genome-wide significance](trans_eqtl_heatmap_Plasma_p5e_8.ipynb): trans results for Plasma, filtered for SNPs at p<5e-8.
