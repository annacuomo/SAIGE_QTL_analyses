# Summary of notebooks here

The notebooks in this folder contain code to plot and count the number of trans eGenes and eQTLs identified when we run SAIGE-QTL genome-wide for three representative cell types, CD4 NC (naive and central memory T cells), B IN (immature and naïve B cells), and Plasma.

* [pruning notebook](pruning_after_pv_filtering.ipynb), after p-value filtering.
* [trans eGenes & eQTLs counter notebook](count_trans_egenes_eqtls.ipynb) across the 3 cell types, all SNPs at MAF>=10%, 2 p-value cutoffs, `p<5e-6` and `p<5e-8`.
* trans eQTLs heatmap plotting (trans + cis as well as trans only), `p<5e-6`:
  * [CD4 NC heatmap notebook p<5e-6](trans_eqtl_heatmap_CD4_NC_p5e_6.ipynb)
  * [B IN heatmap notebook p<5e-6](trans_eqtl_heatmap_B_IN_p5e_6.ipynb)
  * [Plasma heatmap notebook p<5e-6](trans_eqtl_heatmap_Plasma_p5e_6.ipynb)
* trans eQTLs heatmap plotting (trans + cis as well as trans only), `p<5e-8`:
  * [CD4 NC heatmap notebook p<5e-8](trans_eqtl_heatmap_CD4_NC_p5e_8.ipynb)
  * [B IN heatmap notebook p<5e-8](trans_eqtl_heatmap_B_IN_p5e_8.ipynb)
  * [Plasma heatmap notebook p<5e-8](trans_eqtl_heatmap_Plasma_p5e_8.ipynb)
