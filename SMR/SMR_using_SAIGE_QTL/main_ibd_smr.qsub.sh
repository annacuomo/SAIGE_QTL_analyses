#$ -S /bin/bash
#$ -cwd
#$ -o stdout_smr_ct_name
#$ -e stderr_smr_ct_name
#$ -N smr_ct_name
#$ -q short.q
#$ -pe smp 1
#$ -l mem_requested=10G
#$ -r yes
#$ -t 1-22

cd /directflow/SCCGGroupShare/projects/angxue/proj/SAIGE-eQTL/SMR/saige_eQTL/ibd/ct_name/

CELLTYPE="ct_name"
chr=${SGE_TASK_ID};

smr_tool="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/tools/smr_Linux"
plink_dir="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/SMR_analysis/onek1kdata/plink_files_with_newIDs"
gwas_file_path="/directflow/SCCGGroupShare/projects/angxue/proj/SAIGE-eQTL/SMR/gwas_data/ibd.ma"
smr_path="/directflow/SCCGGroupShare/projects/angxue/proj/SAIGE-eQTL/SMR/eQTL_input/smr_data"

${smr_tool} --bfile ${plink_dir}/plink_chr${chr}_renamed_201010 --gwas-summary ${gwas_file_path} \
      --beqtl-summary ${smr_path}/${CELLTYPE}/${CELLTYPE}_chr${chr}_210211 \
      --out ${CELLTYPE}_chr${chr} 2>&1 >> ${CELLTYPE}_chr${chr}.log


####
