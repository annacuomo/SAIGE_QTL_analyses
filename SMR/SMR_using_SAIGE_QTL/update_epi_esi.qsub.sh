## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -o stdout_update_esi_epi_ct_name
#$ -e stderr_update_esi_epi_ct_name
#$ -N update_esi_epi_ct_name
#$ -q short.q
#$ -pe smp 1
#$ -l mem_requested=10G
#$ -r yes
#$ -t 1-22

cd /directflow/SCCGGroupShare/projects/angxue/proj/SAIGE-eQTL/SMR/eQTL_input/

CELLTYPE="ct_name"
chr=${SGE_TASK_ID};

smr_tool="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/tools/smr_Linux"

${smr_tool} --eqtl-summary ./${CELLTYPE}/${CELLTYPE}_chr${chr}.tsv --matrix-eqtl-format --make-besd --out ./smr_data/${CELLTYPE}/${CELLTYPE}_chr${chr}_210211 

epi_dir="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/smr_analysis/smr_data/${CELLTYPE}/${CELLTYPE}_chr${chr}_210211.epi"
esi_dir="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/smr_analysis/smr_data/${CELLTYPE}/${CELLTYPE}_chr${chr}_210211.esi"

${smr_tool} --beqtl-summary ./smr_data/${CELLTYPE}/${CELLTYPE}_chr${chr}_210211 --update-epi ${epi_dir} 

${smr_tool} --beqtl-summary ./smr_data/${CELLTYPE}/${CELLTYPE}_chr${chr}_210211 --update-esi ${esi_dir}

####
