## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -o stdout_convert_besd_ct_name
#$ -e stderr_convert_besd_ct_name
#$ -N convert_besd_ct_name
#$ -q short.q
#$ -pe smp 1
#$ -l mem_requested=10G
#$ -r yes
#$ -t 1-22

cd /directflow/SCCGGroupShare/projects/angxue/proj/SAIGE-eQTL/SMR/eQTL_input/

CELLTYPE="ct_name"
chr=${SGE_TASK_ID};

smr_tool="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/tools/smr_Linux"

touch ./smr_data/${CELLTYPE}/${CELLTYPE}_chr${chr}_210211.log

${smr_tool} --eqtl-summary ./${CELLTYPE}/${CELLTYPE}_chr${chr}.tsv --matrix-eqtl-format --make-besd --out ./smr_data/${CELLTYPE}/${CELLTYPE}_chr${chr}_210211 >> ./smr_data/${CELLTYPE}/${CELLTYPE}_chr${chr}_210211.log


####
