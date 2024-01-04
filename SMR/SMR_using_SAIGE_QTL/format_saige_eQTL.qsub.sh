## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -o stdout_generate_eQTL_input
#$ -e stderr_generate_eQTL_input
#$ -N generate_eQTL_input
#$ -q short.q
#$ -pe smp 1
#$ -l mem_requested=10G
#$ -r yes
#$ -t 1-14

cd /directflow/SCCGGroupShare/projects/angxue/proj/SAIGE-eQTL/SMR/eQTL_input

Rscript=/share/ScratchGeneral/angxue/software/R-4.0.5/bin/Rscript

ct=${SGE_TASK_ID};

touch ct_${ct}.log

$Rscript --vanilla format_saige_eQTL.R ${ct} >> ct_${ct}.log 

####
