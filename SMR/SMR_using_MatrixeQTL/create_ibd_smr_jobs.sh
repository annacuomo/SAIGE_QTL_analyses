##

# new_names=(B_IN B_MEM CD4_NC CD4_ET CD4_SOX4 CD8_NC CD8_ET CD8_S100B DC Mono_C Mono_NC NK_R NK Plasma)

new_names=(BimmNaive Bmem CD4all CD4effCM CD4TGFbStim CD8all CD8eff CD8unknown DC MonoC MonoNC NKact NKmat Plasma)

mkdir -p ibd

for k in {0..13}
do

mkdir -p ./ibd/${new_names[k]}/

# rm -r -f ./${new_names[k]}/

sed 's/ct_name/'${new_names[k]}'/g' main_ibd_smr.qsub.sh > main_ibd_smr_${new_names[k]}.qsub.sh

qsub main_ibd_smr_${new_names[k]}.qsub.sh

rm -f main_ibd_smr_${new_names[k]}.qsub.sh

done

##
