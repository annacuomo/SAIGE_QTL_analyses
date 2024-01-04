##

new_names=(BimmNaive Bmem CD4all CD4effCM CD4TGFbStim CD8all CD8eff CD8unknown DC MonoC MonoNC NKact NKmat Plasma)

for k in {0..13}
do

mkdir -p ./smr_data/${new_names[k]}/

# rm -r -f ./${new_names[k]}/

sed 's/ct_name/'${new_names[k]}'/g' convert_besd_format.qsub.sh > convert_besd_format_${new_names[k]}.qsub.sh

qsub convert_besd_format_${new_names[k]}.qsub.sh

rm -f convert_besd_format_${new_names[k]}.qsub.sh

done

##
