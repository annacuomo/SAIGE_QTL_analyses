####

for(k in c("as","cd","ibd","sle","t1dm")){

new_names <- c("BimmNaive", "Bmem", "CD4all", "CD4effCM", "CD4TGFbStim", "CD8all", "CD8eff", "CD8unknown", "DC", "MonoC", "MonoNC", "NKact", "NKmat", "Plasma", "Erythrocytes","Platelets")

res=list()
for(i in 1:14){
	a=c()
	for(j in 1:22){
		tmp=read.table(paste0("./",k,"/",new_names[i],"/",new_names[i],"_chr",j,".smr"),header=T)
		a=rbind(a,tmp)
	}
	res[[i]]=a
}

## Extract the significant signals

sig = c()
all = c()
for(i in 1:14){

tmp = res[[i]]
tmp$Cell_type = new_names[i]
all = rbind(all, tmp)

tmp = tmp[tmp$p_SMR < (0.05/nrow(tmp)),]
sig = rbind(sig,tmp)

print(nrow(tmp))
}

write.table(all, paste0(k, "_all_14_celltypes_all_smr.txt"),row.names=F,col.names=T,quote=F)
write.table(sig, paste0(k, "_all_14_celltypes_sig_smr.txt"),row.names=F,col.names=T,quote=F)

print(k)
}


####
