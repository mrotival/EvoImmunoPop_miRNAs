# nperm=1
#  df=data.frame(chr=rep(rep(1:22,(nperm+1)),e=5),cond=rep(1:5,(nperm+1)*22),perm=rep(0:nperm,each=22*5))
#  write.table(df,file='/pasteur/homes/mrotival/JobOutput/JobList_sQTL.txt',sep='\t',quote=F,row.names=F,col.names=F)

# df=data.frame(cond=rep(1:5,11),perm=rep(0:10,each=5))
# write.table(df,file='/pasteur/homes/mrotival/JobOutput/JobList_ds.txt',sep='\t',quote=F,row.names=F,col.names=F)



#### extract list of 22264 imputed variants that are multi-allelic (multiple annotated variants at the same position)
# these variants have been removed from Map_imputed bringing it down from 19619457 SNPs to 19597193 

QTLtype='exonRatioQTL'

############# create eQTL list
library(data.table)

QTL_pvalues_allCHR=list()
Res_QTL=list()
for (perm in 1:100){
	cat('\n',perm,':')
	for (CHR in 1:22){
		cat(CHR,'')
		fileName = sprintf('%s/Maxime/evo_immuno_pop_QTLs/%s/permutations/perm%s/Cis-%s_ALL_allCond_chr%s_perm%s_bestP_perFeature.txt', EVO_IMMUNO_POP, QTLtype, perm, QTLtype, CHR, perm)
		QTL_pvalues = fread( fileName )
		if(perm>0){
		QTL_pvalues$perm=perm
			if(CHR==1){
				QTL_pvalues_allCHR[[perm]]=QTL_pvalues
			}else{
				QTL_pvalues_allCHR[[perm]]=rbind(QTL_pvalues_allCHR[[perm]],QTL_pvalues)
				}
			}
		if(perm==0){
		fileName = sprintf('%s/Maxime/evo_immuno_pop_QTLs/%s/permutations/perm%s/Cis-%s_ALL_allCond_chr%s_perm%s_asssoc.txt', EVO_IMMUNO_POP, QTLtype, perm, QTLtype, CHR, perm)
		RESCIS = fread( fileName )
		RESCIS_nodup=RESCIS[order(gene,condition,pvalue),]
		rm(RESCIS);gc()
		RESCIS_nodup=RESCIS_nodup[!duplicated(paste(gene, condition)),]
			Res_QTL[[CHR]]=RESCIS_nodup
			}
		}
	}
Res_QTL=rbindlist(Res_QTL)
QTL_pvalues_allCHR=rbindlist(QTL_pvalues_allCHR)

fwrite(Res_QTL,file=sprintf('%s/Maxime/evo_immuno_pop_QTLs/%s/Res_%s_cond_nodup.txt',EVO_IMMUNO_POP, QTLtype,QTLtype))
fwrite(QTL_pvalues_allCHR,file=sprintf('%s/Maxime/evo_immuno_pop_QTLs/%s/permutations/QTL_pvalues_allperms.txt',EVO_IMMUNO_POP, QTLtype))


QTL_pvalues_observed=list()
for (CHR in 1:22){
	cat(CHR,'')
	fileName = sprintf('%s/Maxime/evo_immuno_pop_QTLs/%s/permutations/perm%s/Cis-%s_ALL_allCond_chr%s_perm%s_bestP_perFeature.txt', EVO_IMMUNO_POP, QTLtype, perm, QTLtype, CHR, perm)
	QTL_pvalues_observed[[CHR]]=fread( fileName )
}
QTL_pvalues_observed=rbindlist(QTL_pvalues_observed)

QTL_pvalues_observed=dcast(QTL_pvalues_observed,feature_id+QTL_type~condition,value.var='pvalue')
QTL_pvalues_observed[,minP:=pmin(NS,LPS,PAM3CSK4,R848,IAV)]
QTL_pvalues_observed=QTL_pvalues_observed[order(minP),]

QTL_pvalues_allCHR=dcast(QTL_pvalues_allCHR,perm+feature_id+QTL_type~condition,value.var='pvalue')
QTL_pvalues_allCHR[,minP:=pmin(NS,LPS,PAM3CSK4,R848,IAV)]

QTL_pvalues_observed$Nexp=sapply(QTL_pvalues_observed$minP, function(pval){
	Nb=QTL_pvalues_allCHR[minP<=pval,.(Nb=length(feature_id)),by=perm]
	mean(Nb$Nb)
}
QTL_pvalues_observed[,Nobs=1:nrow(QTL_pvalues_observed)]
QTL_pvalues_observed[,FDR:=cummin(Nexp/Nobs)]

FDR_compute=function(pval){ y=approxfun(c(0,-log10(QTL_pvalues_observed$pvalue),500),c(pmin(-log10(c(1,QTL_pvalues_observed$FDR)),6),500))(-log10(pval)); 10^-y}

Res_QTL$FDR=FDR_compute(Res_QTL$pvalue)
Res_QTL=Res_QTL[Res_QTL$FDR<0.05,]

library(snpStats)
SNP_genos=list()
for (CHR in 1:22){
	cat(CHR,'')
	Geno=read.plink(paste(EVO_IMMUNO_POP,'/DATA_FREEZE/ERC_Main_Genotyping_24022015/Imputation/EvoImmunoPop_imputation_200x19619457_chr',CHR,'.bed',sep=''))
	GenoNum=as(Geno$genotype,'numeric') 
	SNP_genos[[CHR]]=GenoNum[,na.omit(match(unique(Res_QTL$snps),colnames(GenoNum)))]
}

SNP_genos=cbindlist(SNP_genos)
# change to minor allele number
SNP_genos=t(2-SNP_genos)
fwrite(data.frame(snp_id=rn(SNP_genos),SNP_genos),file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/SNP_genos_',QTLtype,'_bestSNPperCondition.txt',sep=''))










for(i in 1:5){
	w=match(paste(Res_eQTL_nodup$snps,Res_eQTL_nodup$gene,i),paste(Res_eQTL$snps,Res_eQTL$gene,Res_eQTL$cond))
	Res_eQTL_nodup[[paste('Pvalue',condIndex[i],sep='_')]]=Res_eQTL$pvalue[w]
	Res_eQTL_nodup[[paste('Beta',condIndex[i],sep='_')]]=Res_eQTL$beta[w]
}
#Res_eQTL_nodup= Res_eQTL_nodup[,c(1:2,11,3:5,13:25)]
colnames(Res_eQTL_nodup)[5]='MinPvalue'
write.table(Res_eQTL_nodup,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/eQTL_rerun_PopCombined_allCond.txt',sep=''),sep='\t',quote=F,row.names=F)
Res_eQTL_nodup=Res_eQTL_nodup[order(Res_eQTL_nodup$gene,Res_eQTL_nodup$pval),]
Res_eQTL_nodup=Res_eQTL_nodup[!duplicated(Res_eQTL_nodup$gene),]

write.table(Res_eQTL_nodup,file=paste(EVO_IMMUNO_POP,'/Maxime/Splicing/eQTL_rerun/eQTL_rerun_PopCombined.txt',sep=''),sep='\t',quote=F,row.names=F)
