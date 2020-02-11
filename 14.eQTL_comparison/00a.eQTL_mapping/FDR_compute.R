library(data.table)

QTL_pvalues=list()
for(perm in 0:100){
    cat(perm)
    for(CHR in 1:22){
        QTL_pvalues[[paste(CHR,perm)]]=fread(file=paste(EVO_IMMUNO_POP,'/Maxime/miRNA_V2/data/15.eQTL_comparisons/eQTL_mapping/',QTL_type,'/permutations/perm',perm,'/Cis-',QTL_type,'_ALL_allCond_chr',CHR,'_perm',perm,'_bestP_perFeature.txt',sep=''),sep='\t')
        QTL_pvalues[[paste(CHR,perm)]]$perm=perm
        }
    }
    QTL_pvalues=rbindlist(QTL_pvalues)
    QTL_pvalues=QTL_pvalues[,.(pvalue=min(pvalue)),by=.(feature_id,QTL_type,perm)]
    QTL_pvalues=QTL_pvalues[order(pvalue),]

QTL_pvalues$Nb_FP=cumsum(QTL_pvalues$perm>0)/length(setdiff(QTL_pvalues$perm,0))
QTL_pvalues$Nb_Pos=cumsum(QTL_pvalues$perm==0)
QTL_pvalues$FDR=rev(cummin(rev(QTL_pvalues$Nb_FP/QTL_pvalues$Nb_Pos)))

# let's check our FDR
sapply(1:100,function(i){mean(QTL_pvalues$perm[QTL_pvalues$FDR<0.05 & QTL_pvalues$perm%in%c(0,i)]>0)})

QTL_pvalues=QTL_pvalues[perm==0,]

fwrite(QTL_pvalues[perm==0,],file=paste(EVO_IMMUNO_POP,'/Maxime/miRNA_V2/data/15.eQTL_comparisons/eQTL_mapping/',QTL_type,'/FDR_estimates_',QTL_type,'.txt',sep=''),sep='\t')

perm=0
QTL_assoc=list()
for(CHR in 1:22){
    QTL_assoc[[paste(CHR,perm)]]=fread(file=paste(EVO_IMMUNO_POP,'/Maxime/miRNA_V2/data/15.eQTL_comparisons/eQTL_mapping/',QTL_type,'/permutations/perm',perm,'/Cis-',QTL_type,'_ALL_allCond_chr',CHR,'_perm',perm,'_asssoc.txt',sep=''),sep='\t')
    }
QTL_assoc=rbindlist(QTL_assoc)

FDR_compute=function(pval,pvalObs,FDRobs){y=approxfun(c(0,-log10(pvalObs),500),c(pmin(-log10(c(1,FDRobs)),6),500))(-log10(pval)); 10^-y}
QTL_assoc$FDR=FDR_compute(QTL_assoc$pvalue,QTL_pvalues$pvalue,QTL_pvalues$FDR)
QTL_assoc=QTL_assoc[order(gene,pvalue),]
QTL_assoc_bestSNP=QTL_assoc[!duplicated(gene),]
fwrite(QTL_assoc_bestSNP,file=paste(EVO_IMMUNO_POP,'/Maxime/miRNA_V2/data/15.eQTL_comparisons/eQTL_mapping/',QTL_type,'/',QTL_type,'_assoc_bestSNP.txt',sep=''),sep='\t')

QTL_assoc_bestSNP_cond=merge(QTL_assoc_bestSNP[FDR<0.05,mget(c('snps','gene'))],QTL_assoc,by=c('snps','gene'))
condIndex=c("NS", "LPS", "PAM3CSK4", "R848", "IAV")
QTL_assoc_bestSNP_cond$condition=factor(QTL_assoc_bestSNP_cond$condition,levels=condIndex)
fwrite(QTL_assoc_bestSNP_cond,file=paste(EVO_IMMUNO_POP,'/Maxime/miRNA_V2/data/15.eQTL_comparisons/eQTL_mapping/',QTL_type,'/',QTL_type,'_bestSNP_FDR5_5cond.txt',sep=''),sep='\t')

#QTL_assoc_bestSNP_cond=dcast(QTL_assoc_bestSNP_cond,snps+gene+CisDist~condition,value.var=c('beta','pvalue','FDR','statistic'))
#fwrite(QTL_assoc,file=paste(EVO_IMMUNO_POP,'/Maxime/evo_immuno_pop_QTLs/',QTL_type,'/',QTL_type,'.txt',sep=''),sep='\t')
