miR_QTL=fread(sprintf('%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/SupTable4A_mirQTLs_raw.tsv',EVO_IMMUNO_POP))

isomiR_QTL=fread(sprintf('%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/SupTable4C_isomirQTLs_raw.tsv',EVO_IMMUNO_POP))

snps_data = fread(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/08.snp_data/general_informations/all_snps.tsv", sep=""))
snps_data[, snp_name := snp.name]
snps_data[, snp.name := NULL]

snp_data=fread(paste(EVO_IMMUNO_POP, "/Maxime/Evo_Immuno_pop_data/SNP_annotations/Map_imputed_essential_informations.txt", sep=""))
snp_data=merge(snps_data,snp_data)

miR_snps=snps_data[(MAF_EUB>.05 | MAF_AFB>.05) & (in_5p_miRNA | in_3p_miRNA)]


rankTransform = function(x){
        percentile=rank(x,ties.method='average')/(length(x)+1)
        mean_level=mean(x,na.rm=T)
        sd_level=sd(x,na.rm=T)
        qnorm(percentile,mean_level,sd_level)
        }




condIndex=c("NS","LPS","PAM3CSK4","R848","IAV")
options(stringsAsFactors=FALSE,max.print=9999)

pop='ALL'
CisDist=1e6
pvCis=0.05
minFreq=0.05
library(snpStats)
library(MatrixEQTL)

SampleAnnotFile=paste(EVO_IMMUNO_POP,'/Maxime/Evo_Immuno_pop_data/SampleAnnot.txt',sep='')
SampleAnnot=fread( SampleAnnotFile )
SampleAnnot=as.data.frame(SampleAnnot)

library(data.table)

Gene_annot=as.data.frame(fread(sprintf('%s/Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/GeneAnnot_expressed.txt',EVO_IMMUNO_POP)))
colnames(Gene_annot)[1:5]=c('gene_id','seqnames','strand','start','end')
Gene_annot$strand=ifelse(Gene_annot$strand>0,'+','-')

#RESperm = rbindlist(RESperm)
#fwrite(RESperm, file=sprintf('%s/Maxime/miRNA_V2/data/16_trans_effects/perm/allPerm_miRQTL_transEffects.txt',EVO_IMMUNO_POP))

G2S=function(x){Gene_annot[match(x,Gene_annot[,1]),'Associated Gene Name']}

trans_0=fread(sprintf('%s/Maxime/miRNA_V2/data/16_trans_effects/perm/miRQTL_transEffects_%s.txt',EVO_IMMUNO_POP,0))
trans_0$symbol=G2S(trans_0$gene)
miR_QTLs=rbind(miR_QTL[,.(snps,miRNA)],isomiR_QTL[,.(snps,miRNA)])
trans_0_target=merge(trans_0,miR_QTLs[which(!duplicated(miR_QTLs)),],by='snps',allow.cartesian=TRUE)


miRNA_BS=fread(sprintf('%s/Maxime/miRNA_V2/data/19_miR_targets/Comparison_5methods.txt',EVO_IMMUNO_POP))
miRNA_BS=miRNA_BS[,.(target_Score=max(target_Score)),by=.(miRNA,Symbol,algo)]

miRNA_BS_count=miRNA_BS[,.(NbAlgo=.N),by=.(miRNA,Symbol)]

colnames(trans_0_target)[colnames(trans_0_target)=='symbol']='Symbol'
trans_0_target_count=merge(trans_0_target,miRNA_BS_count,by=c('miRNA','Symbol'),all.x=TRUE,fill=0)
trans_0_target_count[is.na(NbAlgo),NbAlgo:=0]
trans_0_target_count=unique(trans_0_target_count[order(snps,gene,-NbAlgo,cond),],by=c('snps','gene','cond'))
trans_0_target_count[,.N,by=NbAlgo][order(NbAlgo),]
#    NbAlgo       N
# 1:      0 5863655
# 2:      1  635595
# 3:      2  397190
# 4:      3  166775
# 5:      4  223990
# 6:      5  121305

library(qvalue)
getBootLow=function(pvalue){
	quantile(replicate(1000,1-pi0est(p = sample(pvalue, replace=TRUE))$pi0),0.025)
}
getBootHigh=function(pvalue){
	quantile(replicate(1000,1-pi0est(p = sample(pvalue, replace=TRUE))$pi0),0.975)
}

pi0estim=trans_0_target_count[,.(pi1=1-pi0est(p = pvalue)$pi0,pi1_low=getBootLow(pvalue),pi1_high=getBootHigh(pvalue)),by=NbAlgo][order(NbAlgo),]
fwrite(pi0estim,file=sprintf('%s/Maxime/miRNA_V2/data/16_trans_effects/pi0estim_byNbAlgo.txt',EVO_IMMUNO_POP))

pi0estim_2=trans_0_target_count[,.(pi1=1-pi0est(p = pvalue)$pi0,pi1_low=getBootLow(pvalue),pi1_high=getBootHigh(pvalue)),by=NbAlgo>0]
pi0estim_2
fwrite(pi0estim_2,file=sprintf('%s/Maxime/miRNA_V2/data/16_trans_effects/pi0estim_anyAlgo.txt',EVO_IMMUNO_POP))

#  pi0estim_2
#    NbAlgo         pi1     pi1_low    pi1_high
# 1:  FALSE 0.004629030 0.002337125 0.006853294
# 2:   TRUE 0.007126305 0.002362799 0.011686766

trans_0_target_count[,fisher.test(table(pvalue<.05,NbAlgo>0))]
# 
# 	Fisher's Exact Test for Count Data
# 
# data:  table(pvalue < 0.05, NbAlgo > 0)
# p-value = 0.5727
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.9896679 1.0057604
# sample estimates:
# odds ratio 
#  0.9976719 



pdf(sprintf('%s/Maxime/miRNA_V2/figures/Revisions/trans_effect_targets.pdf',EVO_IMMUNO_POP),height=2,width=3)
#pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/SourcesOfpopDE_NS_Enrichment.pdf',EVO_IMMUNO_POP),height=5,width=1.8)

p <- ggplot(pi0estim,aes(NbAlgo,pi1,col=NbAlgo)) + geom_pointrange( ymin = pi0estim$pi1_low,ymax= pi0estim$pi1_high) + ylim(c(0,.035)) #+ scale_color_manual(values=colors)
p <- p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),axis.text.y = element_text(angle = 90, hjust = .5, vjust=0.5)) 
p <- p + geom_hline(aes(yintercept=0)) + xlab('')+ theme(legend.position = "none") +  + geom_hline(aes(yintercept=0))
p <- p + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+xlab('Number of database where the interaction is present')
print(p)
dev.off()

pi1=c()
for (perm in 0:100){
    cat(perm,'')
    trans_0=fread(sprintf('%s/Maxime/miRNA_V2/data/16_trans_effects/perm/miRQTL_transEffects_%s.txt',EVO_IMMUNO_POP,perm))
 
}

pi1
#  [1] 5.149773e-03 6.538302e-03 0.000000e+00 1.467121e-03 0.000000e+00
#  [6] 0.000000e+00 9.853807e-04 4.392215e-03 2.087316e-03 1.197675e-02
# [11] 0.000000e+00 0.000000e+00 4.614022e-03 0.000000e+00 0.000000e+00
# [16] 0.000000e+00 0.000000e+00 1.618161e-02 8.102693e-03 1.372052e-03
# [21] 0.000000e+00 2.331855e-02 0.000000e+00 4.563552e-03 9.423678e-04
# [26] 1.584876e-03 2.667240e-03 1.409373e-03 0.000000e+00 5.097479e-04
# [31] 0.000000e+00 9.719979e-03 3.318294e-03 0.000000e+00 1.445233e-02
# [36] 1.354946e-02 0.000000e+00 1.723956e-02 1.605995e-02 2.387260e-03
# [41] 1.264008e-02 0.000000e+00 6.449587e-03 1.426385e-03 0.000000e+00
# [46] 7.472382e-03 1.804666e-03 2.669939e-02 1.362624e-02 1.505778e-02
# [51] 0.000000e+00 1.357643e-02 0.000000e+00 1.132600e-02 0.000000e+00
# [56] 6.564793e-03 0.000000e+00 5.528671e-03 2.045293e-03 1.438578e-03
# [61] 0.000000e+00 3.755734e-03 7.279908e-03 5.112682e-03 0.000000e+00
# [66] 0.000000e+00 0.000000e+00 1.732266e-02 0.000000e+00 2.133541e-02
# [71] 1.957470e-03 0.000000e+00 2.026725e-02 6.214139e-03 1.043848e-02
# [76] 0.000000e+00 4.349861e-03 1.334765e-02 5.680720e-03 2.779909e-05
# [81] 0.000000e+00 5.833369e-03 1.337938e-02 2.288531e-02 5.012387e-03
# [86] 1.792375e-02 7.357817e-03 0.000000e+00 0.000000e+00 0.000000e+00
# [91] 6.249598e-03 0.000000e+00 3.878577e-03 0.000000e+00 0.000000e+00
# [96] 0.000000e+00 5.340314e-03 0.000000e+00 0.000000e+00 1.520660e-02
#[101] 0.000000e+00
> mean(pi1[-1]>pi1[1])
#[1] 0.36

