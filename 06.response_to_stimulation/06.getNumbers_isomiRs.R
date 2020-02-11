
# get isomiR annot

isomirs_Annot_nosubs=fread(sprintf("%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_nosubs_FULL.tsv",EVO_IMMUNO_POP))
isomirs_Annot_FULL=fread(sprintf("%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_FULL_V2.0.tsv",EVO_IMMUNO_POP))



luq(gsub('(hsa.*)_MIMAT.*','\\1',isomirs_DE[,miRNA]))
#435
luq(gsub('(hsa.*)_MIMAT.*','\\1',isomirs_DE[,isomir]))
# 2,027
luq(isomirs_DE[,isomir])
# 2,049


# compute isomiR ratios

isomir_ratio_cast=fread(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_ratios_aggregated_nosubs.GCRL_Batch_lane_corrected.tsv",EVO_IMMUNO_POP))
isomir_ratio_cast=isomir_ratio_cast[-grep('other',isomir_ID)]
dim(isomir_ratio_cast)
isomir_ratio_melted = melt(isomir_ratio_cast, id = c("mirID", "isomir_ID"),
                               variable.name = "sample",
                               value.name = "ratio")
isomir_ratio_melted[, individual := substr(sample, 1,6)]
isomir_ratio_melted[, condition := substr(sample, 8,8)]
isomir_ratio_melted[, population := substr(sample, 1,3)]
inverseNormalRankTransform=function(x){n=length(x);
									qnorm(frank(x,na.last=FALSE,ties.method='random')/(n+1),mean(x,na.rm=T),sd(x,na.rm=T))}
isomir_ratio_melted[, ratio_transformed := inverseNormalRankTransform(ratio),by=.(mirID,isomir_ID)]


# get differential expression

isomirs_DE_bayes=fread(sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/isomirs_differentially_expressed_in_conditions_Likelihoods_withBeta.tsv",EVO_IMMUNO_POP))
isomirs_DE_bayes[,is_canonical:= isomirs_Annot_nosubs[match(isomirs_DE_bayes$isomir, ID), is_cannonical]]

isomirs_DE=fread(sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/isomirs_differentially_expressed_in_conditions.tsv",EVO_IMMUNO_POP))
isomirs_DE[,is_canonical:= isomirs_Annot_nosubs[match(isomirs_DE$isomir, ID), is_cannonical]]
isomirs_DE[,condition_test:=factor(condition_test,levels=c("LPS","PAM3CSK4","R848","IAV"))]

# number of isomiR/miRNA we are working with 

luq(isomirs_DE[,miRNA])
# 473
luq(gsub('(hsa.*)_MIMAT.*','\\1',isomirs_DE[,miRNA]))
# 435
luq(gsub('(hsa.*)_MIMAT.*','\\1',isomirs_DE[,isomir]))
# 2,027

luq(isomirs_DE[,isomir])
# 2,049


# how many miRNA change their isomiRs (counting different loci as 1 or as several miRNAs)
luq(gsub('(hsa.*)_MIMAT.*','\\1',isomirs_DE[fdr<0.01,miRNA])) # 290
luq(isomirs_DE[fdr<0.01,miRNA]) # 316

# how many miRNA change their canonical isomiRs (counting different loci as 1 or as several miRNAs)
luq(gsub('(hsa.*)_MIMAT.*','\\1',isomirs_DE[fdr<0.01 & is_canonical,miRNA])) # 193, 67%
luq(isomirs_DE[fdr<0.01 & is_canonical,miRNA]) # 212 # 67%

# in the remainder of the script each locus is counted as 1 miRNA

isomirs_DE[fdr<0.01 & is_canonical,.(NbCAN_diff=luq(miRNA),NbCAN_up=luq(miRNA[beta>0]),NbCANDown=luq(miRNA[beta<0]), PctCANDown=100*luq(miRNA[beta<0])/luq(miRNA)),by=condition_test]
#   condition_test NbCAN_diff NbCAN_up NbCANDown PctCANDown
#4:            LPS         60       26        34   56.66667
#3:       PAM3CSK4         62       21        41   66.12903
#1:           R848        162       54       108   66.66667
#2:            IAV        164       50       114   69.51220


colERC5=c("#525252AA", "#E31A1CAA", "#33A02CAA", "#1F78B4AA", "#6A3D9AAA")
condIndex=c("NS","LPS","PAM3CSK4","R848","IAV")

pdf(sprintf("%s/Maxime/miRNA_V2/figures/06.response_to_stimulation/isomirs_differentially_expressed_in_conditions_effectSize.pdf",EVO_IMMUNO_POP))
p <- ggplot(isomirs_DE[fdr<0.01,],aes(y=100*abs(beta),x=condition_test,fill=condition_test))
    p <- p + geom_violin() + scale_fill_manual(values=colERC5[-1]) 
    p <- p + scale_x_discrete(name="condition") + scale_y_sqrt(name ="abs(beta)",breaks=c(0,1,5,10,20,30,40,50))
    p <- p + geom_boxplot(fill="#FFFFFF88", outlier.size=0.5, notch=TRUE,width=0.4)
    p <- p + theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank() , legend.position = "none")
    p <- p + theme(axis.text.x = element_text(angle = 60, hjust = 1))
    print(p)
dev.off()


isomirs_DE[fdr<0.01,mean(abs(beta)>0.05),by= condition_test]
#   condition_test         V1
#4:            LPS 0.05214724
#2:       PAM3CSK4 0.06590258
#1:           R848 0.09988519
#3:            IAV 0.11070111


miR_modif_DE_all=fread(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/tests/stimulation_differences.txt",EVO_IMMUNO_POP))
tested_diff=c("is_3p_cannonical","is_3p_extended","is_3p_reduced","is_3p_template_canonical","is_3p_template_extended","is_3p_template_reduced","is_5p_canonical","is_5p_extended","is_5p_reduced","nta_3p_u_miR","nta_3p_a_miR")
miR_modif_DE_list=list()
miR_modif_DE=melt(miR_modif_DE_all[type%in%tested_diff, ], measure.vars=list(paste('Mean',condIndex[-1],sep='_'), paste('P',condIndex[-1],sep='_')), value.name = c("mean", "pvalue"),variable.name = "condition")
miR_modif_DE$condition=condIndex[as.numeric(miR_modif_DE$condition)+1]
miR_modif_DE$fdr=p.adjust(miR_modif_DE$pvalue)
miR_modif_DE$Delta=(miR_modif_DE$mean-miR_modif_DE$Mean_NS)*100
miR_modif_DE$Mean_NS=miR_modif_DE$Mean_NS*100

miR_modif_DE[fdr<0.01,][order(pvalue)]

cols=c("arm","type","Mean_NS","Delta_LPS","fdr_LPS","Delta_PAM3CSK4","fdr_PAM3CSK4","Delta_R848","fdr_R848","Delta_IAV","fdr_IAV")
miR_modif_DE_arm=dcast(miR_modif_DE,arm + type + Mean_NS~ condition, value.var = c("Delta", "fdr"))[,mget(cols)]
fwrite(miR_modif_DE_arm,file=sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/tests/SupTable2C_stimulation_differences.txt",EVO_IMMUNO_POP),sep='\t')

Pct_modif=fread(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/Pct_modifs_perSample.txt",EVO_IMMUNO_POP))

samples_names=unlist(fread(sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/sample_names_977_highQuality.tsv",EVO_IMMUNO_POP)))

i='is_3p_reduced'
Pct_3preduction=fread(sprintf('%s//Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/data/is_3p_reduced_allmiRNAs.txt',EVO_IMMUNO_POP))
ID=sort(unique(isomiR_annot$mirID))
condition = factor(condIndex[as.numeric(substr(samples_names, 8,8)),condIndex)
Pct_3preduction_3P=apply(Pct_3preduction[grep('-3p',ID),],2,mean)
Pct_3preduction_5P=apply(Pct_3preduction[grep('-5p',ID),],2,mean)
boxplot(Pct_3preduction_3P~condition)

boxplot(Pct_3preduction_5P~condition)



# check the isomiRs of miR-194-5p
x=isomirs_Annot_FULL[hsa_ID=='hsa-miR-194-5p', .(ID,mean_isoMiR_NS, mean_isoMiR_IAV)]
ggplot(melt(x),aes(x= ID,y=value,fill=variable))+geom_bar(stat="identity",position=position_dodge())+theme_bw()+theme(axis.text.x = element_text(angle = 90)).