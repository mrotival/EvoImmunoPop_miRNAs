###################################################################################
##Aim of the script be able to plot a visual explanation of a isomir ratio change##
###################################################################################

##############
##Libariries##
##############
require(ggplot2)

###########################
##### get GO info       ###
###########################

source(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/scripts/00a.GOannotation/01.Load_and_format_databases.R", sep=""))
GO_informations = load_BHF_UCL_miRNA()



#################################
##  load isomiR ratio Data     ##
#################################

isomir_ratio_cast=fread(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_ratios_aggregated_nosubs.GCRL_Batch_lane_corrected.tsv",EVO_IMMUNO_POP))
#isomir_ratio_cast=isomir_ratio_cast[-grep('other',isomir_ID)]
#dim(isomir_ratio_cast)
isomir_ratio_melted = melt(isomir_ratio_cast, id = c("mirID", "isomir_ID"),
                               variable.name = "sample",
                               value.name = "ratio")
isomir_ratio_melted[, individual := substr(sample, 1,6)]
isomir_ratio_melted[, condition := substr(sample, 8,8)]
isomir_ratio_melted[, population := substr(sample, 1,3)]
isomir_ratio_melted[,hsa_ID := gsub('(hsa.*)_(MIMAT.*)_(MI[0-9]+)','\\1',mirID)]

# get meanRatio NS
meanRatio=isomir_ratio_melted[condition==1,.(baseRatio=mean(ratio,na.rm=T)),by=isomir_ID]

#################################
##  load isomiR count Data     ##
#################################

isomir_count_cast=fread(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_counts_aggregated_nosubs.log2RPM.GCRL_Batch_lane_corrected.tsv",EVO_IMMUNO_POP))
isomir_count_melted = melt(isomir_count_cast, id = c("V1"),
                               variable.name = "sample",
                               value.name = "count")
isomir_count_melted[, individual := substr(sample, 1,6)]
isomir_count_melted[, condition := substr(sample, 8,8)]
isomir_count_melted[, population := substr(sample, 1,3)]
isomir_count_melted[, hsa_ID := gsub('.*(hsa.*)_(MIMAT.*)_(MI[0-9]+).*','\\1',V1)]
isomir_count_melted[, mirID := gsub('.*(hsa.*)_(MIMAT.*)_(MI[0-9]+).*','\\1_\\2_\\3',V1)]
isomir_count_melted[, isomir_ID := V1]


########
##Data##
########
isomir_annot=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_nosubs_FULL.tsv',EVO_IMMUNO_POP))

isomiR_DE=fread(sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/isomirs_differentially_expressed_in_conditions.tsv",EVO_IMMUNO_POP))
isomiR_DE[, canonical := isomir_annot[match(isomiR_DE$isomir, isomir_annot$ID), is_cannonical]]
isomiR_DE[, baseRatio:= meanRatio[match(isomiR_DE$isomir, meanRatio$isomir), baseRatio]]

isomiR_DE[, pbonf := p.adjust(pvalue, method = "bonferroni")]
isomiR_DE[, fdr := p.adjust(pvalue, method = "fdr")]
isomiR_DE[, minPval_miRNA := min(pvalue),by=miRNA]
isomiR_DE[, minPval_isomiR := min(pvalue),by=isomir]
isomiR_DE[, maxBeta_miRNA := max(abs(beta)),by=miRNA]
isomiR_DE[, maxBeta_isomiR := max(abs(beta)),by=isomir]
isomiR_DE=isomiR_DE[order(minPval_miRNA,pvalue,-beta)]
isomir_DE_diff= unique(isomiR_DE[fdr<.01 ,isomir],by='isomir')
isomir_DE_diff_miR= unique(isomiR_DE[fdr<.01 ,miRNA],by='miRNA')
luq(isomir_DE_diff_miR) # 316
isomir_DE_diff_up= unique(isomiR_DE[fdr<.01 & beta>0,isomir],by='isomir')
isomir_DE_diff_up_miR= unique(isomiR_DE[fdr<.01 & beta>0,miRNA],by='miRNA')

canonical_isomiR_DE = isomiR_DE[which(canonical),]
canonical_isomir_DE_diff_miR= unique(canonical_isomiR_DE[fdr<.01 ,miRNA],by='miRNA')
luq(canonical_isomir_DE_diff_miR) # 212


# isomirs_ratios=fread(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_ratios_aggregated_nosubs.GCRL_Batch_lane_corrected.tsv",EVO_IMMUNO_POP))
# isomirs_ratios=isomirs_ratios[-grep('other',isomir_ID)]
# dim(isomirs_ratios)

# isomirs_differential_expression[, start := isomir_annot[match(isomirs_differential_expression$isomir, isomir_ID), start]]
# isomirs_differential_expression[, end := isomir_annot[match(isomirs_differential_expression$isomir, isomir_ID), end]]


#301 miRNA were wa can do the test
condition_name = c("1"="NS", "2"="LPS", "3" = "PAM3CSK4", "4" ="R848", "5" = "IAV")
canonical_isomiR_DE[, condition := factor(condition_test,levels=condition_name)]

canonical_isomiR_DE[,.(signif=sum(fdr<.01),signifUp=sum(fdr<.01 & beta>0),signifDown=sum(fdr<.01 & beta<0),PctUp=sum(fdr<.01 & beta>0)/sum(fdr<.01),PctDown=sum(fdr<.01 & beta<0)/sum(fdr<.01)) ,by=condition]
# condition signif signifUp signifDown     PctUp
# 1:      R848    162       54        108 0.3333333
# 2:       IAV    164       50        114 0.3048780
# 3:  PAM3CSK4     63       22         41 0.3492063
# 4:       LPS     60       26         34 0.4333333

diff_canonical_isomiR=canonical_isomiR_DE[fdr<.01,isomir]
diff_canonical_miR=canonical_isomiR_DE[fdr<.01,miRNA]
Nb_isomiR=luq(canonical_isomiR_DE[,isomir])

isomiR_DE_lik = fread(sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/isomirs_differentially_expressed_in_conditions_Likelihoods_withBeta.tsv",EVO_IMMUNO_POP), colClasses=c("character","character","character","numeric","numeric","numeric","integer","numeric"))
isomiR_DE_lik[,bestModel:=condition_test]
isomiR_DE_lik[,condition_test:=NULL]

# isomiR_DE_lik = isomiR_DE_lik[isomir%in%diff_canonical_isomiR,]
nullModel = isomiR_DE_lik[bestModel=='00000',mget(c('isomir','miRNA','logLik'))][order(miRNA,isomir)]
bestModel = unique(isomiR_DE_lik[order(miRNA,isomir,-ProbModel)],by='isomir')
bestModel = merge(bestModel,nullModel,suffixes=c('','.null'),by=c('miRNA','isomir'),all.x=TRUE)
bestModel[,logDiff:=logLik-logLik.null]
bestModel[,LR_pvalue:=pchisq(2*logDiff,1,low=F)]
bestModel[,LR_fdr:=p.adjust(LR_pvalue,'fdr')]
bestModel=bestModel[,diff_isomiR:=isomir%in%isomir_DE_diff]
bestModel=bestModel[,diff_miR:=miRNA%in%isomir_DE_diff_miR]
# isomiR_DE_lik = isomiR_DE_lik[order(isomir,-ProbModel)]
# isomiR_DE_lik = isomiR_DE_lik[df>3]
# best_model=unique(isomiR_DE_lik,by= "isomir")


#######################################################################
######       table with all FC + best model 			           ####
#######################################################################

isomir_DE_cast=dcast(isomiR_DE,isomir+miRNA+canonical+baseRatio+minPval_miRNA+minPval_isomiR+maxBeta_miRNA+maxBeta_isomiR~condition_test,value.var=c('beta','pvalue','fdr'))

# add likelihood informations for DE miRs (other miRs are NA)

isomir_DE_table = merge(isomir_DE_cast,bestModel[which(diff_isomiR),mget(c("isomir","miRNA","bestModel","ProbModel","LR_fdr","diff_isomiR","diff_miR"))],by=c('isomir','miRNA'),all.x=TRUE)

cols=c("isomir","miRNA","canonical","baseRatio","minPval_miRNA","minPval_isomiR","beta_LPS", "pvalue_LPS","fdr_LPS", 
			"beta_PAM3CSK4", "pvalue_PAM3CSK4",   "fdr_PAM3CSK4",
			"beta_R848","pvalue_R848",  "fdr_R848",
			 "beta_IAV","pvalue_IAV", "fdr_IAV","bestModel","ProbModel")
isomir_DE_table=isomir_DE_table[order(minPval_miRNA,minPval_isomiR),mget(cols)]


isomir_DE_table2=isomir_DE_table


# order by minP
isomir_DE_table[,minPval_isomiR:=NULL]
isomir_DE_table[,minPval_miRNA:=NULL]


# write
fwrite(isomir_DE_table,file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable2B_isomirDE.tsv",EVO_IMMUNO_POP),sep='\t')

######## detail isomiR infos
isomir_DE_table2=merge(isomir_DE_table2,isomir_annot,by.x="isomir",by.y="ID",all.x=T)

cols2=c('hsa_ID','MIMAT','MI_ID','shift_5p','shift_3p','isomiR_subs','isomir_sequence',"canonical","baseRatio","minPval_miRNA","minPval_isomiR","beta_LPS", "pvalue_LPS","fdr_LPS", 
			"beta_PAM3CSK4", "pvalue_PAM3CSK4",   "fdr_PAM3CSK4",
			"beta_R848","pvalue_R848",  "fdr_R848",
			 "beta_IAV","pvalue_IAV", "fdr_IAV","bestModel","ProbModel",'isomir')

isomir_DE_table2=isomir_DE_table2[order(minPval_miRNA,minPval_isomiR),mget(cols2)]
isomir_DE_table2[,isomiR_subs:=gsub('T','U',isomiR_subs)]
isomir_DE_table2[,isomir_sequence:=gsub('T','U',isomir_sequence)]
isomir_DE_table2[,minPval_isomiR:=NULL]
isomir_DE_table2[,minPval_miRNA:=NULL]

fwrite(isomir_DE_table2,file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable2B_V2_isomirDE.tsv",EVO_IMMUNO_POP),sep='\t')




#############################################################################
####    create a plot with sharing of isomiR changes across conditions   ####
#############################################################################


Nb_isomiR=luq(isomir_DE_table$isomir)

allModels_STIM=cbind(NS=rep(0,8),
				LPS=rep(rep(0:1,e=8),1),
				PAM3CSK4=rep(rep(0:1,e=4),2),
				R848=rep(rep(0:1,e=2),4),
				IAV=rep(rep(0:1,e=1),8))
rownames(allModels_STIM)=apply(allModels_STIM,1,paste,collapse='')


code_values_ordered = c("01000","00100","00010","00001","01100","01010","01001",
                        "00110","00101","00011","01110","01101","01011","00111","01111")

# pdf(sprintf("%s/Maxime/miRNA_V2/figures/06.response_to_stimulation/miRNA_differentially_expressed_in_conditions_Likelihoods_part1.pdf",EVO_IMMUNO_POP),height=3,width=5)
# par(mar=c(4,7,1,1))
# Image(t(allModels_STIM[code_values_ordered,-1]))
# dev.off()

tab=table(isomir_DE_table[, bestModel])
tab=tab[code_values_ordered]
tab[is.na(tab)]=0
names(tab)=code_values_ordered

pdf(sprintf("%s/Maxime/miRNA_V2/figures/06.response_to_stimulation/isomir_differentially_expressed_in_conditions_Likelihoods_part2.pdf",EVO_IMMUNO_POP),height=3,width=5)
par(mar=c(4,8,3,3))
	allModels=allModels_STIM[code_values_ordered,-1]
    Prob_diff=sapply(2:5,function(i){sum(tab[substr(names(tab),i,i)==1])/Nb_isomiR})
    ExpectedPct=exp(allModels%*%log(matrix(Prob_diff,4))+(1-allModels)%*%log(1-matrix(Prob_diff,4)))
	x=barplot(tab,col='white',ylab='Model frequency',las=2)
	for (i in 2:5){
    	barplot(tab*allModels[,i-1],col=gsub('AA$','66',colERC5[i]),las=2,add=T,xaxt='n',yaxt='n')
	}
    x=barplot(tab,col='white',ylab='Model frequency',las=2)
	points(x,ExpectedPct*sum(tab),pch='-',col='grey',cex=3)
	for (i in 2:5){
    	barplot(tab*allModels[,i-1],col=gsub('AA$','66',colERC5[i]),las=2,add=T,xaxt='n',yaxt='n')
	}
    pbinom(as.numeric(tab),Nb_isomiR,as.numeric(ExpectedPct),low=F)
# 1.000000e+00 1.000000e+00 1.000000e+00 4.493919e-01 7.024789e-01 1.000000e+00 9.999934e-01 9.928960e-01 1.000000e+00 5.826827e-02 6.249472e-41 9.677411e-01 3.108674e-01 7.081362e-01 1.768179e-26
 dev.off()

tab=table(isomir_DE_table[which(canonical), bestModel])
tab=tab[code_values_ordered]
tab[is.na(tab)]=0
names(tab)=code_values_ordered


pdf(sprintf("%s/Maxime/miRNA_V2/figures/06.response_to_stimulation/canonical_isomir_differentially_expressed_in_conditions_Likelihoods_part2.pdf",EVO_IMMUNO_POP),height=3,width=5)
par(mar=c(4,8,3,3))
	allModels=allModels_STIM[code_values_ordered,-1]
    Prob_diff=sapply(2:5,function(i){sum(tab[substr(names(tab),i,i)==1])/Nb_isomiR})
    ExpectedPct=exp(allModels%*%log(matrix(Prob_diff,4))+(1-allModels)%*%log(1-matrix(Prob_diff,4)))
	x=barplot(tab,col='white',ylab='Model frequency',las=2)
	for (i in 2:5){
    	barplot(tab*allModels[,i-1],col=gsub('AA$','66',colERC5[i]),las=2,add=T,xaxt='n',yaxt='n')
	}
    x=barplot(tab,col='white',ylab='Model frequency',las=2)
	points(x,ExpectedPct*sum(tab),pch='-',col='grey',cex=3)
	for (i in 2:5){
    	barplot(tab*allModels[,i-1],col=gsub('AA$','66',colERC5[i]),las=2,add=T,xaxt='n',yaxt='n')
	}
    pbinom(as.numeric(tab),Nb_isomiR,as.numeric(ExpectedPct),low=F)
# 1.000000e+00 1.000000e+00 1.000000e+00 4.493919e-01 7.024789e-01 1.000000e+00 9.999934e-01 9.928960e-01 1.000000e+00 5.826827e-02 6.249472e-41 9.677411e-01 3.108674e-01 7.081362e-01 1.768179e-26
 dev.off()

################################
#  load isomiR ratio Data     ##
################################

#isomir_ratio_cast=fread(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_ratios_aggregated_nosubs.GCRL_Batch_lane_corrected.tsv",EVO_IMMUNO_POP))
#isomir_ratio_cast=isomir_ratio_cast[-grep('other',isomir_ID)]
#dim(isomir_ratio_cast)
#isomir_ratio_melted = melt(isomir_ratio_cast, id = c("mirID", "isomir_ID"),
#                               variable.name = "sample",
#                               value.name = "ratio")
#isomir_ratio_melted[, individual := substr(sample, 1,6)]
#isomir_ratio_melted[, condition := substr(sample, 8,8)]
#isomir_ratio_melted[, population := substr(sample, 1,3)]
#isomir_ratio_melted[,hsa_ID := gsub('(hsa.*)_(MIMAT.*)_(MI[0-9]+)','\\1',mirID)]

################################
#  load isomiR count Data     ##
################################

#isomir_count_cast=read.table(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_counts_aggregated_nosubs.log2RPM.GCRL_Batch_lane_corrected.tsv",EVO_IMMUNO_POP),check.names=FALSE)
#isomir_count_cast$mirID=isomir_ratio_cast$mirID
#isomir_count_cast$isomir_ID=isomir_ratio_cast$isomir_ID
#
#isomir_count_cast=as.data.table(isomir_count_cast)
#isomir_count_melted = melt(isomir_count_cast, id = c("mirID","isomir_ID"),
#                               variable.name = "sample",
#                               value.name = "count")
#isomir_count_melted[, individual := substr(sample, 1,6)]
#isomir_count_melted[, condition := substr(sample, 8,8)]
#isomir_count_melted[, population := substr(sample, 1,3)]
#isomir_count_melted[,hsa_ID := gsub('(hsa.*)_(MIMAT.*)_(MI[0-9]+)','\\1',mirID)]
#isomir_count_melted[,isomir_ID := V1]


##################################################
##  plot isomiR response 1 by 1  (function)     ##
##################################################

plot_isomiRNA_cond = function(miRNA,cond=1:5,add.violin=T,count=F){
	if(grepl('MIMAT',miRNA)){
		if(count){
       		temp=isomir_count_melted[mirID==miRNA & condition %in% cond,] 
	    }else{
    	    temp=isomir_ratio_melted[mirID==miRNA & condition %in% cond,]
        	temp[,ratio:=100*ratio]
        	}
        }else{
		if(count){
       		temp=isomir_count_melted[hsa_ID==miRNA & condition %in% cond,] 
	    }else{
    	    temp=isomir_ratio_melted[hsa_ID==miRNA & condition %in% cond,]
        	temp[,ratio:=100*ratio]
        	}
        }
    temp[,condition := factor(condIndex[as.numeric(condition)],levels=condIndex),]
    regexpr_isomiR='chr([0-9XYMT]+)_([-0-9]+)_([-0-9]+)_([+-])_([ATGC]+)_[0-2];(.*)_(hsa.*)_(MIMAT.*)_(MI[0-9]+)'
    regexpr_other='(hsa.*)_(MIMAT.*)_(MI[0-9]+)_other'
    temp[,isomir_ID_simple:=ifelse(!grepl('other',isomir_ID),gsub(regexpr_isomiR,'\\7_\\2;\\3;\\6',isomir_ID),gsub(regexpr_other,'\\1_other',isomir_ID))]
    if(count){
    p <- ggplot(temp,aes(x=isomir_ID_simple,y=count,fill=condition))+theme_bw()
    }else{
    p <- ggplot(temp,aes(x=isomir_ID_simple,y=ratio,fill=condition))+theme_bw()
    }
    p <- p + scale_fill_manual(values=colERC5[cond]) + geom_jitter(width=0.2,size=0.5,colour='darkgrey')+facet_grid(~condition)
  if(add.violin){
        p <- p + geom_violin(scale="width") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot(fill="#FFFFFF88", outlier.size=0, notch=TRUE,width=0.4) 
    }else{
        p <- p +  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot( outlier.size=0, notch=TRUE)
    }
      if(!count){
      p <- p + ylim(0,100) + xlab('') + ylab('Percentage of miRNA reads')
      }else{
      p <- p + ylab('log2(RPM)') + xlab('')
      }
    p <- p + theme(axis.text.x = element_text(angle = 60, hjust = 1))
    print(p)
}


####################################
##    plot isomiR response 1 by 1 ##
####################################

for (i in 1:length(isomir_DE_diff_miR)){
	cat(i)
	mymiR=isomir_DE_diff_miR[i]
    Nb_boxplots=length(unique(isomir_count_melted[mirID==miRNA,isomir_ID]))
	pdf(sprintf("%s/Maxime/miRNA_V2/figures/06.response_to_stimulation/exemples/isomirs/%s_%s-ALL_isomiRs.pdf",EVO_IMMUNO_POP,i,mymiR),height=4,width=4+.5*(Nb_boxplots-1))
	plot_isomiRNA_cond(miRNA=mymiR,cond=1:5)
	dev.off()
	pdf(sprintf("%s/Maxime/miRNA_V2/figures/06.response_to_stimulation/exemples/isomirs/%s_%s-ALL_isomiRs_count.pdf",EVO_IMMUNO_POP,i,mymiR),height=4,width=4+.5*(Nb_boxplots-1))
	plot_isomiRNA_cond(miRNA=mymiR,cond=1:5,count=T)
	dev.off()
}



#########################################################################
## 				 Number of Up/Down regulated Canonical isomiR		   ##
#########################################################################

to_plot = canonical_isomiR_DE[,.(signif=sum(fdr<.01),signifUp=sum(fdr<.01 & beta>0),signifDown=sum(fdr<.01 & beta<0),PctUp=sum(fdr<.01 & beta>0)/sum(fdr<.01)) ,by=condition]


# to_plot = canonical_isomirs_differential_expression[, .(number_up_regulated = sum(type == "upregulated"),
#             number_down_regulated = sum(type == "downregulated")), by = condition]

to_plot_mat=t(as.matrix(to_plot[,mget(c("signifUp","signifDown"))]))
colnames(to_plot_mat)=to_plot[,condition]

pdf(paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/figures/06.response_to_stimulation/number_canonical_isomiR_updownregulated_per_condition_global.pdf",sep=""),height=3.5,width=3)
par(mar=c(7,4,4,4))
barplot(to_plot_mat[,condition_name[-1]],beside=T,col=colERC[-(1:2)],las=2,space=c(0,.5),ylab='Number of miRNAs')
legend('top',fill=grey(c(.8,.3)),legend=c('canonical isomiR is up-regulated','canonical isomiR is down-regulated'),bty='n',inset=-.6,xpd=T)
dev.off()


#########################################################################
########    effect size of isomiR change across conditions	 ############
#########################################################################

isomiR_DE[fdr<0.01,mean(abs(beta)>0.05),by= condition_test]
#   condition_test         V1
#4:            LPS 0.05214724
#2:       PAM3CSK4 0.06590258
#1:           R848 0.09988519
#3:            IAV 0.11070111

colERC5=c("#525252AA", "#E31A1CAA", "#33A02CAA", "#1F78B4AA", "#6A3D9AAA")
condIndex=c("NS","LPS","PAM3CSK4","R848","IAV")
isomiR_DE[,condition_test:=factor(condition_test,condIndex[-1])]
pdf(sprintf("%s/Maxime/miRNA_V2/figures/06.response_to_stimulation/isomirs_differentially_expressed_in_conditions_effectSize.pdf",EVO_IMMUNO_POP))
p <- ggplot(isomiR_DE[fdr<0.01,],aes(y=100*abs(beta),x=condition_test,fill=condition_test))
    p <- p + geom_violin() + scale_fill_manual(values=colERC5[-1]) 
    p <- p + scale_x_discrete(name="condition") + scale_y_sqrt(name ="abs(beta)",breaks=c(0,1,5,10,20,30,40,50))
    p <- p + geom_boxplot(fill="#FFFFFF88", outlier.size=0.5, notch=TRUE,width=0.4)
    p <- p + theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank() , legend.position = "none")
    p <- p + theme(axis.text.x = element_text(angle = 60, hjust = 1))
    print(p)
dev.off()



#########################################################
######		Table miRNA modifs response 			#####
#########################################################

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

# Pct_modif=fread(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/Pct_modifs_perSample.txt",EVO_IMMUNO_POP))
# 
# samples_names=unlist(fread(sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/sample_names_977_highQuality.tsv",EVO_IMMUNO_POP)))
# 
# i='is_3p_reduced'
# Pct_3preduction=fread(sprintf('%s//Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/data/is_3p_reduced_allmiRNAs.txt',EVO_IMMUNO_POP))
# ID=sort(unique(isomiR_annot$mirID))
# condition = factor(condIndex[as.numeric(substr(samples_names, 8,8)),condIndex)
# Pct_3preduction_3P=apply(Pct_3preduction[grep('-3p',ID),],2,mean)
# Pct_3preduction_5P=apply(Pct_3preduction[grep('-5p',ID),],2,mean)
# boxplot(Pct_3preduction_3P~condition)
# 
# boxplot(Pct_3preduction_5P~condition)