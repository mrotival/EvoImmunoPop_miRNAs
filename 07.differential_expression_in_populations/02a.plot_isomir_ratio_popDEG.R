###################################################################################
##Aim of the script be able to plot a visual explanation of a isomir ratio change##
###################################################################################

##############
##Libariries##
##############
require(ggplot2)
luq = function(x){length(unique(x))}
colERC5 = c("#525252AA","#E31A1CAA","#33A02CAA","#1F78B4AA","#6A3D9AAA")
colERC = c("#969696", "#525252", "#FB9A99", "#E31A1C", "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4", "#CAB2D6", "#6A3D9A")
condIndex = c("NS","LPS","PAM3CSK4","R848","IAV" )

names(colERC5) = condIndex
names(colERC)=paste(rep(condIndex,e=2),'-',rep(c('EUB','AFB'),5))

###########################
##### get GO info       ###
###########################

source(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/scripts/00a.GOannotation/01.Load_and_format_databases.R", sep=""))
GO_informations = load_BHF_UCL_miRNA()


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
##@@@@@@                  LOAD ALL DATA (per sample)                        @@@@@@@##
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

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
isomir_ratio_melted[, condition := condIndex[as.numeric(substr(sample, 8,8))]]
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
isomir_count_melted[, condition := condIndex[as.numeric(substr(sample, 8,8))]]
isomir_count_melted[, population := substr(sample, 1,3)]
isomir_count_melted[, hsa_ID := gsub('.*(hsa.*)_(MIMAT.*)_(MI[0-9]+).*','\\1',V1)]
isomir_count_melted[, mirID := gsub('.*(hsa.*)_(MIMAT.*)_(MI[0-9]+).*','\\1_\\2_\\3',V1)]
isomir_count_melted[, isomir_ID := V1]


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
##@@@@@@                         LOAD RESULTS                               @@@@@@@##
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

########
##Data##
########
mir_popDE = fread(sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/miRNA_differentially_expressed_between_populations.tsv",EVO_IMMUNO_POP))
mir_popDE=mir_popDE[order(pvalue,-log2FC_EminusA)]
mir_popDE_diff= unique(mir_popDE[fdr<.01 ,miRNA],by='miRNA')
#mir_popDE_diff_up= unique(mir_popDE[fdr<.01 & log2FC_EminusA>0,miRNA],by='miRNA')
mir_popDE[,minP:=min(fdr),by=miRNA]
mir_popDE[,maxFC:=max(abs(log2FC_EminusA)),by=miRNA]
mir_popDE[,condition_test:=factor(condition_test,levels=condIndex,ordered=TRUE)]
luq(mir_popDE_diff) # 244


### ismomiR
isomir_annot=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_nosubs_FULL.tsv',EVO_IMMUNO_POP))

isomiR_popDE=fread(sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/isomirs_differentially_expressed_between_populations.tsv",EVO_IMMUNO_POP))
isomiR_popDE[, canonical := isomir_annot[match(isomiR_popDE$isomir, isomir_annot$ID), is_cannonical]]

isomiR_popDE[, pbonf := p.adjust(pvalue, method = "bonferroni")]
isomiR_popDE[, fdr := p.adjust(pvalue, method = "fdr")]
isomiR_popDE[, minPval_miRNA := min(pvalue),by=miRNA]
isomiR_popDE[, minPval_isomiR := min(pvalue),by=isomir]
isomiR_popDE[, maxBeta_miRNA := max(abs(log2FC_EminusA)),by=miRNA]
isomiR_popDE[, maxBeta_isomiR := max(abs(log2FC_EminusA)),by=isomir]

isomiR_popDE=isomiR_popDE[order(minPval_miRNA,pvalue,-log2FC_EminusA)]

isomiR_popDE[,hsa_ID:=gsub('(.*)_MIMAT.*','\\1',miRNA)]
isomiR_popDE_diff= unique(isomiR_popDE[fdr<.01 ,isomir],by='isomir')
luq(isomiR_popDE_diff) # 528

isomiR_popDE_diff_miR= unique(isomiR_popDE[fdr<.01 ,hsa_ID],by='hsa_ID')
luq(isomiR_popDE_diff_miR) # 188

luq(c(isomiR_popDE_diff_miR,mir_popDE_diff)) # 351

luq(intersect(isomiR_popDE_diff_miR,mir_popDE_diff)) # 81


isomiR_popDE_diff_HigherInEUB= unique(isomiR_popDE[fdr<.01 & log2FC_EminusA>0,isomir],by='isomir')
luq(isomiR_popDE_diff_HigherInEUB) # 233
isomiR_popDE_diff_HigherInEUB_miR= unique(isomiR_popDE[fdr<.01 & log2FC_EminusA>0,hsa_ID],by='hsa_ID')
luq(isomiR_popDE_diff_HigherInEUB_miR) # 149

canonical_isomiR_popDE = isomiR_popDE[which(canonical),]
canonical_isomiR_popDE_diff_miR= unique(canonical_isomiR_popDE[fdr<.01 ,hsa_ID],by='hsa_ID')
luq(canonical_isomiR_popDE_diff_miR) # 106

canonical_isomiR_popDE[, condition := factor(condition_test,levels=condIndex)]

diff_canonical_isomiR=canonical_isomiR_popDE[fdr<.01,isomir]
diff_canonical_miR=canonical_isomiR_popDE[fdr<.01,miRNA]
Nb_isomiR=luq(canonical_isomiR_popDE[,isomir])

isomiR_popDE_lik = fread(sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/isomirs_differentially_expressed_between_populations_Likelihoods.tsv",EVO_IMMUNO_POP), colClasses=c("character","character","character","numeric","numeric","numeric","integer","numeric"))
isomiR_popDE_lik[,bestModel:=condition_test]
isomiR_popDE_lik[,condition_test:=NULL]

# isomiR_popDE_lik = isomiR_popDE_lik[isomir%in%diff_canonical_isomiR,]
nullModel = isomiR_popDE_lik[bestModel=='00000',mget(c('isomir','miRNA','logLik'))][order(miRNA,isomir)]
bestModel = unique(isomiR_popDE_lik[order(miRNA,isomir,-ProbModel)],by='isomir')
bestModel = merge(bestModel,nullModel,suffixes=c('','.null'),by=c('miRNA','isomir'),all.x=TRUE)
bestModel[,logDiff:=logLik-logLik.null]
bestModel[,LR_pvalue:=pchisq(2*logDiff,1,low=F)]
bestModel[,LR_fdr:=p.adjust(LR_pvalue,'fdr')]
bestModel=bestModel[,diff_isomiR:=isomir%in%isomiR_popDE_diff]
bestModel=bestModel[,diff_miR:=miRNA%in%isomiR_popDE_diff_miR]
# isomiR_popDE_lik = isomiR_popDE_lik[order(isomir,-ProbModel)]
# isomiR_popDE_lik = isomiR_popDE_lik[df>3]
# best_model=unique(isomiR_popDE_lik,by= "isomir")


#######################################################################
######       table with all FC + best model 			           ####
#######################################################################

isomiR_popDE_cast=dcast(isomiR_popDE,isomir+miRNA+canonical+minPval_miRNA+minPval_isomiR+maxBeta_miRNA+maxBeta_isomiR~condition_test,value.var=c('MeanAFB','MeanEUB','log2FC_EminusA','pvalue','fdr'))
isomiR_popDE_cast[,hsa_ID:=gsub('.*(hsa.*)_(MIMAT.*)_(MI[0-9]+)_(MI[0-9]+).*','\\1', miRNA)]
# add likelihood informations for DE miRs (other miRs are NA)

isomiR_popDE_table = merge(isomiR_popDE_cast,bestModel[which(diff_isomiR),mget(c("isomir","miRNA","bestModel","ProbModel","LR_fdr","diff_isomiR","diff_miR"))],by=c('isomir','miRNA'),all.x=TRUE)
isomiR_popDE_table[,hsa_ID:=gsub('(.*)_MIMAT.*','\\1',hsa_ID)]
isomiR_popDE_table=merge(isomiR_popDE_table,isomir_annot,all.x=T,by.x=c('isomir','hsa_ID'),by.y=c('ID','hsa_ID'))

cols=c('hsa_ID','MIMAT','MI_ID','shift_5p','shift_3p','isomiR_subs','isomir_sequence',"canonical","minPval_miRNA","minPval_isomiR",
        "MeanAFB_NS","MeanEUB_NS", "pvalue_NS","fdr_NS",
        "MeanAFB_LPS","MeanEUB_LPS", "pvalue_LPS","fdr_LPS", 
			"MeanAFB_PAM3CSK4","MeanEUB_PAM3CSK4", "pvalue_PAM3CSK4",   "fdr_PAM3CSK4",
			"MeanAFB_R848","MeanEUB_R848","pvalue_R848",  "fdr_R848",
			 "MeanAFB_IAV","MeanEUB_IAV","pvalue_IAV", "fdr_IAV","bestModel","ProbModel","isomir")
isomiR_popDE_table=isomiR_popDE_table[order(minPval_miRNA,minPval_isomiR),mget(cols)]
isomiR_popDE_table[,isomir_sequence:=gsub('T','U',isomir_sequence)]
isomiR_popDE_table[,isomiR_subs:=gsub('T','U',isomiR_subs)]

# order by minP
isomiR_popDE_table[,minPval_isomiR:=NULL]
isomiR_popDE_table[,minPval_miRNA:=NULL]

# write
fwrite(isomiR_popDE_table,file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4B_isomir_popDE_V2.tsv",EVO_IMMUNO_POP),sep='\t')

isomiR_popDE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4B_isomir_popDE_V2.tsv",EVO_IMMUNO_POP),colClass=c(rep('character',7),'logical',rep('numeric',20),'character','numeric','character'))
#isomiR_popDE_cast[,hsa_ID:=gsub('.*(hsa.*)_(MIMAT.*)_(MI[0-9]+).*','\\1', miRNA)]
#merge( isomiR_popDE_cast[maxBeta_isomiR>0.1 & minPval_isomiR<0.0001,],GO_informations,by.x='hsa_ID',by.y='miRNA')
# miR-187-3p is a good candidate
# miR-187-3p is a good candidate. So is miR-146a-3p.

#############################################################################
####    create a plot with sharing of isomiR changes across conditions   ####
#############################################################################

Nb_isomiR=luq(isomiR_popDE_table$isomir)

allModels_STIM=cbind(NS=rep(0:1,16),
				LPS=rep(rep(0:1,e=16),1),
				PAM3CSK4=rep(rep(0:1,e=8),2),
				R848=rep(rep(0:1,e=4),4),
				IAV=rep(rep(0:1,e=2),8))
rownames(allModels_STIM)=apply(allModels_STIM,1,paste,collapse='')

tab_NbCond_isomiR_popDE=table(apply(allModels_STIM,1,sum)[isomiR_popDE_table$bestModel])
#   1   2   3   4   5 
#  12  35  52 107 322 

pieCol=c("#41B6C4","#A1DAB4","#FFFFB2","#FECC5C","#E31A1C")
par(mar=c(3,3,3,3))
tabpct=paste(round(100*tab_NbCond_isomiR_popDE/sum(tab_NbCond_isomiR_popDE), 1),'%')
pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/pie_NbCond_isomiR_popDE.pdf",EVO_IMMUNO_POP),width=4,height=4)
pie(tab_NbCond_isomiR_popDE,col=pieCol,init.angle=90,labels=tabpct) # 5.5 x 5.5 inches
pie(tab_NbCond_isomiR_popDE,col=pieCol,init.angle=90,labels=rep(' ',5)) # 4 x 4 inches
dev.off()


code_values_ordered = c("00000","10000","01000","00100","00010","00001",
                                "11000","10100","10010","10001",
                                "01100","01010","01001","00110","00101","00011",
                                "11100","11010","11001","10110","10101","10011",
                                "01110","01101","01011","00111",
                                "01111","10111","11011","11101","11110","11111")

# allModels_STIM=cbind(NS=rep(0,8),
# 				LPS=rep(rep(0:1,e=8),1),
# 				PAM3CSK4=rep(rep(0:1,e=4),2),
# 				R848=rep(rep(0:1,e=2),4),
# 				IAV=rep(rep(0:1,e=1),8))
# rownames(allModels_STIM)=apply(allModels_STIM,1,paste,collapse='')
# 
# 
# code_values_ordered = c("01000","00100","00010","00001","01100","01010","01001",
#                         "00110","00101","00011","01110","01101","01011","00111","01111")

# pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/miRNA_differentially_expressed_in_conditions_Likelihoods_part1.pdf",EVO_IMMUNO_POP),height=3,width=5)
# par(mar=c(4,7,1,1))
# Image(t(allModels_STIM[code_values_ordered,-1]))
# dev.off()

tab=table(isomiR_popDE_table[, bestModel])
tab=tab[code_values_ordered]
tab[is.na(tab)]=0
names(tab)=code_values_ordered

pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/isomir_differentially_expressed_between_populations_Likelihoods_part2.pdf",EVO_IMMUNO_POP),height=3,width=5)
par(mar=c(4,8,3,3))
	allModels=allModels_STIM[code_values_ordered,]
    Prob_diff=sapply(1:5,function(i){sum(tab[substr(names(tab),i,i)==1])/Nb_isomiR})
    ExpectedPct=exp(allModels%*%log(matrix(Prob_diff,5))+(1-allModels)%*%log(1-matrix(Prob_diff,5)))
	x=barplot(tab,col='white',ylab='Model frequency',las=2)
	for (i in 1:5){
    	barplot(tab*allModels[,i],col=gsub('AA$','66',colERC5[i]),las=2,add=T,xaxt='n',yaxt='n')
	}
    x=barplot(tab,col='white',ylab='Model frequency',las=2)
	points(x,ExpectedPct*sum(tab),pch='-',col='grey',cex=3)
	for (i in 1:5){
    	barplot(tab*allModels[,i],col=gsub('AA$','66',colERC5[i]),las=2,add=T,xaxt='n',yaxt='n')
	}
    pbinom(as.numeric(tab),Nb_isomiR,as.numeric(ExpectedPct),low=F)
# 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 9.880805e-01 9.999767e-01 9.990579e-01 9.998474e-01 9.999751e-01 7.152404e-01 5.307685e-01 9.978964e-01 9.954876e-01 9.749120e-01 4.736626e-14 4.353265e-17 1.058043e-09 1.421647e-04 7.157477e-10 0.000000e+00

 dev.off()

tab=table(isomiR_popDE_table[which(canonical), bestModel])
tab=tab[code_values_ordered]
tab[is.na(tab)]=0
names(tab)=code_values_ordered


tab_NbCond_can_isomiR_popDE=table(apply(allModels_STIM,1,sum)[isomiR_popDE_table[which(canonical), bestModel]])
# 1  2  3  4  5 
# 5  6 15 15 78

pieCol=c("#41B6C4","#A1DAB4","#FFFFB2","#FECC5C","#E31A1C")
par(mar=c(3,3,3,3))
tabpct=paste(round(100*tab_NbCond_can_isomiR_popDE/sum(tab_NbCond_can_isomiR_popDE), 1),'%')
pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/pie_NbCond_can_isomiR_popDE.pdf",EVO_IMMUNO_POP),width=4,height=4)
pie(tab_NbCond_can_isomiR_popDE,col=pieCol,init.angle=90,labels=tabpct) # 5.5 x 5.5 inches
pie(tab_NbCond_can_isomiR_popDE,col=pieCol,init.angle=90,labels=rep(' ',5)) # 4 x 4 inches
dev.off()


pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/canonical_isomir_differentially_expressed_between_populations_Likelihoods_part2.pdf",EVO_IMMUNO_POP),height=3,width=5)
par(mar=c(4,8,3,3))
	allModels=allModels_STIM[code_values_ordered,]
    Prob_diff=sapply(1:5,function(i){sum(tab[substr(names(tab),i,i)==1])/Nb_isomiR})
    ExpectedPct=exp(allModels%*%log(matrix(Prob_diff,5))+(1-allModels)%*%log(1-matrix(Prob_diff,5)))
	x=barplot(tab,col='white',ylab='Model frequency',las=2)
	for (i in 1:5){
    	barplot(tab*allModels[,i],col=gsub('AA$','66',colERC5[i]),las=2,add=T,xaxt='n',yaxt='n')
	}
    x=barplot(tab,col='white',ylab='Model frequency',las=2)
	points(x,ExpectedPct*sum(tab),pch='-',col='grey',cex=3)
	for (i in 1:5){
    	barplot(tab*allModels[,i],col=gsub('AA$','66',colERC5[i]),las=2,add=T,xaxt='n',yaxt='n')
	}
    pbinom(as.numeric(tab),Nb_isomiR,as.numeric(ExpectedPct),low=F)
# 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 9.817707e-01 9.284064e-01 9.886667e-01 9.204912e-01 9.284064e-01 9.378378e-01 9.846543e-01 9.529156e-01 9.887988e-01 9.469260e-01 2.023492e-02 2.177191e-02 1.890967e-01 1.963208e-03 2.017935e-01 1.041580e-04 1.170113e-04 2.187748e-02 2.353361e-02 1.373153e-04 5.623087e-15 5.623087e-15 6.876507e-05 2.394188e-07 2.951638e-07 0.000000e+00
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


plot_isomiRNA_pop = function(miRNA,cond=1:5,add.violin=T,count=F){
	if(grepl('MIMAT',miRNA)){
		if(count){
       		temp=isomir_count_melted[mirID==miRNA & condition %in% condIndex[cond],] 
	    }else{
    	    temp=isomir_ratio_melted[mirID==miRNA & condition %in% condIndex[cond],]
        	temp[,ratio:=100*ratio]
        	}
        }else{
		if(count){
       		temp=isomir_count_melted[hsa_ID==miRNA & condition %in% condIndex[cond],] 
	    }else{
    	    temp=isomir_ratio_melted[hsa_ID==miRNA & condition %in% condIndex[cond],]
        	temp[,ratio:=100*ratio]
        	}
        }
    temp[,condition := factor(condition,levels=condIndex[cond],ordered=TRUE),]
    temp[,cond_pop := factor(paste(condition,'-',population),levels=paste(rep(condIndex,e=2),'-',rep(c('AFB','EUB'),5)))]

    
    regexpr_isomiR='chr([0-9XYMT]+)_([-0-9]+)_([-0-9]+)_([+-])_([ATGC]+)_[0-2];(.*)_(hsa.*)_(MIMAT.*)_(MI[0-9]+)'
    regexpr_other='(hsa.*)_(MIMAT.*)_(MI[0-9]+)_other'
    temp[,isomir_ID_simple:=ifelse(!grepl('other',isomir_ID),gsub(regexpr_isomiR,'\\7_\\2;\\3;\\6',isomir_ID),gsub(regexpr_other,'\\1_other',isomir_ID))]
		if(count){
		p <- ggplot(temp,aes(x=paste(isomir_ID_simple,population),y=count,fill=cond_pop))+theme_bw()+facet_grid(~condition)
		}else{
		p <- ggplot(temp,aes(x=paste(isomir_ID_simple,population),y=ratio,fill=cond_pop))+theme_bw()+facet_grid(~condition)
		}
		p <- p + scale_fill_manual(values=colERC) + geom_jitter(width=0.2,size=0.5,colour='darkgrey')+facet_grid(~condition)
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

for (i in 1:length(isomiR_popDE_diff_miR)){
	cat(i)
	mymiR=isomiR_popDE_diff_miR[i]
    Nb_boxplots=length(unique(isomir_count_melted[mirID==mymiR,isomir_ID]))
	pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/exemples/isomirs/%s_%s-ALL_isomiRs.pdf",EVO_IMMUNO_POP,i,mymiR),height=4,width=8+1*(Nb_boxplots-1))
	plot_isomiRNA_pop(miRNA=mymiR,cond=1:5)
	dev.off()
	pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/exemples/isomirs/%s_%s-ALL_isomiRs_count.pdf",EVO_IMMUNO_POP,i,mymiR),height=4,width=8+1*(Nb_boxplots-1))
	plot_isomiRNA_pop(miRNA=mymiR,cond=1:5,count=T)
	dev.off()
}



#########################################################################
## 				 Number of Up/Down regulated Canonical isomiR		   ##
#########################################################################

Nb_popDE_isomiR=canonical_isomiR_popDE[,.(signif=sum(fdr<.01),HigherEUB=sum(fdr<.01 & log2FC_EminusA>0),HigherAFB=sum(fdr<.01 & log2FC_EminusA<0),PctHigherEUB=sum(fdr<.01 & log2FC_EminusA>0)/sum(fdr<.01),PctPctHigherAFB=sum(fdr<.01 & log2FC_EminusA<0)/sum(fdr<.01)) ,by=condition]
Nb_popDE_isomiR=Nb_popDE_isomiR[order(condition),P:=pbinom(HigherEUB, signif, .5, low=F)]

#   condition signif HigherEUB HigherAFB PctHigherEUB PctPctHigherAFB            P
# 1:        NS     41        34         7    0.8292683       0.1707317 2.436799e-06
# 2:       LPS     53        49         4    0.9245283       0.0754717 2.759792e-12
# 3:  PAM3CSK4     58        48        10    0.8275862       0.1724138 4.479701e-08
# 4:      R848     70        59        11    0.8428571       0.1571429 4.002412e-10
# 5:       IAV     71        58        13    0.8169014       0.1830986 6.738834e-09
# canoncial isomiR are higher in europeans, which suggest a population bias in the definition of miRNA isomiRs


to_plot = Nb_popDE_isomiR
to_plot_mat=t(as.matrix(to_plot[,mget(c("HigherEUB","HigherAFB"))]))
colnames(to_plot_mat)=to_plot[,condition]

pdf(paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/figures/07.differential_expression_in_populations/number_canonical_isomiR_higherEA_per_condition_global.pdf",sep=""),height=3.5,width=3)
par(mar=c(7,4,4,4))
barplot(to_plot_mat[,condIndex],beside=T,col=colERC,las=2,space=c(0,.5),ylab='Number of miRNAs')
legend('top',fill=grey(c(.8,.3)),legend=c('canonical isomiR is more expressed in Europeans','canonical isomiR is more expressed in Africans'),bty='n',inset=-.6,xpd=T)
dev.off()
    

#########################################################################
########    effect size of isomiR change across conditions	 ############
#########################################################################

isomiR_popDE[fdr<0.01,mean(abs(log2FC_EminusA)>0.05),by= condition_test][order(factor(condition_test,levels=condIndex))]
#  condition_test         V1
#1:             NS 0.02325581
#2:            LPS 0.05633803
#3:       PAM3CSK4 0.08300395
#4:           R848 0.05244755
#5:            IAV 0.03470032


isomiR_popDE[,condition_test:=factor(condition_test,levels=condIndex)]
pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/isomirs_differentially_expressed_between_populations_effectSize.pdf",EVO_IMMUNO_POP))
p <- ggplot(isomiR_popDE[fdr<0.01,],aes(y=100*abs(log2FC_EminusA),x=condition_test,fill=condition_test))
    p <- p + geom_violin() + scale_fill_manual(values=colERC5) 
    p <- p + scale_x_discrete(name="condition") + scale_y_sqrt(name ="abs(log2FC_EminusA)",breaks=c(0,1,5,10,20,30,40,50))
    p <- p + geom_boxplot(fill="#FFFFFF88", outlier.size=0.5, notch=TRUE,width=0.4)
    p <- p + theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank() , legend.position = "none")
    p <- p + theme(axis.text.x = element_text(angle = 60, hjust = 1))
    print(p)
dev.off()
# no significant differences between conditions

#########################################################
######		Table miRNA modifs response 			#####
#########################################################

miR_modif_DE_all=fread(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/tests/population_differences.txt",EVO_IMMUNO_POP))
condPopIndex=c("MeanAFB_NS", "MeanEUB_NS","MeanAFB_LPS","MeanEUB_LPS", "MeanAFB_PAM3CSK4","MeanEUB_PAM3CSK4","MeanAFB_R848", "MeanEUB_R848", "MeanAFB_IAV", "MeanEUB_IAV")
colnames(miR_modif_DE_all)[2:11]=condPopIndex
tested_diff=c("is_3p_cannonical","is_3p_extended","is_3p_reduced","is_3p_template_canonical","is_3p_template_extended","is_3p_template_reduced","is_5p_canonical","is_5p_extended","is_5p_reduced","nta_3p_u_miR","nta_3p_a_miR")
miR_modif_DE_list=list()
miR_modif_DE=melt(miR_modif_DE_all[type%in%tested_diff, ], measure.vars=list(paste('MeanAFB',condIndex,sep='_'), paste('MeanEUB',condIndex,sep='_'), paste('P',condIndex,sep='_')), value.name = c("MeanAFB","MeanEUB", "pvalue"),variable.name = "condition")
miR_modif_DE$condition=condIndex[as.numeric(miR_modif_DE$condition)]
miR_modif_DE$fdr=p.adjust(miR_modif_DE$pvalue)
#miR_modif_DE$Delta=(miR_modif_DE$mean-miR_modif_DE$Mean_NS)*100
#miR_modif_DE$Mean_NS=miR_modif_DE$Mean_NS*100

miR_modif_DE[fdr<0.01,][order(pvalue)]

cols=c("arm","type","MeanAFB_NS", "MeanEUB_NS","pvalue_NS","fdr_NS","MeanAFB_LPS", "MeanEUB_LPS","pvalue_LPS","fdr_LPS","MeanAFB_PAM3CSK4", "MeanEUB_PAM3CSK4","pvalue_PAM3CSK4","fdr_PAM3CSK4","MeanAFB_R848", "MeanEUB_R848","pvalue_R848","fdr_R848","MeanAFB_IAV","MeanEUB_IAV","pvalue_IAV","fdr_IAV")
miR_modif_DE_arm=dcast(miR_modif_DE,arm + type ~ condition, value.var = c("MeanAFB","MeanEUB",'pvalue',"fdr"))[,mget(cols)]
fwrite(miR_modif_DE_arm,file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4D_isomiR_modifs_population_differences_V2.txt",EVO_IMMUNO_POP) ,sep='\t')

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