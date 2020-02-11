
library(scales)
require(ggplot2)

condIndex=c("NS","LPS","PAM3CSK4","R848","IAV")
colERC5=c("#525252AA","#E31A1CAA","#33A02CAA","#1F78B4AA","#6A3D9AAA")
names(colERC5)=condIndex

################################################################
####        Number of associated miRNAs per gene            ####
################################################################
glmnetResults=fread(sprintf('%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_withoutTranscription_withoutSVA/gene_expression_glmnet_alpha0.5_allConds.tsv',EVO_IMMUNO_POP))
glmnetResults=glmnetResults[probs_selection>.8 & Q==9,]
Nb_miRNAs_by_cond=glmnetResults[,.(Nb_miRNA=sum(grepl('hsa',selected_variables))),by=.(gene,condition)]
Nb_miRNAs_by_cond[,condition:=factor(condIndex[condition], condIndex)]
Nb_miRNAs_by_cond[,.(Pct=mean(as.numeric(Nb_miRNA)>0)),by=.(condition)]
#   condition       Pct
#1:        NS 0.3171410
#2:       LPS 0.3379711
#3:  PAM3CSK4 0.2491652
#4:      R848 0.4457783
#5:       IAV 0.3985530
glmnetResults=fread('/Volumes/evo_immuno_pop/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_withoutTranscription_withoutSVA/gene_expression_glmnet_alpha0.5_allConds_allQ_ProbOver0.5.tsv')

max(Nb_miRNAs_by_cond$Nb_miRNA)
# [1] 6
Nb_miRNAs_by_cond[,Nb_miRNA:=as.character(Nb_miRNA)]

#pdf(sprintf('%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/variableQ_globalFDR/Number_of_miRNAs_perGene_sqrt_scale.pdf',EVO_IMMUNO_POP),width=3.5,height=2.7)
#p <- ggplot(Nb_miRNAs_by_cond,aes(x= as.character(Nb_miRNA),fill=condition)) 
#p <- p + geom_bar(position="dodge") + scale_fill_manual(values=colERC5) + scale_y_continuous(trans = sqrt_trans(),    breaks = c(0,400,1600,3600,6400),limits=c(0,10000))
##breaks=c(0,10,100,250,1000,2500,5000,10000)) + ylim(0,10000)
#p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
#p <- p + ylab('Number of genes') + xlab('Number of associated miRNAs') + theme(legend.position = c(.8, 0.59))
#print(p)
#dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/variableQ_globalFDR/Number_of_miRNAs_perGene.pdf',EVO_IMMUNO_POP),width=3.5,height=2.7)
Nb_miRNAs_by_cond[,Nb_miRNA_char:=factor(ifelse(Nb_miRNA>3,'4+',as.character(Nb_miRNA)),levels=c('0','1','2','3','4+'))]

p <- ggplot(Nb_miRNAs_by_cond,aes(x= Nb_miRNA_char,fill=condition)) 
p <- p + geom_bar(position="dodge") + scale_fill_manual(values=colERC5) + scale_y_continuous(limits=c(0,10000))
#breaks=c(0,10,100,250,1000,2500,5000,10000)) + ylim(0,10000)
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p <- p + ylab('Number of genes') + xlab('Number of associated miRNAs') + theme(legend.position = c(.7, 0.59))
print(p)
dev.off()




###############################################################################
####        Characterize miR-gene associations  (basal state)              ####
###############################################################################

##### load CAR scores
cond=1
CARSCORE_cond = fread(sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/variance_explained_by_miRNA_CARScore_withoutTranscription_withoutSVA/variance_explained_by_miRNA_CARScore_%s.tsv",EVO_IMMUNO_POP,cond))
CARSCORE_cond = CARSCORE_cond[!is.na(variance_explained_by_pop)]
#CARSCORE_cond[, eval(colnames_to_delete) := NULL]

CARSCORE_cond_melt=melt(CARSCORE_cond,measure.vars=list(miRNA=paste('miRNA',1:10,sep='_'),variance=paste('miRNA_variance',1:10,sep='_'),sign=paste('miRNA_carscore_sign',1:10,sep='_')))
CARSCORE_cond_melt=CARSCORE_cond_melt[!is.na(miRNA),]
CARSCORE_cond_melt=CARSCORE_cond_melt[,code_miRNA:=paste(gene, miRNA)]

##### load miRBS
miRNABS = fread(paste(EVO_IMMUNO_POP, "Martin/miRNA/11.miRNA_binding_sites/results/miRNA_binding_sites_on_protein_coding_3primeUTR_simplified.tsv", sep=""))
miRNABS[, code := paste(EnsemblGeneID, gsub("-", "_", miRNA))]
miRNABS=miRNABS[EnsemblGeneID%in%CARSCORE_cond$gene,]
Count_BS=miRNABS[,.(NbTarget=length(unique(EnsemblGeneID))),by='miRNA']

##### add miRBS informations
CARSCORE_cond_melt_withBS=merge(CARSCORE_cond_melt,unique(miRNABS,by='code'),by.x=c('code_miRNA'),by.y=c('code'),all.x=T,suffixes=c('','.toremove'))
CARSCORE_cond_melt_withBS[,miRNA.toremove:=NULL]
CARSCORE_cond_melt_withBS[,isTarget:=ifelse(!is.na(CARSCORE_cond_melt_withBS$transcripts),'yes','no')]

##### Numbers part 1
nrow(CARSCORE_cond_melt) # 6,009
mean(CARSCORE_cond_melt$sign==1)*100 # 57.24746
CARSCORE_cond_melt_withBS[sign==-1,mean(isTarget=='yes')] # 12%

tab=table(CARSCORE_cond_melt_withBS$sign, CARSCORE_cond_melt_withBS$isTarget)
tab
#       no  yes
#  -1 2252  317
#  1  2942  498

odds.ratio(tab[2:1,])
#            LowerCI       OR  UpperCI alpha          P
#odds ratio 1.021282 1.163595 1.326265  0.05 0.02240232

mirDBV6=fread('/Volumes/evo_immuno_pop/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miR_db_V6/miRDB_v6.0_prediction_result.txt')
mirDBV6=mirDBV6[grep('hsa-',V1),]
colnames(mirDBV6)=c('miRNA','Refseq','DownReg')

corresp_RefSeq=fread('/Volumes/evo_immuno_pop/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miR_db_V6/mart_export.txt')
colnames(corresp_RefSeq)=make.names(colnames(corresp_RefSeq))
corresp_RefSeq=corresp_RefSeq[RefSeq.match.transcript!='',]
corresp_RefSeq[,RefSeq.match.gene:=strsplit(RefSeq.match.transcript,'.')[[1]][1],by= RefSeq.match.transcript]
corresp_RefSeq[,RefSeq.match.gene:=strsplit(RefSeq.match.transcript,'\\.')[[1]][1],by= RefSeq.match.transcript]
corresp_RefSeq=corresp_RefSeq[,.(Gene.stable.ID,Gene.name,RefSeq.match.gene)]

mirDBV6=merge(corresp_RefSeq,mirDBV6,all.y=T,by.x='RefSeq.match.gene',by.y='Refseq')
mirDBV6[,code:=paste(Gene.stable.ID,gsub('-','_',miRNA))]
CARSCORE_cond_melt_withBS_miRDBV6=merge(CARSCORE_cond_melt,unique(mirDBV6,by='code'),by.x=c('code_miRNA'),by.y=c('code'),all.x=T,suffixes=c('','.toremove'))
CARSCORE_cond_melt_withBS_miRDBV6[,isTarget:=ifelse(!is.na(CARSCORE_cond_melt_withBS_miRDBV6$Gene.name),'yes','no')]
CARSCORE_cond_melt_withBS_miRDBV6[,miRNA.toremove:=NULL]
CARSCORE_cond_melt_withBS_miRDBV6[,inMirDBV6:=ifelse(gene%in%unique(mirDBV6$Gene.stable.ID),'yes','no')]

##### Numbers part 1
nrow(CARSCORE_cond_melt_withBS_miRDBV6[inMirDBV6=='yes',]) # 3746
mean(CARSCORE_cond_melt_withBS_miRDBV6[inMirDBV6=='yes',sign]==1)*100 # 55.41911
CARSCORE_cond_melt_withBS_miRDBV6[sign==-1 & inMirDBV6=='yes',mean(isTarget=='yes')] # 12%

tab=table(CARSCORE_cond_melt_withBS_miRDBV6[inMirDBV6=='yes',sign], CARSCORE_cond_melt_withBS_miRDBV6[inMirDBV6=='yes',isTarget])
tab
       no  yes
  -1 1627   43
  1  2014   62
  
odds.ratio(tab[2:1,])
#             LowerCI        OR  UpperCI alpha         P
# odds ratio 0.5646873 0.8585516 1.295013  0.05 0.4864181

###### number are unchanged if we use targets from an independant source (miR DB) 


###############################################################################
####        Characterize miR-transcription associations  (basal state)     ####
###############################################################################
cor_miRNA_intron_noSV_BS=fread(sprintf('%s/Maxime/miRNA_V2/data/13.correlation_transcription_miRNA/cor_miRNA_intron_MatrixEQTL_allCond_noSV_BS_FDR5.txt',EVO_IMMUNO_POP))
cor_miRNA_intron_noSV_BS[,isTarget:=ifelse(!is.na(per1),'yes','no')]


resCor_allCond=list()
Count_allCond=list()

for (cond in 0:5){
    cat(cond)
    if(cond==0){
        Count=unique(cor_miRNA_intron_noSV_BS[FDR<.05,],by='code')[,.(NbCor=length(gene),NbCorTarget=sum(isTarget=='yes'),NbCorPos=sum(r>0),NbCorNeg=sum(r<0),NbCorPosTarget=sum(r>0 & isTarget=='yes'),NbCorNegTarget=sum(r<0 & isTarget=='yes')),by='miRNA']
        }else{
        Count=unique(cor_miRNA_intron_noSV_BS[FDR<.05 & condition==cond, ],by='code')[,.(NbCor=length(gene),NbCorTarget=sum(isTarget=='yes'),NbCorPos=sum(r>0),NbCorNeg=sum(r<0),NbCorPosTarget=sum(r>0 & isTarget=='yes'),NbCorNegTarget=sum(r<0 & isTarget=='yes')),by='miRNA']
        }
    Count=merge(Count_BS,Count)
    Count[,cond:=cond]
    Count_allCond[[1+cond]]=Count

    Count[,NbGene:=nrow(CARSCORE_cond)]
    Count[,nonTarget:=NbGene-NbTarget]
    Count[,nonCor:=NbGene-NbCor]
    Count[,nonTargetCor:=NbCor-NbCorTarget]
    Count[,nonCorTarget:=NbTarget-NbCorTarget]
    Count[,nonTargetnonCor:=nonCor-nonCorTarget]
    Count[,nonCorPos:=NbGene-NbCorPos]
    Count[,nonTargetCorPos:=NbCorPos-NbCorPosTarget]
    Count[,nonCorPosTarget:=NbTarget-NbCorPosTarget]
    Count[,nonTargetnonCorPos:=nonCorPos-nonCorPosTarget]
    Count[,nonCorNeg:=NbGene-NbCorNeg]
    Count[,nonTargetCorNeg:=NbCorNeg-NbCorNegTarget]
    Count[,nonCorNegTarget:=NbTarget-NbCorNegTarget]
    Count[,nonTargetnonCorNeg:=nonCorNeg-nonCorNegTarget]

    resCor=list()
    for (i in 1:nrow(Count)){
        tab=matrix(c(Count[i,nonTargetnonCor],Count[i,nonTargetCor],Count[i,nonCorTarget],Count[i,NbCorTarget]),2)
        resCorTarget=data.table(odds.ratio(tab),miR=Count[i,miRNA],Type='Cor')

        tab=matrix(c(Count[i,nonTargetnonCorPos],Count[i,nonTargetCorPos],Count[i,nonCorPosTarget],Count[i,NbCorPosTarget]),2)
        resCorPosTarget=data.table(odds.ratio(tab),miR=Count[i,miRNA],Type='CorPos')

        tab=matrix(c(Count[i,nonTargetnonCorNeg],Count[i,nonTargetCorNeg],Count[i,nonCorNegTarget],Count[i,NbCorNegTarget]),2)
        resCorNegTarget=data.table(odds.ratio(tab),miR=Count[i,miRNA],Type='CorNeg')
        resCor[[i]]=rbind(resCorTarget,resCorPosTarget,resCorNegTarget)
    }
    resCor=rbindlist(resCor)
    resCor=resCor[Type!='Cor']
    resCor[,FDR:=p.adjust(P,'fdr')]
    resCor[,cond:=cond]
    resCor_allCond[[1+cond]]=resCor	
}
resCor_allCond=rbindlist(resCor_allCond)
Count_allCond=rbindlist(Count_allCond)

  cond    NbCor
1:    0 441.8389
2:    1 100.3070
3:    2 137.9492
4:    3 102.9483
5:    4 173.3458
6:    5 165.4196
Count_allCond[,mean(NbCor),by=cond]

fwrite(resCor_allCond,file=sprintf('%s/Maxime/miRNA_V2/data/13.correlation_transcription_miRNA/Enrichment_miRNABS_in_Genes_with_correlated_transcription_rate.txt',EVO_IMMUNO_POP),sep='\t')
fwrite(Count_allCond,file=sprintf('%s/Maxime/miRNA_V2/data/13.correlation_transcription_miRNA/Counts_miRNABS_in_Genes_with_correlated_transcription_rate.txt',EVO_IMMUNO_POP),sep='\t')

resCor_allCond=fread(sprintf('%s/Maxime/miRNA_V2/data/13.correlation_transcription_miRNA/Enrichment_miRNABS_in_Genes_with_correlated_transcription_rate.txt',EVO_IMMUNO_POP))
Count_allCond=fread(sprintf('%s/Maxime/miRNA_V2/data/13.correlation_transcription_miRNA/Counts_miRNABS_in_Genes_with_correlated_transcription_rate.txt',EVO_IMMUNO_POP))
resCor_allCond_cast=dcast(resCor_allCond, miR+cond~ Type,value.var=c('LowerCI','OR','UpperCI','P','FDR'))

Count_allCond[(FDR_CorNeg<0.05 | FDR_CorPos<0.05) & cond==0,]

rename=function(dt,old,new){
    colnames(dt)[colnames(dt)==old]=new
    dt
}
Count_allCond=rename(Count_allCond,'miRNA','miR')

Count_allCond=merge(Count_allCond,resCor_allCond_cast,by=c('miR','cond'))
Count_allCond_full = merge(resCor_allCond,Count_allCond,by=c('miR','cond'))

Data=Count_allCond_full[cond==0 & ( (OR_CorNeg>1 & FDR_CorNeg<0.05) | (OR_CorPos>1  & FDR_CorPos<0.05)),]

p <- ggplot(Data,aes(x=miR,y= OR,col= Type))+geom_point(position = position_dodge(width=.5))
p <- p + geom_errorbar( ymin = Data$LowerCI,ymax= Data$UpperCI,position = position_dodge(width=.5)) + ylim(c(-1,12)) + scale_color_manual(values=colERC5[c(4,2)])
p <- p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),axis.text.y = element_text(angle = 90, hjust = .5, vjust=0.5)) 
p <- p + geom_hline(aes(yintercept=1)) + xlab('')
#p <- p + geom_violin( scale = "width") + scale_fill_manual(values=rev(c("#E31A1C","#FFFFB2","#A1DAB4","#225EA8")))
#p <- p + geom_boxplot(width = 0.1, fill = '#FFFFFF88', outlier.shape = NA)
p <- p + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
#p <- p + theme(axis.text.x = element_text(angle = 90),axis.text.y = element_text(angle = 90)) + ylab('Percentage of variance explained')
#p <- p + coord_flip() + theme(legend.position="none") 
#p <- p + ylab('Percentage of variance explained') + xlab('')

Data=Count_allCond_full[cond==0 & ( (OR_CorNeg>1 & FDR_CorNeg<0.05) | (OR_CorPos>1  & FDR_CorPos<0.05)),]


p <- ggplot(Data,aes(x=miR,y= OR,col= Type))+geom_point(position = position_dodge(width=.5))

pdf(sprintf('%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/OR_miRNA_transcription_rate.pdf',EVO_IMMUNO_POP), width=5, height=5)
print(p)
dev.off()
p <- ggplot(Data,aes(x=miR,y= NbCor,col= Type))+geom_point(position = position_dodge(width=.5))


############################################################################################
####        Number of associated miRNAs per gene (ajusted on transcription)             ####
############################################################################################

#glmnetResults=fread(sprintf('%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/glmnet_correlation_miRNA_geneRPKM_transcriptionAdjusted_svaAdjusted/gene_expression_glmnet_alpha0.5_allConds.tsv',EVO_IMMUNO_POP))
glmnetResults=fread(sprintf('%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_transcriptionTested_withoutSVA/gene_expression_glmnet_alpha0.5_allConds.tsv',EVO_IMMUNO_POP))
glmnetResults=glmnetResults[probs_selection>.8 & Q==9,]

Nb_miRNAs_by_cond=glmnetResults[,.(Nb_miRNA=sum(grepl('hsa',selected_variables))),by=.(gene,condition)]
Nb_miRNAs_by_cond[,condition:=factor(condIndex[condition], condIndex)]
Nb_miRNAs_by_cond[,.(Pct=mean(as.numeric(Nb_miRNA)>0)),by=.(condition)]
# condition       Pct
#1:        NS 0.2401582
#2:       LPS 0.2204747
#3:  PAM3CSK4 0.1489923
#4:      R848 0.3111697
#5:       IAV 0.2959126
Nb_miRNAs_by_cond[,Nb_miRNA:=as.character(Nb_miRNA)]
###############################################################################
####        Characterize miR-gene associations  (basal state)              ####
###############################################################################

##### load CAR scores
cond=1
CARSCORE_cond = fread(sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/variance_explained_by_miRNA_CARScore_transcriptionTested_withoutSVA/variance_explained_by_miRNA_CARScore_%s.tsv",EVO_IMMUNO_POP,cond))
CARSCORE_cond = CARSCORE_cond[!is.na(variance_explained_by_pop)]
#CARSCORE_cond[, eval(colnames_to_delete) := NULL]

CARSCORE_cond_melt=melt(CARSCORE_cond,measure.vars=list(miRNA=paste('miRNA',1:10,sep='_'),variance=paste('miRNA_variance',1:10,sep='_'),sign=paste('miRNA_carscore_sign',1:10,sep='_')))
CARSCORE_cond_melt=CARSCORE_cond_melt[!is.na(miRNA),]
CARSCORE_cond_melt=CARSCORE_cond_melt[,code_miRNA:=paste(gene, miRNA)]


##### load miRBS

miRNABS = fread(paste(EVO_IMMUNO_POP, "Martin/miRNA/11.miRNA_binding_sites/results/miRNA_binding_sites_on_protein_coding_3primeUTR_simplified.tsv", sep=""))
miRNABS[, code := paste(EnsemblGeneID, gsub("-", "_", miRNA))]
miRNABS=miRNABS[EnsemblGeneID%in%CARSCORE_cond$gene,]
Count_BS=miRNABS[,.(NbTarget=length(unique(EnsemblGeneID))),by='miRNA']

##### add miRBS informations
CARSCORE_cond_melt_withBS=merge(CARSCORE_cond_melt,unique(miRNABS,by='code'),by.x=c('code_miRNA'),by.y=c('code'),all.x=T,suffixes=c('','.toremove'))
CARSCORE_cond_melt_withBS[,miRNA.toremove:=NULL]
CARSCORE_cond_melt_withBS[,isTarget:=ifelse(!is.na(CARSCORE_cond_melt_withBS$transcripts),'yes','no')]


##### Numbers part 1
#nrow(CARSCORE_cond_melt) # 3,252
#
#tab=table(CARSCORE_cond_melt_withBS$sign, CARSCORE_cond_melt_withBS$isTarget)
#tab



###############################################################################
####        Characterize miR-gene associations  (basal state)              ####
###############################################################################

##### load miRBS
miRNABS = fread(paste(EVO_IMMUNO_POP, "Martin/miRNA/11.miRNA_binding_sites/results/miRNA_binding_sites_on_protein_coding_3primeUTR_simplified.tsv", sep=""))
miRNABS[, code := paste(EnsemblGeneID, gsub("-", "_", miRNA))]
# 	miRNABS=miRNABS[EnsemblGeneID%in%CARSCORE_cond$gene,]
# 	Count_BS=miRNABS[,.(NbTarget=length(unique(EnsemblGeneID))),by='miRNA']

######### load CAR scores
#CARSCORE_noTranscription=list()
#CARSCORE_melt_noTranscription=list()
#for (cond in 1:5){
#	CARSCORE_cond = fread(sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/old_with_bad_SVA/variance_explained_by_miRNA_CARScore_without_transcription/variance_explained_by_miRNA_CARScore_%s.tsv",EVO_IMMUNO_POP,cond))
#	CARSCORE_cond = CARSCORE_cond[!is.na(variance_explained_by_pop)]
#	CARSCORE_noTranscription[[cond]]=CARSCORE_cond
#
#	CARSCORE_cond_melt=melt(CARSCORE_cond,measure.vars=list(miRNA=paste('miRNA',1:10,sep='_'),variance=paste('miRNA_variance',1:10,sep='_'),sign=paste('miRNA_carscore_sign',1:10,sep='_')))
#	CARSCORE_cond_melt=CARSCORE_cond_melt[!is.na(miRNA),]
#	CARSCORE_cond_melt=CARSCORE_cond_melt[,code_miRNA:=paste(gene, miRNA)]
#
#
##########	add miRBS informations
#	CARSCORE_cond_melt_withBS=merge(CARSCORE_cond_melt,unique(miRNABS,by='code'),by.x=c('code_miRNA'),by.y=c('code'),all.x=T,suffixes=c('','.toremove'))
#	CARSCORE_cond_melt_withBS[,miRNA.toremove:=NULL]
#	CARSCORE_cond_melt_withBS[,isTarget:=ifelse(!is.na(CARSCORE_cond_melt_withBS$transcripts),'yes','no')]
#	CARSCORE_melt_noTranscription[[cond]]=CARSCORE_cond_melt_withBS
#}
#CARSCORE_noTranscription = rbindlist(CARSCORE_noTranscription)
#CARSCORE_melt_noTranscription = rbindlist(CARSCORE_melt_noTranscription)
#
#
#info_to_add=CARSCORE_melt_noTranscription[, .(variance_explained_by_miRNAs_neg=sum(ifelse(sign==-1,variance,0)), variance_explained_by_miRNAs_binding_neg=sum(ifelse(isTarget=='yes' & sign==-1,variance,0))),by=.(gene,condition)]
#CARSCORE_noTranscription=merge(CARSCORE_noTranscription,info_to_add,all.x=T)
#CARSCORE_noTranscription[is.na(variance_explained_by_miRNAs_binding_neg),]$variance_explained_by_miRNAs_binding_neg=0
#CARSCORE_noTranscription[is.na(variance_explained_by_miRNAs_neg),]$variance_explained_by_miRNAs_neg=0
#
#
#p <- ggplot(to_plot, aes(x = type_of_variance, y = variance_explained, fill = type_of_variance))
#p <- p + geom_violin( scale = "width")
#p <- p + geom_boxplot(width = 0.1, fill = NA, outlier.shape = NA)


##### load CAR scores
CARSCORE=list()
CARSCORE_melt=list()
for (cond in 1:5){
	CARSCORE_cond = fread(sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/variance_explained_by_miRNA_CARScore_transcriptionTested_withoutSVA/variance_explained_by_miRNA_CARScore_%s.tsv",EVO_IMMUNO_POP,cond))
	CARSCORE_cond = CARSCORE_cond[!is.na(variance_explained_by_pop)]
	CARSCORE[[cond]]=CARSCORE_cond

	CARSCORE_cond_melt=melt(CARSCORE_cond,measure.vars=list(miRNA=paste('miRNA',1:10,sep='_'),variance=paste('miRNA_variance',1:10,sep='_'),sign=paste('miRNA_carscore_sign',1:10,sep='_')))
	CARSCORE_cond_melt=CARSCORE_cond_melt[!is.na(miRNA),]
	CARSCORE_cond_melt=CARSCORE_cond_melt[,code_miRNA:=paste(gene, miRNA)]

	##### add miRBS informations
	CARSCORE_cond_melt_withBS=merge(CARSCORE_cond_melt,unique(miRNABS,by='code'),by.x=c('code_miRNA'),by.y=c('code'),all.x=T,suffixes=c('','.toremove'))
	CARSCORE_cond_melt_withBS[,miRNA.toremove:=NULL]
	CARSCORE_cond_melt_withBS[,isTarget:=ifelse(!is.na(CARSCORE_cond_melt_withBS$transcripts),'yes','no')]
	CARSCORE_melt[[cond]]=CARSCORE_cond_melt_withBS
}
CARSCORE = rbindlist(CARSCORE)
CARSCORE_melt = rbindlist(CARSCORE_melt)
CARSCORE_melt=CARSCORE_melt[!is.na(variance),]

info_to_add=CARSCORE_melt[, .(variance_explained_by_miRNAs_neg=sum(ifelse(sign==-1,variance,0)),variance_explained_by_miRNAs_binding_neg=sum(ifelse(isTarget=='yes' & sign==-1,variance,0))),by=.(gene,condition)]

CARSCORE=merge(CARSCORE,info_to_add,all.x=T)
CARSCORE[is.na(variance_explained_by_miRNAs_binding_neg),]$variance_explained_by_miRNAs_binding_neg=0
CARSCORE[is.na(variance_explained_by_miRNAs_neg),]$variance_explained_by_miRNAs_neg=0

CARSCORE[!paste(gene,cond)%in%glmnetResults[selected_variables=='gene_transcription',paste(gene,cond)],variance_explained_by_transcription:=0]

to_plot = CARSCORE[condition==1, .(gene, variance_explained_by_transcription, variance_explained_by_miRNAs, variance_explained_by_miRNAs_neg, variance_explained_by_miRNAs_binding_neg)]
to_plot = melt(to_plot, id = "gene", variable.name = "type_of_variance", value.name = "variance_explained")

to_plot[type_of_variance == "variance_explained_by_transcription", type_of_variance := "Transcription"]
to_plot[type_of_variance == "variance_explained_by_miRNAs", type_of_variance := "miRNAs"]
to_plot[type_of_variance == "variance_explained_by_miRNAs_neg", type_of_variance := "Downregulating miRNAs"]
to_plot[type_of_variance == "variance_explained_by_miRNAs_binding_neg", type_of_variance := "Downregulating miRNAs\n with predicted binding"]
to_plot[, type_of_variance := factor(type_of_variance, levels = rev(c("Transcription", "miRNAs", "Downregulating miRNAs", "Downregulating miRNAs\n with predicted binding")))]
to_plot[,variance_explained:=variance_explained*100]
to_plot[,contributing:=ifelse(variance_explained>0,"yes","no")]

# basal car score informations
pdf(sprintf('%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/CARscores_basal_transcription_tested.pdf',EVO_IMMUNO_POP),width=4.5,height=2.7)
p <- ggplot(to_plot[variance_explained>0], aes(x = type_of_variance, y = variance_explained, fill = type_of_variance))
p <- p + geom_violin( scale = "width") + scale_fill_manual(values=rev(c("#E31A1C","#FFFFB2","#A1DAB4","#225EA8")))
p <- p + geom_boxplot(width = 0.1, fill = '#FFFFFF88', outlier.shape = NA)
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
#p <- p + theme(axis.text.x = element_text(angle = 90),axis.text.y = element_text(angle = 90)) + ylab('Percentage of variance explained')
p <- p + coord_flip() + theme(legend.position="none") + ylim(0,100)
p <- p + ylab('Percentage of variance explained') + xlab('')
print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/CARscores_basal_pies_transcription_tested.pdf',EVO_IMMUNO_POP),width=4.5,height=2.7)
pie_data=to_plot[, .(regulated= mean(variance_explained>0)*100, non_regulated= mean(variance_explained==0)*100),by=type_of_variance]
pie_data
#                                type_of_variance regulated non_regulated  x y
#1:                                  Transcription  86.19151      13.80849 95 4
#2:                                         miRNAs  25.60468      74.39532 95 3
#3:                          Downregulating miRNAs  13.49128      86.50872 95 2
#4: Downregulating miRNAs\n with predicted binding   1.93299      98.06701 95 1
pie_data[,x:=95]
pie_data[,y:=4:1]
par(mar=c(2.1,0,0.3,0.2))
plot(c(0,100),y=.5+c(0,4),axes=F,xlab='',ylab='',col='#00000000')
library(plotrix)
for (i in 1:4){
    floating.pie(pie_data[i,x],pie_data[i,y],c(pie_data[i,regulated],pie_data[i,non_regulated]),r=3,col=c(c("#E31A1C","#FFFFB2","#A1DAB4","#225EA8")[i],"lightgrey"),startpos=pi/2-pie_data[i,regulated]*2*pi/100,lwd=0.5)
}
par(mar=c(5.1,4.1,4.1,2.1))
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/CARscores_basal_bars_transcription_tested.pdf',EVO_IMMUNO_POP),width=4.5,height=2.7)
bar_data=to_plot[, .(regulated = mean(variance_explained)) ,by=type_of_variance]
p <- ggplot(bar_data,aes(x=type_of_variance, y=regulated, fill = type_of_variance))
p <- p + geom_bar(stat="identity") + scale_fill_manual(values=rev(c("#E31A1C","#FFFFB2","#A1DAB4","#225EA8")))
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p <- p + coord_flip() + theme(legend.position="none") + ylim(0,30) + ylab('Percentage of variance explained')+xlab('')
p <- p + geom_text(aes(y=regulated, label=paste(round(regulated,1),'%')), vjust=0.5,hjust=-.1, size=3.5)
print(p)
dev.off()

#Transcriptions
pdf(sprintf('%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/CARscores_transcription_transcription_tested.pdf', EVO_IMMUNO_POP),width=4.5,height=2.7)
to_plot <- CARSCORE[1:.N,]
to_plot[,condition:=factor(condIndex[condition],levels=condIndex),]
to_plot[,variance_explained_by_transcription:=variance_explained_by_transcription*100]
p <- ggplot(to_plot[variance_explained_by_transcription>0,], aes(x = condition, y = variance_explained_by_transcription, fill = condition))
p <- p + geom_violin(scale = "width")
p <- p + geom_boxplot(width = 0.2, fill = '#FFFFFF88', notch = T,  outlier.shape = NA)
p <- p + scale_fill_manual(values = colERC5)
p <- p + ylab('Percentage of variance explained') + xlab('') + ylim(0,100)
p <- p + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/CARscores_transcription_pies_transcription_tested.pdf',EVO_IMMUNO_POP),width=4.5,height=2.7)
pie_data=to_plot[, .(regulated= mean(variance_explained_by_transcription>0)*100, non_regulated= mean(variance_explained_by_transcription==0)*100),by=condition]
pie_data[,x:=95]
pie_data[,y:=5:1]
par(mar=c(2.1,0,0.3,0.2))
plot(c(0,100),y=.5+c(0,5),axes=F,xlab='',ylab='',col='#00000000')
library(plotrix)
for (i in 1:5){
    floating.pie(pie_data[i,x],pie_data[i,y],c(pie_data[i,regulated],pie_data[i,non_regulated]),r=3,col=c(colERC5[i],"lightgrey"),startpos=pi/2-pie_data[i,regulated]*2*pi/100,lwd=0.5)
}
par(mar=c(5.1,4.1,4.1,2.1))
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/CARscores_transcription_bars_transcription_tested.pdf',EVO_IMMUNO_POP),width=3.2,height=4.5)
bar_data=to_plot[, .(regulated = mean(variance_explained_by_transcription)) ,by=condition]
p <- ggplot(bar_data,aes(x=condition, y=regulated, fill = condition))
p <- p + geom_bar(stat="identity") + scale_fill_manual(values=colERC5)
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p <- p + theme(legend.position="none") + ylim(0,30) + ylab('Percentage of variance explained')+xlab('')
p <- p + geom_text(aes(y=regulated, label=paste(round(regulated,1),'%')), vjust=-0.2,hjust=0.5, size=3.5)
p <- p + theme(axis.text.x = element_text(angle = 45, hjust=1)) 
print(p)
dev.off()

#miRNAs
pdf(sprintf('%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/CARscores_miRNAs_transcription_tested.pdf', EVO_IMMUNO_POP),width=4.5,height=2.7)
to_plot <- CARSCORE[1:.N,]
to_plot[,condition:=factor(condIndex[condition],levels=condIndex),]
to_plot[,variance_explained_by_miRNAs:=variance_explained_by_miRNAs*100]
p <- ggplot(to_plot[variance_explained_by_miRNAs>0,], aes(x = condition, y = variance_explained_by_miRNAs, fill = condition))
p <- p + geom_violin(scale = "width")
p <- p + geom_boxplot(width = 0.2, fill = '#FFFFFF88', notch = T,  outlier.shape = NA)
p <- p + scale_fill_manual(values = colERC5)
p <- p + ylab('Percentage of variance explained') + xlab('') + ylim(0,100)
p <- p + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/CARscores_miRNAs_pies_transcription_tested.pdf',EVO_IMMUNO_POP),width=4.5,height=2.7)
pie_data=to_plot[, .(regulated= mean(variance_explained_by_miRNAs>0)*100, non_regulated= mean(variance_explained_by_miRNAs==0)*100),by=condition]
pie_data[,x:=95]
pie_data[,y:=5:1]
par(mar=c(2.1,0,0.3,0.2))
plot(c(0,100),y=.5+c(0,5),axes=F,xlab='',ylab='',col='#00000000')
library(plotrix)
for (i in 1:5){
    floating.pie(pie_data[i,x],pie_data[i,y],c(pie_data[i,regulated],pie_data[i,non_regulated]),r=3,col=c(colERC5[i],"lightgrey"),startpos=pi/2-pie_data[i,regulated]*2*pi/100,lwd=0.5)
}
par(mar=c(5.1,4.1,4.1,2.1))
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/CARscores_miRNAs_bars_transcription_tested.pdf',EVO_IMMUNO_POP),width=3.2,height=4.5)
bar_data=to_plot[, .(regulated = mean(variance_explained_by_miRNAs)) ,by=condition]
p <- ggplot(bar_data,aes(x=condition, y=regulated, fill = condition))
p <- p + geom_bar(stat="identity") + scale_fill_manual(values=colERC5)
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p <- p + theme(legend.position="none") + ylim(0,30) + ylab('Percentage of variance explained')+xlab('')
p <- p + geom_text(aes(y=regulated, label=paste(round(regulated,1),'%')), vjust=-0.2,hjust=0.5, size=3.5)
p <- p + theme(axis.text.x = element_text(angle = 45, hjust=1)) 
print(p)
dev.off()


#p <- ggplot(CARSCORES, aes(x = condition, y = variance_explained_by_miRNAs, fill = condition))
#p <- p + geom_violin(scale = "width")
#p <- p + geom_boxplot(width = 0.1, outlier.shape = NA, notch = T)
#p <- p + scale_fill_manual(values = colERC5)
#p <- p + ylab("expression variation explained by miRNAs")
#p <- p + guides(fill=FALSE)
#p <- p + theme_bw()
#
#pdf(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/", "expression_variation_explaind_by_miRNAs.pdf", sep=""), width = 4, height =4)
#print(p)
#dev.off()
#

##### Numbers part 1
sapply(2:5,function(i){wilcox.test(CARSCORE[condition%in%c(1,i),variance_explained_by_transcription]~CARSCORE[condition%in%c(1,i),condition])$p.value})
#[1] 2.874866e-03 2.511714e-31 1.358159e-23 1.079831e-38
sapply(2:5,function(i){wilcox.test(CARSCORE[condition%in%c(1,i),variance_explained_by_miRNAs]~CARSCORE[condition%in%c(1,i),condition])$p.value})
#[1]  1.234967e-03  4.614175e-16 1.060841e-162  6.070905e-54
