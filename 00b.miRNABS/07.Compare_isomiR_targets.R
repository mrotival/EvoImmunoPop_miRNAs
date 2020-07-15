
DATA_DIR=sprintf("%s/Maxime/miRNA_V2/data/",EVO_IMMUNO_POP)
FIG_DIR=sprintf("%s/Maxime/miRNA_V2/figures/Revisions",EVO_IMMUNO_POP)

###########################################################
###				miRANDA  targets 						###
###########################################################

target_miRANDA=fread(sprintf('%s/Martin/miRNA/11.miRNA_binding_sites/results/miRNA_binding_sites_on_protein_coding_3primeUTR_simplified.tsv',EVO_IMMUNO_POP))
allGenes=unique(target_miRANDA$EnsemblGeneID)
target_miRANDA[,NbTarget:=max(table(transcripts)),by=.(miRNA,EnsemblGeneID,GeneNames)]

target_miRANDA[,per1:=as.numeric(substr(per1,1,3))/100]
target_miRANDA[,per2:=as.numeric(substr(per2,1,3))/100]
target_miRANDA[,score:=miRNABS_length*per1*per2/max(miRNABS_length)]
target_miRANDA[,score:=rank(score)/(.N+1)]
target_miRANDA[,target_score:=max(score),by=.(transcripts,miRNA,EnsemblGeneID,GeneNames)]
target_miRANDA_gene=target_miRANDA[order(miRNA,EnsemblGeneID,GeneNames,-NbTarget,-target_score),][which(!duplicated(paste(miRNA,EnsemblGeneID,GeneNames))),]
target_miRANDA_gene[,score:=NULL]
fwrite(target_miRANDA_gene[,mget(c('miRNA','EnsemblGeneID','GeneNames','target_score'))],file=sprintf('%s/19_miR_targets/target_miRANDA_gene_with_HUGO.txt',DATA_DIR),sep='\t')




target_isomiRs=fread(sprintf('%s/Martin/miRNA/11.miRNA_binding_sites/results/common_isomiR_binding_sites_on_protein_coding_3primeUTR_simplified.tsv',EVO_IMMUNO_POP))
target_isomiRs[,per1:=as.numeric(substr(per1,1,3))/100]
target_isomiRs[,per2:=as.numeric(substr(per2,1,3))/100]
target_isomiRs[,score:=isomiRNABS_length*per1*per2/max(isomiRNABS_length)]
target_isomiRs[,score:=rank(score)/(.N+1)]
target_isomiRs[,target_score:=max(score),by=.(transcripts,isomirID,EnsemblGeneID,GeneNames)]
target_isomiRs[,NbTarget:=max(table(transcripts)),by=.(isomirID,EnsemblGeneID,GeneNames)]
target_isomiRs[,miRNA:=gsub('.*(hsa-.*)','\\1',isomirID)]
target_isomiRs_gene=target_isomiRs[order(isomirID,EnsemblGeneID,GeneNames,-NbTarget,-target_score),][which(!duplicated(paste(isomirID,EnsemblGeneID,GeneNames))),]


fwrite(target_isomiRs_gene,file=sprintf('%s/19_miR_targets/target_isomiR_genes.txt',DATA_DIR),sep='\t')



target_miRANDA_gene=fread(sprintf('%s/19_miR_targets/target_miRANDA_gene_with_HUGO.txt',DATA_DIR),sep='\t')
target_isomiRs_gene=fread(sprintf('%s/19_miR_targets/target_isomiR_genes.txt',DATA_DIR),sep='\t')
allGenes=unique(target_miRANDA_gene$EnsemblGeneID)
	
isomiR_Annot=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_commonOnly_V2.0.tsv',EVO_IMMUNO_POP))
isomiR_Annot[,nsubs_other:=nsubs-nta_3p_a_miR-nta_3p_u_miR]


# common non-canonical isomiRs that differ from the canonical only in their start/end site and NTA
isomiR_ID=isomiR_Annot[nsubs_other==0 & isomiR_type!='CAN',][order(-mean_isoMiR_NS),ID]

# library(seqinr)
# UTRseq=read.fasta(sprintf('%s/Martin/miRNA/11.miRNA_binding_sites/temp_data/temp_transcript.fasta',EVO_IMMUNO_POP))
# mean(sapply(UTRseq,length)) # 917.9 bp

Transcript_Annot=fread(sprintf('%s/Maxime/Evo_Immuno_pop_data/TranscriptAnnotation_hg37_ens70.txt',EVO_IMMUNO_POP))
UTR_length=Transcript_Annot[,.(UTR3_length =max(UTR3_length)),by=Ensembl.Gene.ID]
UTR3_length=UTR_length[match(allGenes,UTR_length$Ensembl.Gene.ID),UTR3_length]
names(UTR3_length)=allGenes
Gene_annot=as.data.frame(fread(sprintf('%s/Maxime/Evo_Immuno_pop_data/GeneAnnotation_hg37_ens70.txt',EVO_IMMUNO_POP)))
G2S=function(x){Gene_annot[match(x,Gene_annot[,1]),'Associated Gene Name']}
#allGOterms=as.data.frame(fread(sprintf('%s/Maxime/Shared/allGOterms_EnsGRC37_13042017.txt',EVO_IMMUNO_POP)))
# allGOterms=as.data.frame(fread(sprintf('%s/Maxime/miRNA_V2/data/19_miR_targets/allGO_terms_070520.txt',EVO_IMMUNO_POP)))
# colnames(allGOterms)[c(1,3)]=c('gene','go')

library(goseq)

gene2cat = goseq:::getgo(allGenes, 'hg19', 'ensGene', fetch.cats = c("GO:CC", "GO:BP", "GO:MF"))
names(gene2cat) = rownames(pwf)
cat2gene = goseq:::reversemapping(gene2cat)
gene2cat = goseq:::reversemapping(cat2gene)
allGOterms=data.table(go=rep(names(cat2gene),sapply(cat2gene,length)),gene=unlist(cat2gene))


i=0
GO_lost=list()
GO_new=list()
numbers=list()
for (myID in isomiR_ID){
 	i=i+1
	cat(i,"")
	if(i%%100==0){	cat(i,'/',length(isomiR_ID),'\n')}
	miR_ID=isomiR_Annot[ID==myID,hsa_ID]
	miRANDA_targets=target_miRANDA_gene[miRNA== miR_ID, EnsemblGeneID]
	isomiR_targets=target_isomiRs_gene[isomirID==myID,EnsemblGeneID]
	common_targets=intersect(miRANDA_targets, isomiR_targets)
	new_targets=setdiff(miRANDA_targets, isomiR_targets)
	lost_targets=setdiff(isomiR_targets,miRANDA_targets)
	numbers[[myID]]=c(length(common_targets),length(new_targets),length(lost_targets))
	if(myID%in%unique(GOres$isomiR)){
	if(length(new_targets)>30){
		GO_new[[myID]]=GOSeq(new_targets,allGenes,bias.data=UTR3_length)
		GO_new[[myID]]$isomiR=rep(myID,nrow(GO_new[[myID]]))
		}
	if(length(lost_targets)>30){
		GO_lost[[myID]]=GOSeq(lost_targets,allGenes,bias.data=UTR3_length)
		GO_lost[[myID]]$isomiR=rep(myID,nrow(GO_lost[[myID]]))
		}
	}
}
GO_new=rbindlist(GO_new)
GO_lost=rbindlist(GO_lost)

# GO_lost=fread(sprintf('%s/19_miR_targets/target_isomiR_GO_lost_3UTR.txt',DATA_DIR))
# GO_new=fread(sprintf('%s/19_miR_targets/target_isomiR_GO_new_3UTR.txt',DATA_DIR))
GO_new[,type:='gained']
GO_lost[,type:='lost']
GOres=rbind(GO_new, GO_lost)
# GOres[,MIMAT:=gsub('.*(MIMAT[0-9]+).*','\\1',isomiR)]
cols=c('ID',"hsa_ID","MIMAT","MI_ID",'shift_3p','shift_5p','isomiR_subs','isomir_sequence')
GOres=merge(GOres, isomiR_Annot[,mget(cols)],by.x='isomiR',by.y='ID')
cols=c("hsa_ID", "MI_ID","MIMAT", "shift_3p", "shift_5p", "isomiR_subs", "isomir_sequence", 
 "type", "category", "ontology", "Term", "numDEInCat", "numInCat","Pvalue","FDR", "FoldEnrich",  "genes" )

fwrite(GOres[,mget(cols)],file=sprintf('%s/19_miR_targets/target_isomiR_GO_both_3UTR_clean.txt',DATA_DIR),sep='\t')



GOres=fread(sprintf('%s/19_miR_targets/target_isomiR_GO_both_3UTR.txt',DATA_DIR))
numbers_DT=as.data.table(do.call(rbind,numbers))
colnames(numbers_DT)=c('common','new','lost')
numbers_DT$ID=isomiR_ID
cols=c('ID',"hsa_ID","MIMAT","MI_ID",'shift_3p','shift_5p','isomiR_subs','isomir_sequence')
numbers_DT=merge(isomiR_Annot[,mget(cols)],numbers_DT,by='ID')
colnames(numbers_DT)[10]='gained'
fwrite(numbers_DT[,-1],file=sprintf('%s/19_miR_targets/target_isomiR_number_different_3UTR_clean.txt',DATA_DIR),sep='\t')
# fwrite(rbindlist(GO_new),file=sprintf('%s/19_miR_targets/target_isomiR_GO_new_3UTR.txt',DATA_DIR),sep='\t')
# fwrite(rbindlist(GO_lost),file=sprintf('%s/19_miR_targets/target_isomiR_GO_lost_3UTR.txt',DATA_DIR),sep='\t')
numbers_DT=fread(sprintf('%s/19_miR_targets/target_isomiR_number_different_3UTR_clean.txt',DATA_DIR),sep='\t')
cols=c("hsa_ID", "MI_ID","MIMAT", "shift_3p", "shift_5p", "isomiR_subs", "isomir_sequence", 
 "type", "category", "ontology", "Term", "numDEInCat", "numInCat","Pvalue","FDR", "FoldEnrich",  "genes" )

fwrite(GOres[,mget(cols)],file=sprintf('%s/19_miR_targets/target_isomiR_GO_both_3UTR_clean.txt',DATA_DIR),sep='\t')


length(unique(GOres[,MIMAT])) #26



isomir_DE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable2B_V2_isomirDE.tsv",EVO_IMMUNO_POP),sep='\t')
isomir_DE_table[isomir=='chr5_1_1_-_AGGCAGTGTATTGCTAGCGGCTGTT_0;_hsa-miR-449c-5p_MIMAT0010251_MI0003823',]
#          hsa_ID        MIMAT     MI_ID shift_5p shift_3p isomiR_subs           isomir_sequence canonical baseRatio   beta_LPS pvalue_LPS   fdr_LPS beta_PAM3CSK4 pvalue_PAM3CSK4 fdr_PAM3CSK4
# 1: hsa-miR-449c-5p MIMAT0010251 MI0003823        1        1             AGGCAGUGUAUUGCUAGCGGCUGUU     FALSE  0.827805 -0.0142207  0.3482661 0.5199183   -0.01501847        0.326978     0.497293
#      beta_R848  pvalue_R848     fdr_R848   beta_IAV   pvalue_IAV     fdr_IAV bestModel ProbModel                                                                         isomir
# 1: -0.06944274 3.418518e-06 1.919053e-05 -0.1210892 4.125069e-13 4.08322e-12        11 0.9909434 chr5_1_1_-_AGGCAGTGTATTGCTAGCGGCTGTT_0;_hsa-miR-449c-5p_MIMAT0010251_MI0003823

luq(merge(numbers_DT,isomiR_DE_table[bestModel!='',],by.x='isomiR',by.y='isomir')$MIMAT.x)
luq(merge(GOres,isomiR_DE_table[bestModel!='',],by.x='isomiR',by.y='isomir')$MIMAT.x) # 10
isomiR_DE_table[isomir%in%isomiR_ID & (bestModel!=''),luq(MIMAT)] # 290


# 10/290 = 3.4%
# allGOterms=fread(sprintf('%s/Maxime/Shared/allGOterms_EnsGRC37_13042017.txt',EVO_IMMUNO_POP))
# geneAnnot=fread(sprintf('%s/Maxime/Shared/GeneAnnotation_hg37_ens70.txt',EVO_IMMUNO_POP))

numbers_DT=fread(file=sprintf('%s/19_miR_targets/target_isomiR_number_different.txt',DATA_DIR),sep='\t')


numbers_DT_2=merge(numbers_DT, isomir_DE_table,by.x='ID',by.y='isomir')
numbers_DT_2[,pct_lost:= lost/(common+lost)]
numbers_DT_2[,pct_diff:= (lost+new)/(common+lost+new)]
numbers_DT_2[,max_beta:= pmax(abs(beta_LPS),abs(beta_PAM3CSK4),abs(beta_R848),abs(beta_IAV))]
fisher.test(table(unique(numbers_DT_2[,hsa_ID])%in%numbers_DT_2[pct_diff>.01,hsa_ID],
				  unique(numbers_DT_2[,hsa_ID])%in%numbers_DT_2[!is.na(numbers_DT_2$bestModel),hsa_ID])))

tab=table(ifelse(unique(numbers_DT_2[,hsa_ID])%in%numbers_DT_2[pct_diff>.01,hsa_ID],'diff_target','same targets'),
      ifelse(unique(numbers_DT_2[,hsa_ID])%in%numbers_DT_2[!is.na(numbers_DT_2$bestModel),hsa_ID],'DE isomiR','non DE'))
tab
#                 DE isomiR non DE
#   diff_target         98     28
#   same targets       166     85    
fisher.test(tab)
98+28/sum(tab) # = 33% des miRNAs ont un isomiR fonctionel

table(ifelse(!is.na(numbers_DT_2$bestModel),'DE isomiR','non DE'),ifelse(numbers_DT_2$pct_diff>.01,'diff_target','same targets'))

str(sapply(GO_lost,length))

NbtargetPerIsomiR=merge(numbers_DT,isomiR_Annot[,mget(c('ID','mirID','hsa_ID','miRNA_arm_num','mean_isoMiR_NS','End_site_shift','End_site_shift_template','Start_site_shift'))],by='ID')
NbtargetPerIsomiR[,Pct_new:=new/(new+common)]
NbtargetPerIsomiR[,Pct_lost:=lost/(lost+common)]
NbtargetPerIsomiR[,Pct_common:= common/(common+lost+new)]
NbtargetPerIsomiR[,Pct_diff:= 1-Pct_common]
NbtargetPerIsomiR[,mean(Pct_diff>0)] # 0.4271293
NbtargetPerIsomiR[,mean(Pct_diff>0.01)] # 0.1949527

NbtargetPerIsomiR[,Start_site_shift:=factor(Start_site_shift,c('< -2','-2','-1','0','1','2','> 2'))]
NbtargetPerIsomiR[,End_site_shift:=factor(End_site_shift,c('< -2','-2','-1','0','1','2','> 2'))]

GO_count=GOres[,.(count=length(category)),by=.(isomiR)]
NbtargetPerIsomiR[,NbGO_diff:=0]
NbtargetPerIsomiR[match(GO_count$isomiR,ID),NbGO_diff:=GO_count$count]
GO_count_lost_gain=GOres[,.(count=length(category)),by=.(isomiR,type)]

NbtargetPerIsomiR[,NbGO_lost:=0]
NbtargetPerIsomiR[match(GO_count_lost_gain[type=='lost',isomiR],ID),NbGO_lost:=GO_count_lost_gain[type=='lost',count]]
NbtargetPerIsomiR[,NbGO_gain:=0]
NbtargetPerIsomiR[match(GO_count_lost_gain[type=='gained',isomiR],ID),NbGO_gain:=GO_count_lost_gain[type=='gained',count]]
NbtargetPerIsomiR[,isomiR_expression_NS:=cut(mean_isoMiR_NS,c(0,1,5,10,50,100,1000,Inf))]
Levels=c('<1','[1,5]',']5,10]',']10,50]',']50,100]',']100,1000]','> 1000')
NbtargetPerIsomiR[,isomiR_expression_NS:=factor(Levels[as.numeric(isomiR_expression_NS)],Levels,ordered=T)]
NbtargetPerIsomiR[,mean(Pct_new,na.rm=T),by= isomiR_expression_NS][order(isomiR_expression_NS),]
#  isomiR_expression_NS         V1
# 1:                   <1 0.24017669
# 2:                [1,5] 0.15503283
# 3:               ]5,10] 0.11270162
# 4:              ]10,50] 0.11780161
# 5:             ]50,100] 0.08647810
# 6:           ]100,1000] 0.05344067
# 7:               > 1000 0.06624363
NbtargetPerIsomiR[,mean(Pct_lost,na.rm=T),by= isomiR_expression_NS][order(isomiR_expression_NS),]
#    isomiR_expression_NS         V1
# 1:                   <1 0.24464584
# 2:                [1,5] 0.14803697
# 3:               ]5,10] 0.10487283
# 4:              ]10,50] 0.11128138
# 5:             ]50,100] 0.07683372
# 6:           ]100,1000] 0.06096289
# 7:               > 1000 0.06317190
NbtargetPerIsomiR[,mean(Pct_diff>.01,na.rm=T),by= isomiR_expression_NS][order(isomiR_expression_NS),]
#    isomiR_expression_NS         V1
# 1:                   <1 0.37500000
# 2:                [1,5] 0.24842767
# 3:               ]5,10] 0.17297297
# 4:              ]10,50] 0.18098160
# 5:             ]50,100] 0.12500000
# 6:           ]100,1000] 0.09954751
# 7:               > 1000 0.10204082

range(NbtargetPerIsomiR[Start_site_shift!=0,Pct_diff])
# 0.4592793 0.9555679
range(NbtargetPerIsomiR[Start_site_shift!=0,Pct_lost])
# 0.2834746 0.9381906

NbtargetPerIsomiR[,.(Pct_lost=mean(Pct_lost,na.rm=T),Pct_new=mean(Pct_new,na.rm=T),
					Pct_any_lost=mean(Pct_lost>.01,na.rm=T),Pct_any_new=mean(Pct_new>.01,na.rm=T),
					Pct_GO_lost=mean(NbGO_lost>0,na.rm=T),Pct_GO_new=mean(NbGO_new>0,na.rm=T)),by= isomiR_expression_NS][order(isomiR_expression_NS),]

NbIsomiRwithChangesInTargets_signed=NbtargetPerIsomiR[,.(Pct_lost_noGO=mean(Pct_lost>.01 & NbGO_lost==0,na.rm=T),Pct_new_noGO=mean(Pct_new>.01 & NbGO_gain==0,na.rm=T),
													Pct_lost_GO=mean(Pct_lost>.01 & NbGO_lost>0,na.rm=T),Pct_new_GO=mean(Pct_new>.01 & NbGO_gain>0,na.rm=T)),by= isomiR_expression_NS][order(isomiR_expression_NS),]

NbIsomiRwithChangesInTargets=NbtargetPerIsomiR[,.(Pct_diff_noGO=mean(Pct_diff>.01 & NbGO_diff==0,na.rm=T),Pct_diff_GO=mean(Pct_diff>.01 & NbGO_diff>0,na.rm=T)),by= isomiR_expression_NS][order(isomiR_expression_NS),]

NbIsomiRwithChangesInTargets=melt(NbIsomiRwithChangesInTargets,id.vars='isomiR_expression_NS')
NbIsomiRwithChangesInTargets[,change:=gsub('Pct_(.*)_(.*)','\\1',variable)]
NbIsomiRwithChangesInTargets[,GO:=gsub('Pct_(.*)_(.*)','\\2',variable)]
dcast(NbIsomiRwithChangesInTargets,isomiR_expression_NS+change~GO)
NbIsomiRwithChangesInTargets[change=='new',change:='Function gained']
NbIsomiRwithChangesInTargets[change=='lost',change:='Function lost']


dir.create(sprintf('%s/isomiR_targets_3UTR/',FIG_DIR))

	pdf(sprintf('%s/isomiR_targets/target_lost_5p_end_shift.pdf',FIG_DIR),width=2.5,height=3)
		p <- ggplot(NbtargetPerIsomiR[End_site_shift==0,],aes(x=Start_site_shift, y=Pct_lost,fill=Start_site_shift)) + scale_fill_brewer(palette='RdBu')
		p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
		p <- p + xlab("Shift of 5' end")+ylab('Percentage of targets lost')+theme(legend.position="none")
		print(p)
		dev.off()

	pdf(sprintf('%s/isomiR_targets/target_lost_3p_end_shift.pdf',FIG_DIR),width=2.5,height=3)
		p <- ggplot(NbtargetPerIsomiR[Start_site_shift==0,],aes(x=End_site_shift, y=Pct_lost,fill=End_site_shift)) + scale_fill_brewer(palette='RdBu')
		p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
		p <- p + xlab("Shift of 3' end")+ylab('Percentage of targets lost')+theme(legend.position="none")
		print(p)
		dev.off()

	pdf(sprintf('%s/isomiR_targets/target_gained_5p_end_shift.pdf',FIG_DIR),width=2.5,height=3)
		p <- ggplot(NbtargetPerIsomiR[End_site_shift==0,],aes(x=Start_site_shift, y=Pct_new,fill=Start_site_shift)) + scale_fill_brewer(palette='RdBu')
		p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
		p <- p + xlab("Shift of 3' end")+ylab('Percentage of targets that are isomiR-specific')+theme(legend.position="none")
		print(p)
		dev.off()

	pdf(sprintf('%s/isomiR_targets/target_gained_3p_end_shift.pdf',FIG_DIR),width=2.5,height=3)
		p <- ggplot(NbtargetPerIsomiR[Start_site_shift==0,],aes(x=End_site_shift, y=Pct_new,fill=End_site_shift)) + scale_fill_brewer(palette='RdBu')
		p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
		p <- p + xlab("Shift of 3' end")+ylab('Percentage of targets that are isomiR-specific')+theme(legend.position="none")
		print(p)
		dev.off()

	pdf(sprintf('%s/isomiR_targets/target_diff_5p_end_shift.pdf',FIG_DIR),width=2.5,height=3)
		p <- ggplot(NbtargetPerIsomiR[End_site_shift==0,],aes(x=Start_site_shift, y=Pct_diff,fill=Start_site_shift)) + scale_fill_brewer(palette='RdBu')
		p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
		p <- p + xlab("Shift of 5' end")+ylab('Percentage of targets that are isomiR-specific')+theme(legend.position="none")
		print(p)
		dev.off()
	pdf(sprintf('%s/isomiR_targets/target_diff_3p_end_shift.pdf',FIG_DIR),width=2.5,height=3)
		p <- ggplot(NbtargetPerIsomiR[Start_site_shift==0,],aes(x=End_site_shift, y=Pct_diff,fill=End_site_shift)) + scale_fill_brewer(palette='RdBu')
		p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
		p <- p + xlab("Shift of 3' end")+ylab('Percentage of targets that are isomiR-specific')+theme(legend.position="none")
		print(p)
		dev.off()


	pdf(sprintf('%s/isomiR_targets/target_lost_Expression.pdf',FIG_DIR),width=2.5,height=3)
		p <- ggplot(NbtargetPerIsomiR,aes(x=isomiR_expression_NS, y=Pct_lost,fill=isomiR_expression_NS)) + scale_fill_brewer(palette='Blues')
		p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
		p <- p + xlab("steady-state isomiR expression (CPM)")+ylab('Percentage of targets lost')+theme(legend.position="none")
		print(p)
		dev.off()
# 
# 	pdf(sprintf('%s/isomiR_targets/target_gained_Expression.pdf',FIG_DIR),width=2.5,height=3)
# 		p <- ggplot(NbtargetPerIsomiR,aes(x=isomiR_expression_NS, y=Pct_new,fill=isomiR_expression_NS)) + scale_fill_brewer(palette='Blues')
# 		p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 		p <- p + xlab("steady-state isomiR expression (CPM)")+ylab('Percentage of targets that are isomiR-specific')+theme(legend.position="none")
# 		print(p)
# 		dev.off()
# 
# 	pdf(sprintf('%s/isomiR_targets/target_diff_Expression.pdf',FIG_DIR),width=2.5,height=3)
# 		p <- ggplot(NbtargetPerIsomiR,aes(x=isomiR_expression_NS, y=Pct_diff,fill=isomiR_expression_NS)) + scale_fill_brewer(palette='Blues')
# 		p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 		p <- p + xlab("steady-state isomiR expression (CPM)")+ylab('Percentage of targets that are isomiR-specific')+theme(legend.position="none")
# 		print(p)
# 		dev.off()

#

# 	pdf(sprintf('%s/isomiR_targets/target_Expression_perGO_3UTR.pdf',FIG_DIR),width=5,height=3)
# 	p <- ggplot(melt(NbIsomiRwithChangesInTargets),aes(x=isomiR_expression_NS,y=value*100,fill=as.factor(GO)))+facet_wrap(~change)
# 	p <- p + theme_classic()+geom_bar(stat='Identity') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 	p <- p + scale_fill_manual(values=c("#E31A1CAA","#1F78B4AA"))+xlab('IsomiR expression at steady state (CPM)')+ylab('Percentage of isomiRs with a significant change of targets')	
# 	print(p)
# 	dev.off()

	pdf(sprintf('%s/isomiR_targets/target_Expression_perGO_3UTR_diff.pdf',FIG_DIR),width=5,height=3)
	p <- ggplot(NbIsomiRwithChangesInTargets,aes(x=isomiR_expression_NS,y=value*100,fill=as.factor(GO)))
	p <- p + theme_classic()+geom_bar(stat='Identity') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
	p <- p + scale_fill_manual(values=c("#E31A1CAA","#1F78B4AA"))+xlab('IsomiR expression at steady state (CPM)')+ylab('Percentage of isomiRs with a significant change of targets')	
	print(p)
	dev.off()

### Details on miR-6503-3p isomiRs
isomiR_ID=isomiR_Annot[nsubs_other==0 & isomiR_type!='CAN' & grepl('miR-6503-3p',hsa_ID),][order(-mean_isoMiR_NS),ID]

for (myID in isomiR_ID){
	i=i+1
	cat(i,'')
	if(i%%100==0){	cat('\n')}
	miR_ID=isomiR_Annot[ID==myID,hsa_ID]
	miRANDA_targets=target_miRANDA_gene[miRNA== miR_ID, EnsemblGeneID]
	isomiR_targets=target_isomiRs_gene[isomirID==myID,EnsemblGeneID]
	common_targets=intersect(miRANDA_targets, isomiR_targets)
	new_targets=setdiff(miRANDA_targets, isomiR_targets)
	lost_targets=setdiff(isomiR_targets,miRANDA_targets)
	numbers[[myID]]=c(length(common_targets),length(new_targets),length(lost_targets))
	if(length(new_targets)>30){
		GO_new[[myID]]=GOSeq(new_targets,allGenes,bias.data=UTR3_length,addGenes=TRUE,addCI=TRUE)
		GO_new[[myID]]$isomiR=rep(myID,nrow(GO_new[[myID]]))
		}
	if(length(lost_targets)>30){
		GO_lost[[myID]]=GOSeq(lost_targets,allGenes,bias.data=UTR3_length,addGenes=TRUE,addCI=TRUE)
		GO_lost[[myID]]$isomiR=rep(myID,nrow(GO_lost[[myID]]))
		}
}


isomiR_Annot[,editing:=ifelse(grepl('.*([0-9]A..G).*', isomiR_subs),gsub('.*([0-9]A..G).*','\\1', isomiR_subs),'')]
isomiR_6503=isomiR_Annot[grepl('miR-6503',ID),  .(ID,editing, Start_site_shift,NS=mean_isoMiR_NS/sum(mean_isoMiR_NS),R848=mean_isoMiR_R848/sum(mean_isoMiR_R848))]


isomiR_Annot[grepl('miR-6503-3p',ID),  .(NS=sum(mean_isoMiR_NS),R848=sum(mean_isoMiR_R848),IAV=sum(mean_isoMiR_IAV)),by=.(Start_site_shift!=0,editing!='')] 
# upon stimulation by IAV there is an increase in the 5' shifted forms

isomiR_Annot[grepl('miR-6503-3p',ID),  .(NS=sum(mean_isoMiR_NS),R848=sum(mean_isoMiR_R848),IAV=sum(mean_isoMiR_IAV)),by=.(Start_site_shift,editing!='')] 
