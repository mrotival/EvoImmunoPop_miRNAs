print("loading gene data")

print("loading miRNA data")
corrected_count = fread(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""), sep="\t")

mir_count_melted = melt(corrected_count, id = c("ID"),
                               variable.name = "sample",
                               value.name = "count")
mir_count_melted[, individual := substr(sample, 1,6)]
mir_count_melted[, condition := condIndex[as.numeric(substr(sample, 8,8))]]
mir_count_melted[, population := substr(sample, 1,3)]
inverseNormalRankTransform=function(x){n=length(x);
									qnorm(frank(x,na.last=FALSE,ties.method='random')/(n+1),mean(x,na.rm=T),sd(x,na.rm=T))}
mir_count_melted[, count_transformed := inverseNormalRankTransform(count),by=.(ID)]
mir_count_melted[, condition := factor(condition,levels=condIndex,ordered=TRUE)]

miRNA_coordinate = fread(paste(EVO_IMMUNO_POP, "ERCPilot_SharedDBs/mirbase20/miRNA_mature_coordinates_strandinfo.bed", sep=""))
names(miRNA_coordinate) = c("chromosome", "start", "end", "miRNA_name", "V5", "V6", "V7", "V8")

## filter miRNAs with multiple sites on the genome
 miRNA_coordinate = miRNA_coordinate[miRNA_name %in% corrected_count_transformed$ID] #not enough
 duplicated_miRNAs = miRNA_coordinate[duplicated(miRNA_name), unique(miRNA_name)]
 miRNA_coordinate = miRNA_coordinate[!(miRNA_name %in% duplicated_miRNAs)]


SampleAnnot=fread(sprintf('%s/Maxime/Evo_Immuno_pop_data/SampleAnnot.txt',EVO_IMMUNO_POP))
Million_mRNA_reads=mean(SampleAnnot$sum_all_runs-SampleAnnot$unmapped)/1e6

MiR_Counts=fread(sprintf('%s/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/covariates/all_continuous_variables_V2.0_MR.tsv',EVO_IMMUNO_POP))
Million_miRNA_reads=mean(MiR_Counts$miRNA_aligned_read)/1e6

T_annot=fread(sprintf('%s/Maxime/Evo_Immuno_pop_data/02_TranscriptFPKM_cufflinks/TranscriptAnnot_expressed.txt',EVO_IMMUNO_POP))
T_annot$Nbread_NS= (T_annot$length/100* T_annot$NS_mean * Million_mRNA_reads)

pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/Counts_miRNA_VS_genes.pdf',EVO_IMMUNO_POP))
layout(1:3)
hist(log10(T_annot$Nbread_NS),main='Nb read per transcripts', xlab='log10(Nb read)',col='grey',xlim=c(-6,7),las=1,br=c(-100,seq(-6,7,by=0.5)))
G_annot=T_annot[,.(Nbread_NS=sum(Nbread_NS)),by=gene_id]
hist(log10(G_annot$Nbread_NS),main='Nb read per gene', xlab='log10(Nb read)',col='grey',xlim=c(-6,7),las=1,br=c(-100,seq(-6,7,by=0.5)))
M_annot=mir_count_melted[condition=='NS',.(Nbread_NS=(2^mean(count)-1)*Million_miRNA_reads),by=ID]
M_annot=M_annot[ID %in% miRNA_coordinate[chromosome!='chrX', miRNA_name],]
hist(log10(unlist(M_annot$Nbread_NS)),main='Nb read per miRNA', xlab='log10(Nb read)',col='grey',xlim=c(-6,7),las=1,br=c(-100,seq(-6,7,by=0.5)))
dev.off()


########### Compute percentage of genes with an eQTL
eQTL_cis_annot=fread('/Volumes/evo_immuno_pop/Maxime/miRNA_V2/data/15.eQTL_comparisons/eQTL_mapping/joint_eQTL/FDR_estimates_joint_eQTL.txt')
all_miR=M_annot$ID
high_count_miRs=M_annot[Nbread_NS>50,ID]

GeneAnnot=fread(sprintf('%s/Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/GeneAnnot_expressed.txt',EVO_IMMUNO_POP))
colnames(GeneAnnot)=make.names(colnames(GeneAnnot))
all_genes=GeneAnnot$Ensembl.Gene.ID[ GeneAnnot$Chromosome.Name%in%1:22 ]
non_codingRNAs=GeneAnnot$Ensembl.Gene.ID[GeneAnnot$Gene.Biotype!='protein_coding' & GeneAnnot$Chromosome.Name%in%1:22]

load(paste(HOME,'/Annotation/GOterms/GOterms_12578_genes.Rdata',sep=''))
TF_genes=intersect(allGOterms$gene[allGOterms$go=='GO:0003700'],all_genes)

GeneConservation=fread('/Volumes/@Home/Annotation/Conservation/GeneConservation_hg37_ens70_V4.txt')
LOF_intolerant=intersect(GeneConservation[pLI>.9,gene_id],all_genes)
LOF_Recessive=intersect(GeneConservation[pRec>.9,gene_id],all_genes)
LOF_Null=intersect(GeneConservation[pNull>.9,gene_id],all_genes)

percentage_all_genes_with_eQTL = replicate(1000,mean(sample(all_genes %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))
percentage_all_TF_with_eQTL = replicate(1000,mean(sample(TF_genes %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))
percentage_all_lncRNA_with_eQTL = replicate(1000,mean(sample(non_codingRNAs %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))
percentage_all_protein_coding_with_eQTL = replicate(1000,mean(sample(setdiff(all_genes,non_codingRNAs) %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))
percentage_all_LOF_intolerant_with_eQTL = replicate(1000,mean(sample(LOF_intolerant %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))
percentage_all_LOF_Recessive_with_eQTL = replicate(1000,mean(sample(LOF_Recessive %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))
percentage_all_LOF_Null_with_eQTL = replicate(1000,mean(sample(LOF_Null %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))

percentage_all_miR_with_eQTL=replicate(1000,mean(sample(all_miR %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))
percentage_high_count_miRs_with_eQTL=replicate(1000,mean(sample(high_count_miRs %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))


M_annot[,expression_group:=cut(Nbread_NS,c(0,10,50,100,200,1000,10000,Inf))]
percentage_miR_with_eQTL_by_expression=M_annot[,.(Pct_eQTL=replicate(1000,mean(sample(ID %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))),by=expression_group]

DT=data.table(class=rep(c('All','protein-coding','LOF Intolerant','LOF Recessive','LOF Null','lncRNA','TF','miRNA'),each=1000),Pct_with_eQTL=c(percentage_all_genes_with_eQTL,percentage_all_protein_coding_with_eQTL,percentage_all_LOF_intolerant_with_eQTL,percentage_all_LOF_Recessive_with_eQTL,percentage_all_LOF_Null_with_eQTL,percentage_all_lncRNA_with_eQTL,percentage_all_TF_with_eQTL,percentage_all_miR_with_eQTL))
fwrite(DT,file=sprintf("%s/Maxime/miRNA_V2/data/15.eQTL_comparisons/SourceData_Pct_eQTLs_byGeneClass_v2.tsv",EVO_IMMUNO_POP),sep='\t')

#### make the plot (needs improvement)
DT=fread(sprintf("%s/Maxime/miRNA_V2/data/15.eQTL_comparisons/SourceData_Pct_eQTLs_byGeneClass_v2.tsv",EVO_IMMUNO_POP),sep='\t')
pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/Pct_eQTLs_byGeneClass_v2.pdf',EVO_IMMUNO_POP),width=3.6,height=4.5)
library(RColorBrewer)
library(ggplot2)

colorSet=c(brewer.pal(8,'Set2'),brewer.pal(8,'Pastel2'),brewer.pal(8,'Dark2'),brewer.pal(8,'Set1'),brewer.pal(8,'Pastel1'))
#colorGeneType=c('grey',colorSet[c(1,3,38:37,33,6,25)]) # color scheme1
colorGeneType=c('grey',colorSet[c(1,3,6,2,25,36,28)]) # color scheme2

names(colorGeneType)=c('All','lncRNA','protein-coding','LOF Null','LOF Recessive','LOF Intolerant','TF','miRNA')
DT[,class:=factor(class,levels=c('All','lncRNA','protein-coding','LOF Null','LOF Recessive','LOF Intolerant','TF','miRNA'))]
p <- ggplot(DT[class!='All'],aes(y=Pct_with_eQTL,x=class,fill=class))+geom_violin(scale='width')+geom_boxplot(notch=TRUE,fill='#FFFFFF88')
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
p <- p + scale_fill_manual(values=colorGeneType)
p <- p + xlab('') + ylab('percentage of eQTL (bootsrapped)')
print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/Pct_miR_QTLs_byExpression.pdf',EVO_IMMUNO_POP),width=3.6,height=4.5)
p <- ggplot(percentage_miR_with_eQTL_by_expression,aes(y=Pct_eQTL,x=expression_group,fill=expression_group))+geom_violin(scale='width')+geom_boxplot(notch=TRUE,fill='#FFFFFF88')
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
p <- p + scale_fill_brewer(palette='Blues')
p <- p + xlab('') + ylab('percentage of eQTL (bootsrapped)')
print(p)
dev.off()

G_annot[,expression_group:=cut(Nbread_NS,c(0,10,50,100,200,1000,10000,Inf))]
percentage_eQTL_by_expression=G_annot[,.(Pct_eQTL=replicate(1000,mean(sample(gene_id %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))),by=expression_group]

pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/Pct_eQTLs_byExpression.pdf',EVO_IMMUNO_POP),width=3.6,height=4.5)
p <- ggplot(percentage_eQTL_by_expression,aes(y=Pct_eQTL,x=expression_group,fill=expression_group))+geom_violin(scale='width')+geom_boxplot(notch=TRUE,fill='#FFFFFF88')
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
p <- p + scale_fill_brewer(palette='Greens')
p <- p + xlab('') + ylab('percentage of eQTL (bootsrapped)')
print(p)
dev.off()



##### V3, resticting to gene and miRNAs with > 50 reads

all_genes = intersect(GeneAnnot$Ensembl.Gene.ID[ GeneAnnot$Chromosome.Name%in%1:22 ], G_annot[Nbread_NS>50,gene_id])
non_codingRNAs = intersect(GeneAnnot$Ensembl.Gene.ID[GeneAnnot$Gene.Biotype!='protein_coding' & GeneAnnot$Chromosome.Name%in%1:22], G_annot[Nbread_NS>50,gene_id])

load(paste(HOME,'/Annotation/GOterms/GOterms_12578_genes.Rdata',sep=''))
TF_genes=intersect(allGOterms$gene[allGOterms$go=='GO:0003700'],all_genes)

GeneConservation=fread('/Volumes/@Home/Annotation/Conservation/GeneConservation_hg37_ens70_V4.txt')
LOF_intolerant=intersect(GeneConservation[pLI>.9,gene_id],all_genes)
LOF_Recessive=intersect(GeneConservation[pRec>.9,gene_id],all_genes)
LOF_Null=intersect(GeneConservation[pNull>.9,gene_id],all_genes)

percentage_all_genes_with_eQTL = replicate(1000,mean(sample(all_genes %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))
percentage_all_TF_with_eQTL = replicate(1000,mean(sample(TF_genes %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))
percentage_all_lncRNA_with_eQTL = replicate(1000,mean(sample(non_codingRNAs %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))
percentage_all_protein_coding_with_eQTL = replicate(1000,mean(sample(setdiff(all_genes,non_codingRNAs) %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))
percentage_all_LOF_intolerant_with_eQTL = replicate(1000,mean(sample(LOF_intolerant %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))
percentage_all_LOF_Recessive_with_eQTL = replicate(1000,mean(sample(LOF_Recessive %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))
percentage_all_LOF_Null_with_eQTL = replicate(1000,mean(sample(LOF_Null %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))

percentage_high_count_miRs_with_eQTL=replicate(1000,mean(sample(high_count_miRs %in% eQTL_cis_annot[FDR<.05,feature_id],replace=T)))


DT=data.table(class=rep(c('All','protein-coding','LOF Intolerant','LOF Recessive','LOF Null','lncRNA','TF','miRNA'),each=1000),Pct_with_eQTL=c(percentage_all_genes_with_eQTL,percentage_all_protein_coding_with_eQTL,percentage_all_LOF_intolerant_with_eQTL,percentage_all_LOF_Recessive_with_eQTL,percentage_all_LOF_Null_with_eQTL,percentage_all_lncRNA_with_eQTL,percentage_all_TF_with_eQTL,percentage_high_count_miRs_with_eQTL))
fwrite(DT,file=sprintf("%s/Maxime/miRNA_V2/data/15.eQTL_comparisons/SourceData_Pct_eQTLs_byGeneClass_v3.tsv",EVO_IMMUNO_POP),sep='\t')

#### make the plot (needs improvement)
DT=fread(sprintf("%s/Maxime/miRNA_V2/data/15.eQTL_comparisons/SourceData_Pct_eQTLs_byGeneClass_v3.tsv",EVO_IMMUNO_POP),sep='\t')
pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/Pct_eQTLs_byGeneClass_v3.pdf',EVO_IMMUNO_POP),width=3.6,height=4.5)
library(RColorBrewer)
library(ggplot2)

colorSet=c(brewer.pal(8,'Set2'),brewer.pal(8,'Pastel2'),brewer.pal(8,'Dark2'),brewer.pal(8,'Set1'),brewer.pal(8,'Pastel1'))
#colorGeneType=c('grey',colorSet[c(1,3,38:37,33,6,25)]) # color scheme1
colorGeneType=c('grey',colorSet[c(1,3,6,2,25,36,28)]) # color scheme2

names(colorGeneType)=c('All','lncRNA','protein-coding','LOF Null','LOF Recessive','LOF Intolerant','TF','miRNA')
DT[,class:=factor(class,levels=c('All','lncRNA','protein-coding','LOF Null','LOF Recessive','LOF Intolerant','TF','miRNA'))]
p <- ggplot(DT[class!='All'],aes(y=Pct_with_eQTL,x=class,fill=class))+geom_violin(scale='width')+geom_boxplot(notch=TRUE,fill='#FFFFFF88')
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
p <- p + scale_fill_manual(values=colorGeneType)
p <- p + xlab('') + ylab('percentage of eQTL (bootsrapped)')
print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/Measurement_Error_VS_readCount.pdf',EVO_IMMUNO_POP),width=3.6,height=4.5)
plot(sapply(5:1000,function(i){x=rpois(100,i);sd(x)/mean(x)}),pch=16,cex=.3,ylab='Measurement error due to sampling (CV=SD/mean)',xlab='average read count')
abline(v=50,col='grey')
dev.off()