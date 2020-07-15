

eQTL_sharing=fread('/Volumes/evo_immuno_pop/Maxime/miRNA_V2/data/15.eQTL_comparisons/sharing_eQTL_conditions_likelihoods.tsv')
eQTL_sharing=eQTL_sharing[order(gene,-ProbModel),][!duplicated(gene),]

geneAnnot=fread(sprintf('%s/Maxime/Shared/GeneAnnotation_hg37_ens70.txt',EVO_IMMUNO_POP))
geneAnnot=geneAnnot[which(select),]
geneAnnot[,logFC_LPS:=log2(1+LPS_mean)-log2(1+NS_mean)]
geneAnnot[,logFC_PAM3CSK4:=log2(1+PAM3_mean)-log2(1+NS_mean)]
geneAnnot[,logFC_R848:=log2(1+R848_mean)-log2(1+NS_mean)]
geneAnnot[,logFC_IAV:=log2(1+Flu_mean)-log2(1+NS_mean)]
geneAnnot[,maxFC:=pmax(abs(logFC_LPS),abs(logFC_PAM3CSK4),abs(logFC_R848),abs(logFC_IAV))]

eQTL_DE=merge(eQTL_sharing,geneAnnot,by.x='gene',by.y='Ensembl.Gene.ID')
eQTL_DE[,condition_dependent:=(condition_test!=11111)]
merge(eQTL_sharing,geneAnnot,by.x='gene',by.y='Ensembl.Gene.ID')[,mean(condition_test==11111),by=cut(maxFC,c(0,0.05,0.2,0.5,1,Inf))]
range(merge(eQTL_sharing,geneAnnot,by.x='gene',by.y='Ensembl.Gene.ID')$maxFC)
eQTL_DE[,type:='gene']

#          cut        V1
# 1:   (1,Inf] 0.4522452
# 2: (0.2,0.5] 0.5426009
# 3:   (0.5,1] 0.5281501
# 4:   (0,0.2] 0.5263158

miRDE=fread('/Volumes/evo_immuno_pop/Maxime/miRNA_V2/data/00_tables_publication/SupTable2A_mirDE.tsv')
miRDE[,maxFC:=pmax(abs(log2FC_LPS),abs(log2FC_PAM3CSK4),abs(log2FC_R848),abs(log2FC_IAV))]
miRQTL=fread('/Volumes/evo_immuno_pop/Maxime/miRNA_V2/data/00_tables_publication/SupTable3B_resp-mirQTLs_Annoted_V2.tsv')
miRQTL[,miRNA:=NULL]
miRQTL_DE=merge(miRQTL,miRDE,by='miRNA',suffix=c('.QTL','.DE'))
miRQTL_DE[,condition_dependent:=(bestModel.QTL!=11111)]
miRQTL_DE[,condition_dependent:=(bestModel.QTL!=11111)]
miRQTL_DE[,type:='miR']
QTL_DE=rbind(miRQTL_DE[,.(type,maxFC,condition_dependent)],eQTL_DE[,.(type,maxFC,condition_dependent)])

mod_both=glm(condition_dependent~maxFC+type,family=binomial,data=QTL_DE)

mod_both=glm(condition_dependent~maxFC+type,family=binomial,data=QTL_DE[maxFC>0.1 & maxFC<0.3,])


mod_miR=glm(condition_dependent~maxFC,family=binomial,data=miRQTL_DE)
mod_gene=glm(condition_dependent~maxFC,family=binomial,data=eQTL_DE)
plot(predict.glm(mod_miR,newdata=data.frame(maxFC=seq(0,4,by=.1))))

mod_miRgene=glm(condition_dependent~maxFC,family=binomial,data=miRQTL_DE)
plot(seq(0,4,by=.1),invlogit(predict.glm(mod_miR,newdata=data.frame(maxFC=seq(0,4,by=.1)))),type='l')
lines(seq(0,4,by=.1),invlogit(predict.glm(mod_gene,newdata=data.frame(maxFC=seq(0,4,by=.1)))),col='red')

pdf('QTL_sharing_VS_foldChange.pdf')
layout(1:2)
par(mar=c(4,4,.2,.2))
hh=hist(merge(miRQTL,miRDE,by='miRNA')$maxFC,br=seq(0,8,by=0.1),col='#00000088',freq=F,xlim=c(0,3),xlab='maximum absolute log Fold Change',main='',las=1,add=F,axes=F)
axis(1,las=1,at=0:3);axis(2,las=1)
hist(merge(eQTL_sharing,geneAnnot,by.x='gene',by.y='Ensembl.Gene.ID')$maxFC,br=hh$br,add=T,col='#FF000088',freq=F)

legend('topright',legend=c('miRNA','gene'),fill=c('#00000088','#FF000088'),bty='n')

plot(seq(0,4,by=.1),invlogit(predict.glm(mod_both,newdata=data.frame(maxFC=seq(0,4,by=.1),type='gene'))),col='red',type='l',ylab='probability of \n condition-dependent QTL',ylim=c(0,1),,xlab='maximum absolute log Fold Change',axes=F,xlim=c(0,3))
axis(1,las=1,at=0:3);axis(2,las=1)
legend('bottomright',legend=c('miRNA','gene'),col=c('#00000088','#FF000088'),bty='n',lty=1)
lines(seq(0,4,by=.1),invlogit(predict.glm(mod_both,newdata=data.frame(maxFC=seq(0,4,by=.1),type='miR'))),col='black')
dev.off()

miRNA_count=fread('/Volumes/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.RPM.GCRL_Batch_corrected_V2.0_MR.tsv')

miRNA_logcount=fread('/Volumes/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv')
layout(matrix(1:2,1))
plot(log2(apply(miRNA_count[,-1],1,mean)),log2(apply(miRNA_count[,-1],1,sd)),xlab='mean expression',ylab='Standard deviation',axes=F)
axis(1,at=log2(c(1,10,100,1000,10000)),labels=c(1,10,100,1000,10000),las=1)
axis(2,at=log2(c(1,10,100,1000,10000)),labels=c(1,10,100,1000,10000),las=1)
plot(apply(miRNA_logcount[,-1],1,mean),apply(miRNA_logcount[,-1],1,sd),xlab='mean expression',ylab='Standard deviation',axes=F)
axis(1,at=log2(c(1,10,100,1000,10000)),labels=c(1,10,100,1000,10000),las=1)
axis(2,las=1)
