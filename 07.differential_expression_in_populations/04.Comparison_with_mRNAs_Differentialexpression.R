library(ggplot2)


miR_QTL=fread(sprintf('%s/Maxime/miRNA_V2/data/09.mirQTL/cis_mirQTLs_with_FDR_filtered_best_miRNA_snp_association.tsv',EVO_IMMUNO_POP))
isomiR_QTL=fread(sprintf('%s/Maxime/miRNA_V2/data/10.isomiRQTL/cis_isomiRQTLs_with_FDR_filtered_best_isomiRNA_snp_association.tsv',EVO_IMMUNO_POP))

#isomiR_popDE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable3B_isomir_popDE.tsv",EVO_IMMUNO_POP),colClass=c('character','character','logical',rep('numeric',16),'character','numeric'))
#miR_popDE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable3A_mir_popDE.tsv",EVO_IMMUNO_POP),colClass=c('character',rep('numeric',16),'character','numeric'))
miR_popDE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4A_miR_popDE_V2.tsv",EVO_IMMUNO_POP),colClass=c('character',rep('numeric',20),'character','numeric'))
isomiR_popDE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4B_isomir_popDE_V2.tsv",EVO_IMMUNO_POP),colClass=c(rep('character',7),'logical',rep('numeric',20),'character','numeric','character'))
isomiR_popDE_table[,miRNA:=hsa_ID]

miR_DE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable2A_mirDE.tsv",EVO_IMMUNO_POP),colClass=c('character',rep('numeric',13),'character','numeric'))
isomiR_DE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable2B_isomirDE.tsv",EVO_IMMUNO_POP),colClass=c('character','character','logical',rep('numeric',13),'character','numeric'))

miR_DE_table[,AbslFC:=pmax(abs(log2FC_LPS),abs(log2FC_PAM3CSK4),abs(log2FC_R848),abs(log2FC_IAV))]
miR_popDE_table[,AbslFC:=pmax(abs(MeanAFB_LPS-MeanEUB_LPS),abs(MeanAFB_PAM3CSK4-MeanEUB_PAM3CSK4),abs(MeanAFB_R848-MeanEUB_R848),abs(MeanAFB_IAV-MeanEUB_IAV))]

isomiR_DE_table[,AbslFC:=pmax(abs(beta_LPS),abs(beta_PAM3CSK4),abs(beta_R848),abs(beta_IAV))]
isomiR_popDE_table[,AbslFC:=pmax(abs(MeanAFB_LPS-MeanEUB_LPS),abs(MeanAFB_PAM3CSK4-MeanEUB_PAM3CSK4),abs(MeanAFB_R848-MeanEUB_R848),abs(MeanAFB_IAV-MeanEUB_IAV))]

gene_popDE=fread(sprintf('/Volumes/@Home/03_Analysis/DiffExpression/DiffExp_population_ByCondition_geneLevel.txt',EVO_IMMUNO_POP))
gene_popDE=gene_popDE[,mget(c("Ensembl.Gene.ID","Chromosome.Name","Strand","Gene.Start..bp.","Gene.End..bp.","Associated.Gene.Name","Description","Gene.Biotype",
'Mean_AFB_1','Mean_EUB_1',"Beta_NS_EUBvsAFB","Pval_NS_EUBvsAFB",
'Mean_AFB_2','Mean_EUB_2',"Beta_LPS_EUBvsAFB","Pval_LPS_EUBvsAFB",
'Mean_AFB_3','Mean_EUB_3',"Beta_PAM3CSK4_EUBvsAFB","Pval_PAM3CSK4_EUBvsAFB",
'Mean_AFB_4','Mean_EUB_4',"Beta_R848_EUBvsAFB","Pval_R848_EUBvsAFB",
'Mean_AFB_5','Mean_EUB_5',"Beta_IAV_EUBvsAFB","Pval_IAV_EUBvsAFB"))]

gene_DE=fread(sprintf('/Volumes/@Home/03_Analysis/DiffExpression/DiffExp_population_ByCondition_geneLevel.txt',EVO_IMMUNO_POP))[,mget(c("Ensembl.Gene.ID","NS_mean","LPS_mean","PAM3_mean","R848_mean","Flu_mean"))]
colnames(gene_DE)=c("Ensembl.Gene.ID","Base_mean","LPS_mean","PAM3CSK4_mean","R848_mean","IAV_mean")
gene_DE[,log2FC_LPS:=log2(LPS_mean+1)-log2(Base_mean+1)]
gene_DE[,log2FC_LPS:=log2(PAM3CSK4_mean+1)-log2(Base_mean+1)]
gene_DE[,log2FC_R848:=log2(R848_mean+1)-log2(Base_mean+1)]
gene_DE[,log2FC_IAV:=log2(IAV_mean+1)-log2(Base_mean+1)]

gene_popDE_melt=melt(gene_popDE,id.vars=c('Ensembl.Gene.ID',"Gene.Biotype"),
                                patterns("^Mean_AFB_","^Mean_EUB_","^Pval_"),value.name=c("MeanAFB","MeanEUB",'pvalue'),variable.name='condition')
gene_popDE_melt[,condition:=condIndex[as.numeric(condition)]]
gene_popDE_melt[,MeanAFB:=log2(1+MeanAFB)]
gene_popDE_melt[,MeanEUB:=log2(1+MeanEUB)]

colnames(gene_popDE_melt)[1]='feature'

miR_popDE_melt=melt(miR_popDE_table,id.vars=c('miRNA'),measure=patterns("^MeanAFB_","^MeanEUB_","^pvalue_"),value.name=c("MeanAFB","MeanEUB",'pvalue'),variable.name='condition')
miR_popDE_melt[,condition:=condIndex[as.numeric(condition)]]
miR_popDE_melt[,Gene.Biotype:="mature_miRNA"]
colnames(miR_popDE_melt)[1]='feature'

popDE_melt=rbind(gene_popDE_melt,miR_popDE_melt)
popDE_melt[,Abs_logFC:=abs(MeanAFB-MeanEUB)]
popDE_melt[,condition:=factor(condition,levels=condIndex)]

eQTL_cis_annot=fread('/Volumes/evo_immuno_pop/Maxime/miRNA_V2/data/15.eQTL_comparisons/eQTL_mapping/joint_eQTL/FDR_estimates_joint_eQTL.txt')

count_biotype=table(popDE_melt$Gene.Biotype)/5
frequent_biotype=names(count_biotype)[count_biotype>300]
frequent_biotype=c('mature_miRNA','protein_coding')

popDE_melt=popDE_melt[feature%in% eQTL_cis_annot$feature_id & Gene.Biotype %in%frequent_biotype,]
popDE_melt[,hasQTL:=ifelse(feature%in%eQTL_cis_annot[FDR<.05,feature_id],'QTL','no QTL')]
popDE_melt[,FDR:=p.adjust(pvalue,'fdr')]
popDE_cast=dcast(popDE_melt,Gene.Biotype+feature+hasQTL~condition,value.var=c('Abs_logFC','FDR'))

wilcox.test(popDE_cast[Gene.Biotype=='mature_miRNA',Abs_logFC_NS],popDE_cast[Gene.Biotype=='protein_coding',Abs_logFC_NS])$p.value

wilcox.test(popDE_cast[Gene.Biotype=='mature_miRNA',Abs_logFC_NS],popDE_cast[Gene.Biotype=='mature_miRNA',Abs_logFC_LPS])$p.value
wilcox.test(popDE_cast[Gene.Biotype=='mature_miRNA',Abs_logFC_NS],popDE_cast[Gene.Biotype=='mature_miRNA',Abs_logFC_R848])$p.value
wilcox.test(popDE_cast[Gene.Biotype=='mature_miRNA',Abs_logFC_NS],popDE_cast[Gene.Biotype=='mature_miRNA',Abs_logFC_PAM3CSK4])$p.value
wilcox.test(popDE_cast[Gene.Biotype=='mature_miRNA',Abs_logFC_NS],popDE_cast[Gene.Biotype=='mature_miRNA',Abs_logFC_IAV])$p.value

wilcox.test(popDE_cast[Gene.Biotype=='protein_coding',Abs_logFC_NS],popDE_cast[Gene.Biotype=='protein_coding',Abs_logFC_LPS])$p.value
wilcox.test(popDE_cast[Gene.Biotype=='protein_coding',Abs_logFC_NS],popDE_cast[Gene.Biotype=='protein_coding',Abs_logFC_R848])$p.value
wilcox.test(popDE_cast[Gene.Biotype=='protein_coding',Abs_logFC_NS],popDE_cast[Gene.Biotype=='protein_coding',Abs_logFC_PAM3CSK4])$p.value
wilcox.test(popDE_cast[Gene.Biotype=='protein_coding',Abs_logFC_NS],popDE_cast[Gene.Biotype=='protein_coding',Abs_logFC_IAV])$p.value

popDE_melt[,.(popDE=any(Abs_logFC>.2)),by=.(Gene.Biotype, feature)][,.(Pct_popDE=mean(popDE)),by= Gene.Biotype]
#     Gene.Biotype  Pct_popDE
#1: protein_coding  0.5086060
#2:   mature_miRNA  0.2700186
popDE_melt[,mean(Abs_logFC>.2),by= .(Gene.Biotype, condition)]
#     Gene.Biotype condition  Pct_popDE
# 1: protein_coding        NS 0.1154719
# 2: protein_coding       LPS 0.2628145
# 3: protein_coding  PAM3CSK4 0.1777946
# 4: protein_coding      R848 0.1681483
# 5: protein_coding       IAV 0.2530736
# 6:   mature_miRNA        NS 0.1675978
# 7:   mature_miRNA       LPS 0.1694600
# 8:   mature_miRNA  PAM3CSK4 0.1545624
# 9:   mature_miRNA      R848 0.1843575
#10:   mature_miRNA       IAV 0.1675978
popDE_melt[,.(popDE=any(FDR <0.01 & Abs_logFC>.2)),by=.(Gene.Biotype, feature)][,mean(popDE),by= Gene.Biotype]
popDE_melt[,mean(Abs_logFC>.2 & FDR <0.01),by= .(Gene.Biotype, condition)]

popDE_melt[,outlier:= Abs_logFC > median(Abs_logFC)+1.5*IQR(Abs_logFC),by=.(Gene.Biotype,condition)]  
pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/Comparison_protein_coding/popDE_Delta_by_cond_small.pdf',EVO_IMMUNO_POP),height=4.2,width=6)
   p <- ggplot(popDE_melt[Gene.Biotype%in%frequent_biotype,],aes(x=Gene.Biotype,y=pmin(Abs_logFC,.5),fill=condition))
   p <- p + geom_violin(scale='width') + geom_boxplot(notch=TRUE,fill='#FFFFFF88',outlier.size=0)
   p <- p + geom_jitter(data = function(x) dplyr::filter_(x, ~ outlier),col='#00000011',size=.5,width=.2) 
   p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
   p <- p + scale_fill_manual(values=colERC5)
   p <- p + xlab('feature type') + ylab('Absolute log2FC')
   p <- p +facet_grid(~condition)
   print(p)
  dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/Comparison_protein_coding/popDE_Delta_by_cond_QTL.pdf',EVO_IMMUNO_POP),height=4)
   p <- ggplot(popDE_melt[Gene.Biotype%in%frequent_biotype,],aes(x=paste(hasQTL,Gene.Biotype),y=pmin(Abs_logFC,.5),fill=condition))
   p <- p +geom_violin(scale='width') + geom_boxplot(notch=TRUE,fill='#FFFFFF88')
   p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
   p <- p + scale_fill_manual(values=colERC5)
   p <- p + xlab('feature type') + ylab('Absolute log2FC') 
   p <- p +facet_grid(~condition)
   print(p)
  dev.off()
  
  
Pct_DE_1pctFDR=popDE_melt[Gene.Biotype %in%frequent_biotype ,mean(FDR<0.01),by=.(Gene.Biotype,condition)][order(condition, Gene.Biotype),]
Pct_DE_logFC.2=popDE_melt[Gene.Biotype %in%frequent_biotype ,mean(FDR<0.01 & Abs_logFC>.2),by=.(Gene.Biotype,condition)][order(condition, Gene.Biotype),]

Pct_DE_1pctFDR_any=popDE_cast[Gene.Biotype %in%frequent_biotype ,mean(FDR_NS<0.01 | FDR_LPS<0.01 | FDR_PAM3CSK4<0.01 | FDR_R848<0.01 | FDR_IAV<0.01),by=.(Gene.Biotype)][order(Gene.Biotype),]
Pct_DE_logFC.2_any=popDE_cast[Gene.Biotype %in%frequent_biotype ,mean(FDR_NS<0.01 & Abs_logFC_NS>.2 | FDR_LPS<0.01 & Abs_logFC_LPS>.2 | FDR_PAM3CSK4<0.01 & Abs_logFC_PAM3CSK4>.2 | FDR_R848<0.01 & Abs_logFC_R848>.2 | FDR_IAV<0.01 & Abs_logFC_IAV>.2),by=.(Gene.Biotype)][order(Gene.Biotype),]


pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/Comparison_protein_coding/popDE_FDR1pct_by_cond.pdf',EVO_IMMUNO_POP))
 p <- ggplot(Pct_DE_1pctFDR,aes(x=Gene.Biotype,y=V1,fill=condition))+geom_bar(stat='Identity',position='dodge')
 p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
 p <- p + scale_fill_manual(values=colERC5)
 p <- p + xlab('feature type') + ylab('% differential 1%FDR')
 print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/Comparison_protein_coding/popDE_logFC.2_by_cond.pdf',EVO_IMMUNO_POP))
 p <- ggplot(Pct_DE_logFC.2,aes(x=Gene.Biotype,y=V1,fill= condition))+geom_bar(stat='Identity',position='dodge')
 p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
 p <- p + scale_fill_manual(values=colERC5)
 p <- p + xlab('feature type') + ylab('% differential logFC>.2') 
 print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/Comparison_protein_coding/popDE_any_FDR1pct.pdf',EVO_IMMUNO_POP))
  p <- ggplot(Pct_DE_1pctFDR_any, aes(x=Gene.Biotype,y=V1))+geom_bar(stat='Identity',position='dodge')
  p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
  p <- p + xlab('feature type') + ylab('% differential 1%FDR')
 print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/Comparison_protein_coding/popDE_any_logFC.2.pdf',EVO_IMMUNO_POP))
 p <- ggplot(Pct_DE_logFC.2_any, aes(x=Gene.Biotype,y=V1))+geom_bar(stat='Identity',position='dodge')
 p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
 p <- p + xlab('feature type') + ylab('% differential logFC>.2') 
 print(p)
dev.off()




Pct_DE_1pctFDR=popDE_melt[Gene.Biotype %in%frequent_biotype ,mean(FDR<0.01),by=.(Gene.Biotype,condition,hasQTL)][order(condition, Gene.Biotype),]
Pct_DE_logFC.2=popDE_melt[Gene.Biotype %in%frequent_biotype ,mean(FDR<0.01 & Abs_logFC>.2),by=.(Gene.Biotype,condition,hasQTL)][order(condition, Gene.Biotype),]

Pct_DE_1pctFDR_any=popDE_cast[Gene.Biotype %in%frequent_biotype ,mean(FDR_NS<0.01 | FDR_LPS<0.01 | FDR_PAM3CSK4<0.01 | FDR_R848<0.01 | FDR_IAV<0.01),by=.(Gene.Biotype,hasQTL)][order(Gene.Biotype,hasQTL),]
Pct_DE_logFC.2_any=popDE_cast[Gene.Biotype %in%frequent_biotype ,mean(FDR_NS<0.01 & Abs_logFC_NS>.2 | FDR_LPS<0.01 & Abs_logFC_LPS>.2 | FDR_PAM3CSK4<0.01 & Abs_logFC_PAM3CSK4>.2 | FDR_R848<0.01 & Abs_logFC_R848>.2 | FDR_IAV<0.01 & Abs_logFC_IAV>.2),by=.(Gene.Biotype,hasQTL)][order(Gene.Biotype,hasQTL),]


pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/Comparison_protein_coding/popDE_FDR1pct_by_cond_and_QTL.pdf',EVO_IMMUNO_POP))
 p <- ggplot(Pct_DE_1pctFDR,aes(x=paste(Gene.Biotype,hasQTL),y=V1,fill=condition))+geom_bar(stat='Identity',position='dodge')
 p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
 p <- p + scale_fill_manual(values=colERC5)
 p <- p + xlab('feature type') + ylab('% differential 1%FDR')
 print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/Comparison_protein_coding/popDE_logFC.2_by_cond_and_QTL.pdf',EVO_IMMUNO_POP))
 p <- ggplot(Pct_DE_logFC.2,aes(x=paste(Gene.Biotype,hasQTL),y=V1,fill= condition))+geom_bar(stat='Identity',position='dodge')
 p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
 p <- p + scale_fill_manual(values=colERC5)
 p <- p + xlab('feature type') + ylab('% differential logFC>.2') 
 print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/Comparison_protein_coding/popDE_any_FDR1pct_by_QTL.pdf',EVO_IMMUNO_POP))
  p <- ggplot(Pct_DE_1pctFDR_any, aes(x=paste(Gene.Biotype,hasQTL),y=V1))+geom_bar(stat='Identity',position='dodge')
  p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
  p <- p + xlab('feature type') + ylab('% differential 1%FDR')

 print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/Comparison_protein_coding/popDE_any_logFC.2_by_QTL.pdf',EVO_IMMUNO_POP))
 p <- ggplot(Pct_DE_logFC.2_any, aes(x=paste(Gene.Biotype,hasQTL),y=V1))+geom_bar(stat='Identity',position='dodge')
 p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
 p <- p + xlab('feature type') + ylab('% differential logFC>.2') 
 print(p)
dev.off()


#gene_popDE_Pval[,fdr:=p.adjust(pvalue)]
#gene_popDE_Pval[,variable:=gsub('Pval_','',variable)]
#gene_popDE_Pval=dcast(gene_popDE_Pval,Ensembl.Gene.ID~variable,value.var=c('Pvalue','FDR'))

colnames(gene_popDE)=c("Ensembl.Gene.ID","Chromosome.Name","Strand","Gene.Start..bp.","Gene.End..bp.","Associated.Gene.Name","Description","Gene.Biotype",
'Mean_AFB_NS','Mean_EUB_NS',"Delta_NS_EUBvsAFB","Pval_NS_EUBvsAFB",
'Mean_AFB_LPS','Mean_EUB_LPS',"Delta_LPS_EUBvsAFB","Pval_LPS_EUBvsAFB",
'Mean_AFB_PAM3CSK4','Mean_EUB_PAM3CSK4',"Delta_PAM3CSK4_EUBvsAFB","Pval_PAM3CSK4_EUBvsAFB",
'Mean_AFB_R848','Mean_EUB_R848',"Delta_R848_EUBvsAFB","Pval_R848_EUBvsAFB",
'Mean_AFB_IAV','Mean_EUB_IAV',"Delta_IAV_EUBvsAFB","Pval_IAV_EUBvsAFB")

gene_eQTL=fread(sprintf('%s/Maxime/miRNA_V2/data/09.mirQTL/cis_mirQTLs_with_FDR_filtered_best_miRNA_snp_association.tsv',EVO_IMMUNO_POP))


miR_DE_melt=melt(miR_DE_table,id.vars=c('miRNA','baseMean'),measure=patterns("^log2FC_", "^pvalue_",'^fdr_'),value.name=c('log2FC','pvalue','fdr'),variable.name='condition')
miR_DE_melt[,condition:=condIndex[1+as.numeric(condition)]]



DiffResponse_population_ByCondition_geneLevel_withGxE.txt

/Volumes/@Home/03_Analysis/DiffExpression/DiffExp_population_ByCondition_geneLevel.txt

