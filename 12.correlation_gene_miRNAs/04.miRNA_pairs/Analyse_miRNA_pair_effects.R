miR_pair_DIR=paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/19_miR_targets/mir_pairs_colocalization/CorrelationGeneExpression/",sep='')
FilesList=dir(miR_pair_DIR,pat='withInt')

genes_expression = fread(paste(EVO_IMMUNO_POP, "Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/FPKM_matrix.txt", sep=""))
genes_expression = melt(genes_expression, id = "ID", variable.name = "sample", value.name = "expression")
genes_expression[, individual := substr(sample, 1,6)]
genes_expression[, condition := substr(sample, 8,8)]
genes_expression[, population := substr(sample, 1,3)]
# genes_expression = genes_expression[condition == cond_to_investigate]
names(genes_expression)[1] = "gene"

genes_names = genes_expression[,unique(gene)]

miR_pair_gene=fread(sprintf('%s/19_miR_targets/mir_pairs_colocalization/miR_pairs_genecount.txt',DATA_DIR))
miR_pair_gene=miR_pair_gene[Padj_Coloc_enrich<.01 & NbColoc_gene>100,]

miR_pair=fread(sprintf('%s/19_miR_targets/mir_pairs_colocalization/miR_pairs.txt',DATA_DIR))
miR_pair=miR_pair[paste(miRNA.x,miRNA.y)%chin%paste(miR_pair_gene$miRNA.x,miR_pair_gene$miRNA.y),]
miR_pair=miR_pair[EnsemblGeneID %chin% genes_names,]

miR_pair=miR_pair[order(miRNA.x,miRNA.y,EnsemblGeneID,-NbColoc),]
miR_pair=miR_pair[!duplicated(paste(miRNA.x,miRNA.y,EnsemblGeneID)),]

FilesList[i]

if(!grepl(FilesList[i]))
miR1=gsub("JointEffect_geneRPKM_miR1_(.*)_miR2_(.*)_cond(.*).tsv",'\\1',FilesList[i])
miR2=gsub("JointEffect_geneRPKM_miR1_(.*)_miR2_(.*)_cond(.*).tsv",'\\2',FilesList[i])
cond=gsub("JointEffect_geneRPKM_miR1_(.*)_miR2_(.*)_cond(.*).tsv",'\\3',FilesList[i])

miR_pair_cor=list()
WithOutInt=grep('cond.\\.tsv',FilesList)
WithInt=grep('cond._withInt\\.tsv',FilesList)

for (i in WithOutInt){
	miR_pair_cor[[i]]=fread(paste(miR_pair_DIR,FilesList[i],sep='/'), sep="\t")
	}
	miR_pair_cor=rbindlist(miR_pair_cor)
	
	miR_pair_cor_annot=merge(miR_pair_cor,miR_pair[,mget(c('EnsemblGeneID','miRNA.x','miRNA.y','NbColoc'))],by.x=c('gene','miR1','miR2'),by.y=c('EnsemblGeneID','miRNA.x','miRNA.y'),all.x=T)
	miR_pair_cor_annot[is.na(NbColoc),NbColoc:=0]

table(miR_pair_cor_annot$NbColoc)

miR_pair_cor_test=miR_pair_cor_annot[,.(miRNA_effect=wilcox.test(pvalue[NbColoc==0],pvalue[NbColoc>0],alt='greater')$p.value),by=.(miR1,miR2,condition)]

library(qvalue)
pi0_est=miR_pair_cor_annot[,.(pi1=replicate(100,1-pi0est(sample(pvalue,replace=T))$pi0)),by=.(ifelse(NbColoc>0,'target','non-target'),condition)]
pi0_est=miR_pair_cor_annot[,.(pi1=replicate(100,1-mean(sample(pvalue,replace=T)>.05)/0.95)),by=.(ifelse(NbColoc>0,'target','non-target'),condition)]

condIndex=c("NS", "LPS", "PAM3CSK4", "R848", "IAV")
Bootstrap_pi0=miR_pair_cor_annot[,.(pi1=replicate(100,1-pi0est(sample(pvalue,replace=T))$pi0)),by=.(ifelse(NbColoc>0,'target','non-target'),condition)]
Bootstrap_pi0[,condition:=condIndex[condition]]
fwrite(Bootstrap_pi0,file=paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/19_miR_targets/mir_pairs_colocalization/Bootstrap_pi0_target_vs_nontargets.txt",sep=''))

colBoxplot=c("#969696", "#525252", "#FB9A99", "#E31A1C", "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4", "#CAB2D6", "#6A3D9A")
names(colBoxplot)=paste(condIndex[rep(1:5,each=2)],c('non-target','target'))
pdf(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/19_miR_targets/mir_pairs_colocalization/Bootstrap_pi0_target_vs_nontargets.pdf",sep=''),width=4, heigh=4)
    p <- ggplot(Bootstrap_pi0,aes(x=ifelse,y=pi1,fill=paste(cond,ifelse)))+geom_violin(scale='width')+geom_boxplot(width=.5,fill='#FFFFFF88',notch=T)
    p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
    p <- p + xlab('') + ylab('% of genes associated\n with expression of the miRNA pair')
#     p <- p + scale_fill_manual(values=colBoxplot)+ facet_wrap(~cond)
    print(p)
dev.off()

Bootstrap_pi0=miR_pair_cor_annot[,.(pi1=replicate(100,1-pi0est(sample(pvalue,replace=T))$pi0)),by=.(ifelse(NbColoc>0,'target','non-target'),condition)]
Bootstrap_pi0[,condition:=condIndex[condition]]
fwrite(Bootstrap_pi0,file=paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/19_miR_targets/mir_pairs_colocalization/Bootstrap_pi0_target_vs_nontargets_allCond.txt",sep=''))

pdf(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/19_miR_targets/mir_pairs_colocalization/Bootstrap_pi0_target_vs_nontargets_allCond.pdf",sep=''),width=4, heigh=4)
    p <- ggplot(Bootstrap_pi0,aes(x=ifelse,y=pi1,fill=paste(condition,ifelse)))+geom_violin(scale='width')+geom_boxplot(width=.5,fill='#FFFFFF88',notch=T)
    p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
    p <- p + xlab('') + ylab('% of genes associated\n with expression of the miRNA pair')
    p <- p + scale_fill_manual(values=colBoxplot)+ facet_grid(~factor(condition,condIndex)) + theme(legend.position='top')
    print(p)
dev.off()

Bootstrap_pi0=miR_pair_cor_annot[,.(pi1=replicate(100,1-mean(sample(pvalue,replace=T)>.05)/0.95)),by=.(ifelse(NbColoc>0,'target','non-target'),condition)]
Bootstrap_pi0[,condition:=condIndex[condition]]
fwrite(Bootstrap_pi0,file=paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/19_miR_targets/mir_pairs_colocalization/Bootstrap_pi0_target_vs_nontargets_allCond_altMethod.txt",sep=''))

pdf(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/19_miR_targets/mir_pairs_colocalization/Bootstrap_pi0_target_vs_nontargets_allCond_altMethod.pdf",sep=''),width=6, heigh=3.5)
    p <- ggplot(Bootstrap_pi0,aes(x=ifelse,y=pi1,fill=paste(condition,ifelse)))+geom_violin(scale='width')+geom_boxplot(width=.5,fill='#FFFFFF88',notch=T)
    p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
    p <- p + xlab('') + ylab('% of genes associated with\nexpression of the miRNA pair')
    p <- p + scale_fill_manual(values=colBoxplot)+ facet_grid(~factor(condition,condIndex)) + theme(legend.position='top')
    print(p)
dev.off()

miR_pair_cor_annot[,is_target:=ifelse(NbColoc>0,'target','non-target')]
miR_pair_cor_annot[,condition:=condIndex[condition]]
pdf(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/19_miR_targets/mir_pairs_colocalization/InteractionEffect_target_vs_nontargets_allCond.pdf",sep=''),width=6, heigh=3.5)
    p <- ggplot(miR_pair_cor_annot,aes(x=is_target,y=miR12_car,fill=paste(condition,is_target)))+geom_violin(scale='width')+geom_boxplot(width=.5,fill='#FFFFFF88',notch=T)
    p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
    p <- p + xlab('') + ylab('Interaction Effect on gene expression')
    p <- p + scale_fill_manual(values=colBoxplot)+ facet_grid(~factor(condition,condIndex)) + theme(legend.position='top')
    print(p)
dev.off()

pdf(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/19_miR_targets/mir_pairs_colocalization/InteractionEffect_target_vs_nontargets_allCond.pdf",sep=''),width=6, heigh=3.5)
    p <- ggplot(miR_pair_cor_annot,aes(x=is_target,y=miR12_car,fill=paste(condition,is_target)))+geom_violin(scale='width')+geom_boxplot(width=.5,fill='#FFFFFF88',notch=T)
    p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
    p <- p + xlab('') + ylab('Interaction Effect on gene expression')
    p <- p + scale_fill_manual(values=colBoxplot)+ facet_grid(~factor(condition,condIndex)) + theme(legend.position='top')
    print(p)
dev.off()

pdf(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/19_miR_targets/mir_pairs_colocalization/pvalue_hist_effect_on_targets.pdf",sep=''),width=6, heigh=3.5)
hist(miR_pair_cor_annot[NbColoc==0 & condition=='NS',pvalue],col='grey',br=seq(0,1,l=20))
hist(miR_pair_cor_annot[NbColoc>0 & condition=='NS',pvalue],col='grey',br=seq(0,1,l=20))
hist(miR_pair_cor_annot[NbColoc==0 & condition=='IAV',pvalue],col='grey',br=seq(0,1,l=20))
hist(miR_pair_cor_annot[NbColoc>0 & condition=='IAV',pvalue],col='grey',br=seq(0,1,l=20))
dev.off()


miR_pair_cor_annot[pvalue<0.05,mean(miR12_car<0),by=.(condition,is_target)]


pdf(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/19_miR_targets/mir_pairs_colocalization/Sign_InteractionEffect_target_vs_nontargets_allCond.pdf",sep=''),width=6, heigh=3.5)
p <- ggplot(miR_pair_cor_annot,aes(x=as.factor(is_target),fill=as.factor(sign(miR12_car))))
p <- p + theme_classic()+geom_bar(stat='count',position="fill")
p <- p + scale_fill_manual(values=c("#1F78B4AA","#E31A1CAA"))+xlab('Does the gene containa pair of colocalized targets')+ylab('Effect of miRNA interaction on gene expression')	
p <- p + theme(legend.position="top")+ theme(legend.position="top", axis.text.x = element_text(angle = 45, hjust = 1)) + facet_grid(~factor(condition,condIndex))
print(p)
dev.off()


miR_pair_cor_annot[,.(pi1=replicate(100,1-pi0est(sample(pvalue,replace=T))$pi0)),by=.(ifelse(NbColoc>0,'target','non-target'),condition)
