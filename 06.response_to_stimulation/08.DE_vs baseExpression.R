DATA_DIR=sprintf("%s/Maxime/miRNA_V2/data/",EVO_IMMUNO_POP)
FIG_DIR=sprintf("%s/Maxime/miRNA_V2/figures/Revisions",EVO_IMMUNO_POP)

condIndex=c("NS","LPS","PAM3CSK4","R848","IAV")

library(lme4)
library(lmerTest)

print("loading data")
corrected_count = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""), sep="\t")

mir_count_melted = melt(corrected_count, id = c("ID"),
                               variable.name = "sample",
                               value.name = "count")
mir_count_melted[, individual := substr(sample, 1,6)]
mir_count_melted[, condition := substr(sample, 8,8)]
mir_count_melted[, population := substr(sample, 1,3)]
inverseNormalRankTransform=function(x){n=length(x);
									qnorm(frank(x,na.last=FALSE,ties.method='random')/(n+1),mean(x,na.rm=T),sd(x,na.rm=T))}
mir_count_melted[, count_transformed := inverseNormalRankTransform(count),by=.(ID)]

total_count=mir_count_melted[,.(count=sum(2^count-1)),by= .(sample, individual, condition ,population)]
ggplot(total_count,aes(x= condition,y=log2(count),fill=condition))+geom_violin(scale='width')+geom_boxplot(width=.5,fill='#FFFFF88',notch=T)


miR_DE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable2A_mirDE.tsv",EVO_IMMUNO_POP),colClass=c('character',rep('numeric',13),'character','numeric'))
isomiR_DE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable2B_V2_isomirDE.tsv",EVO_IMMUNO_POP),colClass=c('character','character','character','numeric','numeric','character','character','logical',rep('numeric',13),'character','numeric','character'))

Expr_group=cut(2^miR_DE_table$baseMean-1,c(0,1,10,100,1000,Inf))
levels(Expr_group)=c('< 1','1-10','10-100','100-1000','> 1000')
miR_DE_table[,Expr_group:=Expr_group]
table(miR_DE_table[,Expr_group],ifelse(miR_DE_table[,bestModel]!='','DE',''))

#      (0,1]      (1,10]    (10,100] (100,1e+03] (1e+03,Inf] 
#          33         328         158          84          55 

# tab=table(miR_DE_table[,Expr_group],ifelse(miR_DE_table[,bestModel]!='','DE',''))
# tab=table(miR_DE_table[,Expr_group],ifelse(miR_DE_table[,log2FC_LPS]>.2,'DE up',ifelse(miR_DE_table[,log2FC_LPS]< -.2,'DE down','')))

miR_DE_melt=melt(miR_DE_table,id.vars=c('miRNA','baseMean','Expr_group','bestModel','ProbModel'),
				patterns(logFC="^log2FC_", FDR="^fdr_",pvalue='^pvalue_'),value.factor=TRUE)
miR_DE_melt[,condition:=sort(condIndex[-1])[as.numeric(variable)]]
names(colERC5)=condIndex

	pdf(sprintf('%s/DE_vs_Expr_logFC.pdf',FIG_DIR),width=5,height=6)
		p <- ggplot(miR_DE_melt,aes(x=Expr_group,y=logFC,fill=condition)) + scale_fill_manual(values=colERC5[-1])+facet_wrap(~factor(condition,condIndex[-1]))
		p <- p + theme_bw()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
		p <- p + xlab("miRNA expression level")+ylab('log FC')+theme(legend.position="none")
		print(p)
		dev.off()

	pdf(sprintf('%s/DE_vs_Expr_AbslogFC.pdf',FIG_DIR),width=2.5,height=3)
		p <- ggplot(miR_DE_melt,aes(x=Expr_group,y=logFC,fill=condition)) + scale_fill_manual(values=colERC5[-1])
		p <- p + theme_bw()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
		p <- p + xlab("miRNA expression level")+ylab('log FC')+theme(legend.position="none")+scale_y_continuous(trans='sqrt')
		print(p)		
		dev.off()


	pdf(sprintf('%s/DE_vs_Expr_scatter.pdf',FIG_DIR),width=2.5,height=3)
		p <- ggplot(miR_DE_melt,aes(x=baseMean,y=logFC,col=ifelse(FDR<.01,condition,'NS'))) + scale_colour_manual(values=colERC5)
		p <- p + theme_classic()+geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
		p <- p + xlab("miRNA expression level")+ylab('log FC')+theme(legend.position="none")
		p <- p + geom_hline(yintercept=0,col='lightgrey')
		print(p)
		dev.off()
		
	pdf(sprintf('%s/DE_vs_Expr_barplot.pdf',FIG_DIR),width=2.5,height=3)
		miR_DE_melt[,direction:=ifelse(logFC>0,'up','down')]
		PctDE=miR_DE_melt[,.(up=mean(FDR<.01 & logFC>0)*100, down=mean(FDR<.01 & logFC<0)*100),by=.(condition,Expr_group)]
		PctDE=melt(PctDE,variable.name='direction',value.name='PctDE')
		p <- ggplot(PctDE,aes(x=Expr_group,y=PctDE,fill=direction)) + scale_fill_manual(values=c("#E31A1CAA","#1F78B4AA"))+facet_wrap(~factor(condition,condIndex[-1]))
		p <- p + theme_bw()+geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
		p <- p + xlab("basal expression")+ylab('PctDE miRNA')+theme(legend.position="bottom")
		print(p)
		dev.off()
