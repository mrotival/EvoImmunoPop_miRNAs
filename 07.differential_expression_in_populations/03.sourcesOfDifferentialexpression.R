
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
#

############################################################### 
###### test enrichments of popDE miRNAs in QTL and popDE ###### 
###############################################################

DF_DE=data.table(miRNA=miR_DE_table$miRNA,is_miR_DE=(miR_DE_table$bestModel!=''),Mean_NS=miR_DE_table$baseMean,AbslFC_DE=miR_DE_table$AbslFC, bestModel=miR_DE_table$bestModel)
DF_popDE=data.table(miRNA=miR_popDE_table$miRNA,is_miR_popDE=(miR_popDE_table$bestModel!=''),AbslFC_popDE=miR_popDE_table$AbslFC,pop_bestModel=miR_popDE_table$bestModel)
#DF_popDE=data.table(miRNA=miR_popDE_table$miRNA,is_miR_popDE=(miR_popDE_table$fdr_NS<.05),AbslFC_popDE=miR_popDE_table$AbslFC,pop_bestModel=miR_popDE_table$bestModel)

DF=merge(DF_DE,DF_popDE,by='miRNA')
DF[,has_QTL:= miRNA%in%miR_QTL$gene]
DF[,is_miR_DE_20:= is_miR_DE & AbslFC_DE>.2]
DF[,is_miR_DE_50:= is_miR_DE & AbslFC_DE>.5]

DF[,is_miR_popDE_20:= is_miR_popDE & AbslFC_popDE>.2]
DF[,is_miR_popDE_50:= is_miR_popDE & AbslFC_popDE>.5]

DF[,is_miR_popDE_NS_20:= is_miR_popDE & AbslFC_popDE>.2 & substr(pop_bestModel,1,1)==1]

DF[,bestModel_20:= ifelse(AbslFC_DE>.2,bestModel,'')]
DF[,is_miR_DE_LPS_20:= is_miR_DE & AbslFC_DE>.2 & substr(bestModel,2,2)==1]
DF[,is_miR_DE_PAM3CSK4_20:= is_miR_DE & AbslFC_DE>.2 & substr(bestModel,3,3)==1]
DF[,is_miR_DE_R848_20:= is_miR_DE & AbslFC_DE>.2 & substr(bestModel,4,4)==1]
DF[,is_miR_DE_IAV_20:= is_miR_DE & AbslFC_DE>.2 & substr(bestModel,5,5)==1]

DF[,is_miR_DE_LPS:= is_miR_DE & AbslFC_DE>.2 & substr(bestModel,2,2)==1]
DF[,is_miR_DE_PAM3CSK4:= is_miR_DE & AbslFC_DE>.2 & substr(bestModel,3,3)==1]
DF[,is_miR_DE_R848:= is_miR_DE & AbslFC_DE>.2 & substr(bestModel,4,4)==1]
DF[,is_miR_DE_IAV:= is_miR_DE & AbslFC_DE>.2 & substr(bestModel,5,5)==1]

table(DF$has_QTL, DF$is_miR_popDE)
#        FALSE TRUE
#  FALSE   349  187
#  TRUE     65   57

NSpopDE_DE_PAM3 = DF[is_miR_DE_PAM3CSK4_20 & is_miR_popDE_NS_20 & !has_QTL,miRNA]
NSpopDE_DE_IAV = DF[is_miR_DE_IAV_20 & is_miR_popDE_NS_20 & !has_QTL,miRNA]
 
 
coeff2odds=function(coeff){
    name=rownames(coeff)
    OR=exp(coeff[,1])
    OR_inf=exp(coeff[,1]-1.96*coeff[,2])
    OR_sup=exp(coeff[,1]+1.96*coeff[,2])
    p=coeff[,4]
    data.table(name,OR,OR_inf,OR_sup,p)
}

# no minimal FC, adjusted on Mean NS
coeff=list()
coeff[['is_miR_DE_LPS']]=summary(glm(is_miR_popDE ~ is_miR_DE_LPS + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['is_miR_DE_PAM3CSK4']]=summary(glm(is_miR_popDE ~ is_miR_DE_PAM3CSK4 + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['is_miR_DE_R848']]=summary(glm(is_miR_popDE ~ is_miR_DE_R848 + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['is_miR_DE_IAV']]=summary(glm(is_miR_popDE ~ is_miR_DE_IAV + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['has_QTL']]=summary(glm(is_miR_popDE ~ has_QTL + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['Mean_NS']]=summary(glm(is_miR_popDE ~ Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff=do.call(rbind, coeff)
OR_table=coeff2odds(coeff)
OR_table$name=c('DE miRNA (LPS)','DE miRNA (PAM3CSK4)','DE miRNA (R848)','DE miRNA (IAV)','miR-QTL','basal miRNA expression (2 fold increase)')

OR_table$name=factor(OR_table$name,levels=OR_table$name)
colors=c(colERC5[2:5],colERC[2:1])
names(colors)=OR_table$name

pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/SourcesOfpopDE_Enrichment.pdf',EVO_IMMUNO_POP),height=5,width=1.8)
#pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/SourcesOfpopDE_NS_Enrichment.pdf',EVO_IMMUNO_POP),height=5,width=1.8)

p <- ggplot(OR_table,aes(name,log2(OR),col=name)) + geom_pointrange( ymin = log2(OR_table$OR_inf),ymax= log2(OR_table$OR_sup)) + ylim(c(-1,6)) + scale_color_manual(values=colors)
p <- p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),axis.text.y = element_text(angle = 90, hjust = .5, vjust=0.5)) 
p <- p + geom_hline(aes(yintercept=0)) + xlab('')+ theme(legend.position = "none")
p <- p + theme(panel.grid.minor = element_blank())
print(p)
dev.off()
fwrite(OR_table,file=sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/SourcesOfpopDE_Enrichment.txt',EVO_IMMUNO_POP),sep='\t')
#fwrite(OR_table,file=sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/SourcesOfpopDE_NS_Enrichment.txt',EVO_IMMUNO_POP),sep='\t')
OR_table_miR=OR_table
# minimal FC 0.2, adjusted on Mean NS
coeff=list()
coeff[['is_miR_DE_LPS_20']]=summary(glm(is_miR_popDE_20 ~ is_miR_DE_LPS_20 + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['is_miR_DE_PAM3CSK4_20']]=summary(glm(is_miR_popDE_20 ~ is_miR_DE_PAM3CSK4_20 + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['is_miR_DE_R848_20']]=summary(glm(is_miR_popDE_20 ~ is_miR_DE_R848_20 + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['is_miR_DE_IAV_20']]=summary(glm(is_miR_popDE_20 ~ is_miR_DE_IAV_20 + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['has_QTL']]=summary(glm(is_miR_popDE_20 ~ has_QTL + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['Mean_NS']]=summary(glm(is_miR_popDE_20 ~ Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff=do.call(rbind, coeff)
OR_table=coeff2odds(coeff)
OR_table$name=c('DE miRNA (LPS)','DE miRNA (PAM3CSK4)','DE miRNA (R848)','DE miRNA (IAV)','miR-QTL','basal miRNA expression (2 fold increase)')

OR_table$name=factor(OR_table$name,levels=OR_table$name)
colors=c(colERC5[2:5],colERC[2:1])
names(colors)=OR_table$name

pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/SourcesOfpopDE_Enrichment_strong.pdf',EVO_IMMUNO_POP),height=5,width=1.8)
p <- ggplot(OR_table,aes(name,log2(OR),col=name)) + geom_pointrange( ymin = log2(OR_table$OR_inf),ymax= log2(OR_table$OR_sup)) + ylim(c(-1,6)) + scale_color_manual(values=colors)
p <- p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),axis.text.y = element_text(angle = 90, hjust = .5, vjust=0.5)) 
p <- p + geom_hline(aes(yintercept=0)) + xlab('')+ theme(legend.position = "none")
p <- p + theme(panel.grid.minor = element_blank())
print(p)
dev.off()
fwrite(OR_table,file=sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/SourcesOfpopDE_Enrichment_strong.txt',EVO_IMMUNO_POP),sep='\t')
OR_table_miR_strong=OR_table

#####################################
###### same thing with isomiRs ######
#####################################

DF_DE=data.table(miRNA=isomiR_DE_table$miRNA,isomiR=isomiR_DE_table$isomir,canonical=isomiR_DE_table$canonical,is_isomiR_DE=(isomiR_DE_table$bestModel!=''),Mean_NS=isomiR_DE_table$baseRatio,AbslFC_DE=isomiR_DE_table$AbslFC,bestModel=isomiR_DE_table$bestModel)
DF_popDE=data.table(miRNA=isomiR_popDE_table$miRNA,isomiR=isomiR_popDE_table$isomir,is_isomiR_popDE=(isomiR_popDE_table$bestModel!=''),AbslFC_popDE=isomiR_popDE_table$AbslFC)
#DF_popDE=data.table(miRNA=isomiR_popDE_table$miRNA,isomiR=isomiR_popDE_table$isomir,is_isomiR_popDE=(isomiR_popDE_table$fdr_NS<0.05),AbslFC_popDE=isomiR_popDE_table$AbslFC)
DF=merge(DF_DE,DF_popDE,by='isomiR')
DF[,has_QTL:= isomiR%in%isomiR_QTL$gene]

DF[,is_isomiR_DE_1:= is_isomiR_DE & AbslFC_DE>.01]
DF[,is_isomiR_DE_5:= is_isomiR_DE & AbslFC_DE>.05]

DF[,is_isomiR_popDE_1:= is_isomiR_popDE & AbslFC_popDE>.01]
DF[,is_isomiR_popDE_5:= is_isomiR_popDE & AbslFC_popDE>.05]

DF[,bestModel_1:= ifelse(AbslFC_DE>.01,bestModel,'')]
DF[,bestModel_5:= ifelse(AbslFC_DE>.05,bestModel,'')]
DF[,is_isomiR_DE_LPS_1:= is_isomiR_DE & AbslFC_DE>.01 & substr(bestModel,2,2)==1]
DF[,is_isomiR_DE_PAM3CSK4_1:= is_isomiR_DE & AbslFC_DE>.01 & substr(bestModel,3,3)==1]
DF[,is_isomiR_DE_R848_1:= is_isomiR_DE & AbslFC_DE>.01 & substr(bestModel,4,4)==1]
DF[,is_isomiR_DE_IAV_1:= is_isomiR_DE & AbslFC_DE>.01 & substr(bestModel,5,5)==1]

DF[,is_isomiR_DE_LPS:= is_isomiR_DE & substr(bestModel,2,2)==1]
DF[,is_isomiR_DE_PAM3CSK4:= is_isomiR_DE &  substr(bestModel,3,3)==1]
DF[,is_isomiR_DE_R848:= is_isomiR_DE &  substr(bestModel,4,4)==1]
DF[,is_isomiR_DE_IAV:= is_isomiR_DE & substr(bestModel,5,5)==1]


coeff2odds=function(coeff){
    name=rownames(coeff)
    OR=exp(coeff[,1])
    OR_inf=exp(coeff[,1]-1.96*coeff[,2])
    OR_sup=exp(coeff[,1]+1.96*coeff[,2])
    p=coeff[,4]
    data.table(name,OR,OR_inf,OR_sup,p)
}

# no minimal ratio, adjusted on Mean NS
coeff=list()
coeff[['is_isomiR_DE_LPS']]=summary(glm(is_isomiR_popDE ~ is_isomiR_DE_LPS + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['is_isomiR_DE_PAM3CSK4']]=summary(glm(is_isomiR_popDE ~ is_isomiR_DE_PAM3CSK4 + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['is_isomiR_DE_R848']]=summary(glm(is_isomiR_popDE ~ is_isomiR_DE_R848 + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['is_isomiR_DE_IAV']]=summary(glm(is_isomiR_popDE ~ is_isomiR_DE_IAV + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['has_QTL']]=summary(glm(is_isomiR_popDE ~ has_QTL + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['Mean_NS']]=summary(glm(is_isomiR_popDE ~ Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff=do.call(rbind, coeff)
OR_table=coeff2odds(coeff)
OR_table$name=c('DE isomiR (LPS)','DE isomiR (PAM3CSK4)','DE isomiR (R848)','DE isomiR (IAV)','isomiR-QTL','basal isomiRNA (2 fold increase)')

OR_table$name=factor(OR_table$name,levels=OR_table$name)
colors=c(colERC5[2:5],colERC[2:1])
names(colors)=OR_table$name

pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/SourcesOfpopDE_isomiR_Enrichment.pdf',EVO_IMMUNO_POP),height=5,width=1.8)
#pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/SourcesOfpopDE_NS_isomiR_Enrichment.pdf',EVO_IMMUNO_POP),height=5,width=1.8)
p <- ggplot(OR_table,aes(name,log2(OR),col=name)) + geom_pointrange( ymin = log2(OR_table$OR_inf),ymax= log2(OR_table$OR_sup)) + ylim(c(-1,6)) + scale_color_manual(values=colors)
p <- p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),axis.text.y = element_text(angle = 90, hjust = .5, vjust=0.5)) 
p <- p + geom_hline(aes(yintercept=0)) + xlab('')+ theme(legend.position = "none")
p <- p + theme(panel.grid.minor = element_blank())
print(p)
dev.off()
fwrite(OR_table,file=sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/SourcesOfpopDE_isomiR_Enrichment.txt',EVO_IMMUNO_POP),sep='\t')
#fwrite(OR_table,file=sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/SourcesOfpopDE_NS_isomiR_Enrichment.txt',EVO_IMMUNO_POP),sep='\t')
OR_table_isomiR=OR_table

# minimal Delta ratio = 1%, adjusted on Mean NS
coeff=list()
coeff[['is_isomiR_DE_LPS_1']]=summary(glm(is_isomiR_popDE_1 ~ is_isomiR_DE_LPS_1 + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['is_isomiR_DE_PAM3CSK4_1']]=summary(glm(is_isomiR_popDE_1 ~ is_isomiR_DE_PAM3CSK4_1 + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['is_isomiR_DE_R848_1']]=summary(glm(is_isomiR_popDE_1 ~ is_isomiR_DE_R848_1 + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['is_isomiR_DE_IAV_1']]=summary(glm(is_isomiR_popDE_1 ~ is_isomiR_DE_IAV_1 + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['has_QTL']]=summary(glm(is_isomiR_popDE_1 ~ has_QTL + Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff[['Mean_NS']]=summary(glm(is_isomiR_popDE_1 ~ Mean_NS,data=DF,family=binomial))$coeff[2,]
coeff=do.call(rbind, coeff)
OR_table=coeff2odds(coeff)
OR_table$name=c('DE isomiR (LPS)','DE isomiR (PAM3CSK4)','DE isomiR (R848)','DE isomiR (IAV)', 'isomiR-QTL', 'basal isomiRNA (2 fold increase)')

OR_table$name=factor(OR_table$name,levels=OR_table$name)
colors=c(colERC5[2:5],colERC[2:1])
names(colors)=OR_table$name

pdf(sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/SourcesOfpopDE_isomiR_Enrichment_strong.pdf',EVO_IMMUNO_POP),height=5,width=1.8)
p <- ggplot(OR_table,aes(name,log2(OR),col=name)) + geom_pointrange( ymin = log2(OR_table$OR_inf),ymax= log2(OR_table$OR_sup)) + ylim(c(-1,6)) + scale_color_manual(values=colors)
p <- p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),axis.text.y = element_text(angle = 90, hjust = .5, vjust=0.5))
p <- p + geom_hline(aes(yintercept=0)) + xlab('')+ theme(legend.position = "none")
p <- p + theme(panel.grid.minor = element_blank())
print(p)
dev.off()
fwrite(OR_table,file=sprintf('%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/SourcesOfpopDE_isomiR_Enrichment_strong.txt',EVO_IMMUNO_POP),sep='\t')
OR_table_isomiR_strong=OR_table

############ all isomiRs ############
COEFF=summary(glm(is_isomiR_popDE ~ has_QTL + is_isomiR_DE + Mean_NS,data=DF,family=binomial))$coeff;coeff2odds(COEFF)
#                        OR    OR_inf    OR_sup            p
#has_QTLTRUE      18.919523 6.3020223 56.798968 1.587470e-07
#is_isomiR_DETRUE  3.074797 2.4589924  3.844817 6.793847e-23
#Mean_NS           0.910889 0.5874634  1.412375 6.766184e-01

COEFF=summary(glm(is_isomiR_popDE_1 ~ has_QTL + is_isomiR_DE_1 + Mean_NS,data=DF,family=binomial))$coeff;coeff2odds(COEFF)
#                          OR    OR_inf    OR_sup            p
#has_QTLTRUE        37.125428 13.458076 102.41415 2.925675e-12
#is_isomiR_DE_1TRUE  6.374915  4.665607   8.71045 2.867377e-31
#Mean_NS             7.201326  4.184053  12.39446 1.029812e-12

COEFF=summary(glm(is_isomiR_popDE_1 ~ has_QTL +  I(bestModel_1=='00001') + I(bestModel_1=='00010') + I(bestModel_1=='00011') + I(bestModel_1=='01110') + I(bestModel_1=='01111') + Mean_NS,data=DF,family=binomial))$coeff;coeff2odds(COEFF)
#                                     OR    OR_inf    OR_sup            p
#has_QTLTRUE                   30.976330 11.393579 84.216999 1.720690e-11
#I(bestModel_1 == "00001")TRUE  3.923814  2.350393  6.550530 1.710078e-07
#I(bestModel_1 == "00010")TRUE  5.904962  3.582245  9.733720 3.311007e-12
#I(bestModel_1 == "00011")TRUE  4.766519  3.268821  6.950427 4.868779e-16
#I(bestModel_1 == "01110")TRUE  6.484088  3.532836 11.900750 1.603088e-09
#I(bestModel_1 == "01111")TRUE  1.853572  0.619829  5.543026 2.695163e-01
#Mean_NS                        8.654073  5.082547 14.735326 1.904870e-15

COEFF=summary(glm(is_isomiR_popDE_5 ~ has_QTL +  I(bestModel_5=='00001') + I(bestModel_5=='00010') + I(bestModel_5=='00011') + I(bestModel_5=='01110') + I(bestModel_5=='01111') + Mean_NS,data=DF,family=binomial))$coeff;coeff2odds(COEFF)
#                                     OR     OR_inf     OR_sup            p
#has_QTLTRUE                   49.275586 18.8634476 128.718963 1.783078e-15
#I(bestModel_5 == "00001")TRUE  1.796088  0.3708539   8.698662 4.668715e-01
#I(bestModel_5 == "00010")TRUE  9.699796  1.1530021  81.600936 3.652605e-02
#I(bestModel_5 == "00011")TRUE  1.013117  0.1325005   7.746438 9.899814e-01
#I(bestModel_5 == "01110")TRUE 14.670074  3.8163915  56.391249 9.247930e-05
#I(bestModel_5 == "01111")TRUE 30.970789  4.7368851 202.493777 3.389131e-04
#Mean_NS                       14.052049  4.3876827  45.003275 8.581103e-06

############ Canonical isomiRs only ############

COEFF=summary(glm(is_isomiR_popDE ~ has_QTL + is_isomiR_DE + Mean_NS,data=DF[which(canonical),],family=binomial))$coeff;coeff2odds(COEFF)
#                           OR    OR_inf   OR_sup            p
#has_QTLTRUE      4.821635e+07 0.0000000      Inf 9.831712e-01
#is_isomiR_DETRUE 4.929417e+00 2.9118133 8.345023 2.865076e-09
#Mean_NS          7.874682e-01 0.3322697 1.866274 5.873168e-01

COEFF=summary(glm(is_isomiR_popDE_1 ~ has_QTL + is_isomiR_DE_1 + Mean_NS,data=DF[which(canonical),],family=binomial))$coeff;coeff2odds(COEFF)
#                             OR    OR_inf   OR_sup            p
#has_QTLTRUE        5.087127e+07 0.0000000      Inf 9.835628e-01
#is_isomiR_DE_1TRUE 3.775445e+00 2.2774472 6.258757 2.584051e-07
#Mean_NS            1.086850e+00 0.4343608 2.719495 8.587401e-01
 
COEFF=summary(glm(is_isomiR_popDE_1 ~ has_QTL +  I(bestModel_1=='00001') + I(bestModel_1=='00010') + I(bestModel_1=='00011') + I(bestModel_1=='01110') + I(bestModel_1=='01111') + Mean_NS,data=DF[which(canonical),],family=binomial))$coeff;coeff2odds(COEFF)
#                                        OR    OR_inf    OR_sup            p
#has_QTLTRUE                   4.461067e+07 0.0000000       Inf 9.833769e-01
#I(bestModel_1 == "00001")TRUE 2.307799e+00 1.0468932  5.087371 3.811368e-02
#I(bestModel_1 == "00010")TRUE 6.086400e+00 2.6951387 13.744843 1.389616e-05
#I(bestModel_1 == "00011")TRUE 4.258289e+00 2.2392877  8.097674 9.941117e-06
#I(bestModel_1 == "01110")TRUE 2.879154e+00 0.9191296  9.018889 6.948605e-02
#I(bestModel_1 == "01111")TRUE 3.398043e+00 0.7756353 14.886755 1.046099e-01
#Mean_NS                       1.105432e+00 0.4352689  2.807415 8.330501e-01
 
COEFF=summary(glm(is_isomiR_popDE_5 ~ has_QTL +  I(bestModel_5=='00001') + I(bestModel_5=='00010') + I(bestModel_5=='00011') + I(bestModel_5=='01110') + I(bestModel_5=='01111') + Mean_NS,data=DF[which(canonical),],family=binomial))$coeff;coeff2odds(COEFF)
#                                        OR    OR_inf    OR_sup            p
#has_QTLTRUE                   4.696519e+01 9.0441648 243.88420 4.646209e-06
#I(bestModel_5 == "00001")TRUE 1.005250e-07 0.0000000       Inf 9.958432e-01
#I(bestModel_5 == "00010")TRUE 1.366145e-07 0.0000000       Inf 9.983413e-01
#I(bestModel_5 == "00011")TRUE 1.118400e-07 0.0000000       Inf 9.945446e-01
#I(bestModel_5 == "01110")TRUE 1.241481e+01 1.1547517 133.47237 3.764097e-02
#I(bestModel_5 == "01111")TRUE 1.781692e+01 1.4611620 217.25365 2.399503e-02
#Mean_NS                       2.853420e+00 0.3304127  24.64192 3.404741e-01
