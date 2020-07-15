
DATA_DIR=sprintf("%s/Maxime/miRNA_V2/data/",EVO_IMMUNO_POP)
FIG_DIR=sprintf("%s/Maxime/miRNA_V2/figures/Revisions",EVO_IMMUNO_POP)

library(rtracklayer)

REfSeq_gff <- rtracklayer::import(sprintf("%s/19_miR_targets/GRCh37_latest_genomic_RefSeq.gff.gz",DATA_DIR))
names(elementMetadata(REfSeq_gff)@listData)
# REfSeq_gff=as.data.frame(REfSeq_gff)
# what types of region are present
table(as.character(REfSeq_gff$type[REfSeq_gff$gbkey%in%c('mRNA','ncRNA','Gene','exon','misc_RNA','precursor_RNA')]))
REfSeq_gff=REfSeq_gff[REfSeq_gff$gbkey%in%c('mRNA','ncRNA','Gene','exon','misc_RNA','precursor_RNA'),]

REfSeq_gff[,c('transcript_id','gene','Dbxref')]

ID_matching=data.table(DB_ID = unlist(REfSeq_gff$Dbxref),
				  Refseq = rep(REfSeq_gff$transcript_id,sapply(REfSeq_gff$Dbxref,length)),
				  Symbol = rep(REfSeq_gff$gene,sapply(REfSeq_gff$Dbxref,length)))

ID_matching[,DB_ID:=gsub('HGNC:HGNC','HGNC',DB_ID)]
ID_matching=ID_matching[!duplicated(ID_matching),]
ID_matching[,DB:=gsub('(.*):(.*)','\\1',DB_ID)]
ID_matching[,ID:=gsub('(.*):(.*)','\\2',DB_ID)]
ID_matching[,RefSeq_stable:=gsub('(.*)\\..*','\\1',Refseq)]

fwrite(ID_matching,file=sprintf('%s/19_miR_targets/Refseq_Entrez_HUGO_IDmatching.txt',DATA_DIR),sep='\t')
ID_matching=fread(sprintf('%s/19_miR_targets/Refseq_Entrez_HUGO_IDmatching.txt',DATA_DIR))
###########################################################
###				miR DB V6 targets 						###
###########################################################

miRDB=fread(sprintf("%s/19_miR_targets/miRDB_v6.0_prediction_result.txt-1",DATA_DIR))
colnames(miRDB)=c('miR_ID','target_Refseq','target_Score')
miRDB=miRDB[substr(miR_ID,1,3)=='hsa',]

# count miRNA targets at various thresholds
miRDB_target_number=miRDB[,.(nbTot=length(target_Score),over70=sum(target_Score>70),over80=sum(target_Score>80),over90=sum(target_Score>90),over95=sum(target_Score>95)),by=miR_ID]
melt(miRDB_target_number,value.name='nb')[,.(mean=mean(nb),median=median(nb),pct5=quantile(nb,.05),pct95= quantile(nb,.95)),by=variable]
#    variable       mean median  pct5   pct95
# 1:    nbTot 1270.98682 1017.5 69.75 3267.00
# 2:   over70  551.10354  387.5 18.00 1675.25
# 3:   over80  305.81325  190.0  7.00 1043.75
# 4:   over90  106.99812   51.0  1.00  405.00
# 5:   over95   36.51995   11.0  0.00  151.25

# get number of targets per miR in miR DB V6
p <- ggplot(melt(miRDB_target_number,value.name='nb'),aes(x=log2(1+nb),fill=variable)) + geom_histogram()
p <- p + scale_color_brewer(palette="Greens")+ scale_fill_brewer(palette="Greens") + theme_classic() 
p <- p + ylab('Number of microRNAs')+ xlab('log2(Number of targets +1)')
pdf(sprintf('%s/QC_miR_targets/NumberOtTargets_mirDB_V6.pdf',FIG_DIR),width=3.5,height=2.5)
print(p)
dev.off()

# add gene names to miR_DB
nondup=which(!duplicated(paste(ID_matching$RefSeq_stable,ID_matching$Symbol)))
miRDB=merge(miRDB,ID_matching[nondup,.(RefSeq_stable,Symbol)],by.x='target_Refseq',by.y='RefSeq_stable',all.x=TRUE)
fwrite(miRDB[which(!is.na(Symbol)),],file=sprintf('%s/19_miR_targets/miRDB_v6.0_prediction_result_with_HUGO.txt',DATA_DIR),sep='\t')

###########################################################
###				targetScan targets 	(default)			###
###########################################################
target_TS=fread(sprintf('unzip -p %s/19_miR_targets/targetScan/v7.1/Predicted_Targets_Context_Scores.default_predictions.txt.zip',DATA_DIR))
colnames(target_TS)=gsub('.','_',make.names(gsub('++','Plus',colnames(target_TS),fixed=T)),fixed=TRUE)
target_TS=target_TS[Gene_Tax_ID==9606,]
target_TS[,target_score:=sum(weighted_contextPlus_score_percentile/100),by=.(Transcript_ID,miRNA,Gene_ID,Gene_Symbol)]
target_TS_gene=target_TS[order(miRNA,Gene_ID,Gene_Symbol,-target_score),][which(!duplicated(paste(miRNA,Gene_ID,Gene_Symbol))),]

fwrite(target_TS_gene[,.(miRNA,Gene_ID,Gene_Symbol,target_score)],file=sprintf('%s/19_miR_targets/targetScan/v7.1/DefaultSites_gene_withHUGO.txt',DATA_DIR),sep='\t')

# count miRNA targets at various thresholds default
TScan_target_number=target_TS[,.(nbTot=length(unique(Gene_ID)),
								over70=length(unique(Gene_ID[weighted_contextPlus_score_percentile>70])),
								over80=length(unique(Gene_ID[weighted_contextPlus_score_percentile>80])),
								over90=length(unique(Gene_ID[weighted_contextPlus_score_percentile>90])),
								over95=length(unique(Gene_ID[weighted_contextPlus_score_percentile>95]))),by=miRNA]
								
# get number of targets per miR in TargetScan default							
melt(TScan_target_number,value.name='nb')[,.(mean=mean(nb),median=median(nb),pct5=quantile(nb,.05),pct95= quantile(nb,.95)),by=variable]
#    variable     mean median pct5 pct95
# 1:    nbTot 631.1558    497   34  1388
# 2:   over70 493.5545    375   31  1198
# 3:   over80 424.1277    303   26  1069
# 4:   over90 277.0935    205   16   706
# 5:   over95 146.4829    113    9   377

p <- ggplot(melt(TScan_target_number,value.name='nb'),aes(x=log2(1+nb),fill=variable)) + geom_histogram()
p <- p + scale_color_brewer(palette="Greens")+ scale_fill_brewer(palette="Greens") + theme_classic() 
p <- p + ylab('Number of microRNAs')+ xlab('log2(Number of targets +1)')
pdf(sprintf('%s/QC_miR_targets/NumberOtTargets_TargetScanDefault.pdf',FIG_DIR),width=3.5,height=2.5)
print(p)
dev.off()

# TScan_info=fread(sprintf('unzip -p %s/19_miR_targets/targetScan/v7.1/Predicted_Targets_Info.default_predictions.txt.zip',DATA_DIR))
# colnames(TScan_info)=gsub('.','_',make.names(colnames(TScan_info)),fixed=TRUE)
# TScan_info=TScan_info[Species_ID==9606,]

###########################################################
###				targetScan FULL 	(default)			###
###########################################################
target_TSfull_C=fread(sprintf('unzip -p %s/19_miR_targets/targetScan/v7.1/Conserved_Site_Context_Scores.txt.zip',DATA_DIR))
colnames(target_TSfull_C)=gsub('.','_',make.names(gsub('++','Plus',colnames(target_TSfull_C),fixed=T)),fixed=TRUE)
target_TSfull_C=target_TSfull_C[Gene_Tax_ID==9606,]

target_TSfull_NC=fread(sprintf('unzip -p %s/19_miR_targets/targetScan/v7.1/NonConserved_Site_Context_Scores.txt.zip | grep -e 9606',DATA_DIR))
colnames(target_TSfull_NC)=colnames(target_TSfull_C)
target_TSfull_NC=target_TSfull_NC[,conserved:='no']

target_TSfull_C=target_TSfull_C[,conserved:='yes']
target_TSfull_C=rbind(target_TSfull_C,target_TSfull_NC)
rm(target_TSfull_NC);gc()

fwrite(target_TSfull_C,file=sprintf('%s/19_miR_targets/targetScan/v7.1/allSites_Context_Scores_withHUGO.txt',DATA_DIR),sep='\t')

target_TSfull_C=target_TSfull_C[which(!is.na(weighted_contextPlus_score_percentile)),]
target_TSfull_C=target_TSfull_C[which(Gene_Tax_ID==9606),]


target_TSfull_C[,target_score:=max(weighted_contextPlus_score_percentile/100),by=.(Transcript_ID,miRNA,Gene_ID,Gene_Symbol)]
target_TSfull_C[,nbTarget:=length(weighted_contextPlus_score_percentile),by=.(Transcript_ID,miRNA,Gene_ID,Gene_Symbol)]
target_TSfull_gene=target_TSfull_C[order(miRNA,Gene_ID,Gene_Symbol,-target_score),][which(!duplicated(paste(miRNA,Gene_ID,Gene_Symbol))),]

fwrite(target_TSfull_gene[,.(miRNA,Gene_ID,Gene_Symbol,nbTarget,target_score)],file=sprintf('%s/19_miR_targets/targetScan/v7.1/allSites_gene_withHUGO.txt',DATA_DIR),sep='\t')

# count miRNA targets at various thresholds default
TScan_target_number=target_TSfull_C[,.(nbTot=length(unique(Gene_ID)),
								over70=length(unique(Gene_ID[weighted_contextPlus_score_percentile>70])),
								over80=length(unique(Gene_ID[weighted_contextPlus_score_percentile>80])),
								over90=length(unique(Gene_ID[weighted_contextPlus_score_percentile>90])),
								over95=length(unique(Gene_ID[weighted_contextPlus_score_percentile>95]))),by=miRNA]

TScan_target_number_C=target_TSfull_C[conserved=='yes',.(nbTot=length(unique(Gene_ID)),
								over70=length(unique(Gene_ID[weighted_contextPlus_score_percentile>70])),
								over80=length(unique(Gene_ID[weighted_contextPlus_score_percentile>80])),
								over90=length(unique(Gene_ID[weighted_contextPlus_score_percentile>90])),
								over95=length(unique(Gene_ID[weighted_contextPlus_score_percentile>95]))),by=miRNA]

TScan_target_number_NC=target_TSfull_C[conserved=='no',.(nbTot=length(unique(Gene_ID)),
								over70=length(unique(Gene_ID[weighted_contextPlus_score_percentile>70])),
								over80=length(unique(Gene_ID[weighted_contextPlus_score_percentile>80])),
								over90=length(unique(Gene_ID[weighted_contextPlus_score_percentile>90])),
								over95=length(unique(Gene_ID[weighted_contextPlus_score_percentile>95]))),by=miRNA]

# get number of targets per miR in TargetScan default							
melt(TScan_target_number,value.name='nb')[,.(mean=mean(nb),median=median(nb),pct5=quantile(nb,.05),pct95= quantile(nb,.95)),by=variable]
#    variable      mean median   pct5   pct95
# 1:    nbTot 3882.5345 3934.0 720.75 6260.75
# 2:   over70 2147.8035 2127.0 435.50 3637.50
# 3:   over80 1636.6140 1651.0 335.25 2746.00
# 4:   over90  866.3638  868.0 161.50 1534.00
# 5:   over95  387.8937  380.5  68.00  712.50
melt(TScan_target_number_C,value.name='nb')[,.(mean=mean(nb),median=median(nb),pct5=quantile(nb,.05),pct95= quantile(nb,.95)),by=variable]
#   variable     mean median pct5   pct95
# 1:    nbTot 579.1330  449.0    3 1382.00
# 2:   over70 451.1182  330.0    3 1075.25
# 3:   over80 387.0419  284.0    3  898.75
# 4:   over90 252.7734  189.0    2  640.25
# 5:   over95 133.5517  102.5    2  349.25
melt(TScan_target_number_NC,value.name='nb')[,.(mean=mean(nb),median=median(nb),pct5=quantile(nb,.05),pct95= quantile(nb,.95)),by=variable]
#   variable      mean median   pct5   pct95
# 1:    nbTot 3823.6032 3874.5 720.75 6243.00
# 2:   over70 2087.9290 2058.5 432.50 3599.00
# 3:   over80 1582.0794 1578.0 329.50 2718.00
# 4:   over90  828.6285  814.0 156.25 1511.25
# 5:   over95  367.4616  357.0  64.25  698.75

p <- ggplot(melt(TScan_target_number,value.name='nb'),aes(x=log2(1+nb),fill=variable)) + geom_histogram()
p <- p + scale_color_brewer(palette="Greens")+ scale_fill_brewer(palette="Greens") + theme_classic() 
p <- p + ylab('Number of microRNAs')+ xlab('log2(Number of targets +1)')
pdf(sprintf('%s/QC_miR_targets/NumberOtTargets_TargetScan_All.pdf',FIG_DIR),width=3.5,height=2.5)
print(p)
dev.off()


###########################################################
###				MIRZA-G targets 						###
###########################################################


# add gene names to miRZaG
miRzaG=fread(sprintf("%s/19_miR_targets/mirza-g_all_mirnas_per_gene_scores.tab",DATA_DIR))
colnames(miRzaG)=c('Entrez','miRNA','MIRZA','MIRZAG')
nondup=which(!duplicated(paste(ID_matching[DB=='GeneID',ID],ID_matching[DB=='GeneID',Symbol])))
nondup=which(ID_matching$DB=='GeneID')[nondup]
miRzaG=merge(miRzaG,ID_matching[nondup,.(ID=as.numeric(ID),Symbol)],by.x='Entrez',by.y='ID',all.x=TRUE)
fwrite(miRzaG,file=sprintf('%s/19_miR_targets/mirza-g_all_mirnas_per_gene_scores_with_HUGO.txt',DATA_DIR),sep='\t')

MIRZA_target_number=miRzaG[,.(nbTot=length(Entrez),
								MIRZAover.2=sum(MIRZA>.2), 
								MIRZAover.5=sum(MIRZA>.5), 
								MIRZAover1=sum(MIRZA>1),
								MIRZAGover.2=sum(MIRZAG>.2,na.rm=T),
								MIRZAGover.5=sum(MIRZAG>.5,na.rm=T),
								MIRZAGover1=sum(MIRZAG>1,na.rm=T)),by=miRNA]
melt(MIRZA_target_number,value.name='nb')[,.(mean=mean(nb),median=median(nb),pct5=quantile(nb,.05),pct95= quantile(nb,.95)),by=variable]

# get number of targets per miR in miRANDA
#        variable       mean median  pct5  pct95
# 1:        nbTot 1850.31933   1769 371.0 3465.5
# 2:  MIRZAover.2 1217.67328   1109 174.5 2554.0
# 3:  MIRZAover.5  236.97939    159  11.0  744.5
# 4:   MIRZAover1   24.71217      8   0.0  103.5
# 5: MIRZAGover.2  943.69934    823  93.5 2129.5
# 6: MIRZAGover.5  174.04006    106   5.0  569.0
# 7:  MIRZAGover1   18.42824      5   0.0   77.0

# get number of targets per miR in targetScan
p <- ggplot(melt(MIRZA_target_number,value.name='nb'),aes(x=log2(1+nb),fill=variable)) + geom_histogram()
p <- p + scale_color_brewer(palette="Greens")+ scale_fill_brewer(palette="Greens") + theme_classic() 
p <- p + ylab('Number of microRNAs')+ xlab('log2(Number of targets +1)')
pdf(sprintf('%s/QC_miR_targets/NumberOtTargets_MIRZAG.pdf',FIG_DIR),width=3.5,height=2.5)
print(p)
dev.off()

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

miR_count_gene=target_miRANDA_gene[,.(miRNA,EnsemblGeneID,GeneNames,NbTarget)]
miR_count_gene=miR_count_gene[,.(NbTarget = sum(NbTarget>0)),by= miRNA]

###############################################################
###		compare targets between methods, and score them		###
###############################################################
targets=list()
targets[['miRANDA']]=fread(sprintf('%s/19_miR_targets/target_miRANDA_gene_with_HUGO.txt',DATA_DIR),sep='\t')
targets[['miRANDA']]=targets[['miRANDA']][,.(miRNA,GeneNames,target_score)]
targets[['miRANDA']][,algo:='miRANDA']

targets[['MIRZA']]=fread(sprintf('%s/19_miR_targets/mirza-g_all_mirnas_per_gene_scores_with_HUGO.txt',DATA_DIR),sep='\t')
targets[['MIRZAG']]=targets[['MIRZA']]
targets[['MIRZA']]=targets[['MIRZA']][,.(miRNA,Symbol,MIRZA)]
targets[['MIRZAG']]=targets[['MIRZAG']][,.(miRNA,Symbol,MIRZAG)]
targets[['MIRZA']][,algo:='MIRZA']
targets[['MIRZAG']][,algo:='MIRZAG']
targets[['MIRZAG']]=targets[['MIRZAG']][which(!is.na(MIRZAG)),]

targets[['targetScan_default']]=fread(sprintf('%s/19_miR_targets/targetScan/v7.1/DefaultSites_gene_withHUGO.txt',DATA_DIR),sep='\t')
targets[['targetScan_default']]=targets[['targetScan_default']][,.( miRNA,Gene_Symbol,target_score)]
targets[['targetScan_default']][,algo:='targetScan_default']

targets[['targetScan']]=fread(sprintf('%s/19_miR_targets/targetScan/v7.1/allSites_gene_withHUGO.txt',DATA_DIR),sep='\t')
targets[['targetScan']]=targets[['targetScan']][,.( miRNA,Gene_Symbol,target_score)]
targets[['targetScan']][,algo:='targetScan']

targets[['miRDB']]=fread(sprintf('%s/19_miR_targets/miRDB_v6.0_prediction_result_with_HUGO.txt',DATA_DIR),sep='\t')
targets[['miRDB']]=targets[['miRDB']][,.(miR_ID,Symbol,target_Score)]
targets[['miRDB']][,algo:='miRDB']
targets[['miRDB']][,target_Score:=target_Score/100]

for(i in names(targets)){
	colnames(targets[[i]])=c('miRNA','Symbol','target_Score','algo')
}
targets=rbindlist(targets)

count_miRNA=targets[algo!='targetScan_default',.(miRNA=unique(miRNA)),by= algo][,.N,by=miRNA]
common_miRNA=count_miRNA[N==5,miRNA]
count_gene=targets[algo!='targetScan_default',.(Symbol=unique(Symbol)),by= algo][,.N,by=Symbol]
common_gene=count_gene[N==5,Symbol]

targets_comparable=targets[algo!='targetScan_default' & miRNA %chin% common_miRNA & Symbol %chin% common_gene,]
fwrite(targets_comparable,file=sprintf('%s/19_miR_targets/Comparison_5methods.txt',DATA_DIR),sep='\t')
targets_comparable=fread(file=sprintf('%s/19_miR_targets/Comparison_5methods.txt',DATA_DIR),sep='\t')


targets_compared=dcast(targets_comparable,miRNA+Symbol~ algo,value.var='target_Score',fun.aggregate=max)
targets_compared[,nbAlgo:=(is.finite(MIRZA)+is.finite(MIRZAG)+is.finite(miRANDA)+is.finite(miRDB)+is.finite(targetScan))]
targets_compared[, code := paste(Symbol, gsub("-", "_", miRNA))]

pdf(sprintf('%s/QC_miR_targets/Score_VS_numberOfMethods_MIRZA.pdf',FIG_DIR),width=3.5,height=2.5)
p <- ggplot(targets_compared[is.finite(MIRZA),],aes(x= as.factor(nbAlgo),y= MIRZA,fill= as.factor(nbAlgo)))+theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T)+scale_fill_brewer('Blues')+scale_y_sqrt()+xlab('Number of databases where interaction is predicted')
print(p)
dev.off()

pdf(sprintf('%s/QC_miR_targets/Score_VS_numberOfMethods_MIRZAG.pdf',FIG_DIR),width=3.5,height=2.5)
p <- ggplot(targets_compared[is.finite(MIRZAG),],aes(x= as.factor(nbAlgo),y= MIRZAG,fill= as.factor(nbAlgo)))+theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T)+scale_fill_brewer('Blues')+scale_y_sqrt()+xlab('Number of databases where interaction is predicted')
print(p)
dev.off()

pdf(sprintf('%s/QC_miR_targets/Score_VS_numberOfMethods_miRANDA.pdf',FIG_DIR),width=3.5,height=2.5)
p <- ggplot(targets_compared[is.finite(miRANDA),],aes(x= as.factor(nbAlgo),y= miRANDA,fill= as.factor(nbAlgo)))+theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T)+scale_fill_brewer('Blues')+scale_y_sqrt()+xlab('Number of databases where interaction is predicted')
print(p)
dev.off()

pdf(sprintf('%s/QC_miR_targets/Score_VS_numberOfMethods_miRDB.pdf',FIG_DIR),width=3.5,height=2.5)
p <- ggplot(targets_compared[is.finite(miRDB),],aes(x= as.factor(nbAlgo),y= miRDB,fill= as.factor(nbAlgo)))+theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T)+scale_fill_brewer('Blues')+xlab('Number of databases where interaction is predicted')
print(p)
dev.off()

pdf(sprintf('%s/QC_miR_targets/Score_VS_numberOfMethods_targetScan.pdf',FIG_DIR),width=3.5,height=2.5)
p <- ggplot(targets_compared[is.finite(targetScan),],aes(x= as.factor(nbAlgo),y= targetScan,fill= as.factor(nbAlgo)))+theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T)+scale_fill_brewer('Blues')+xlab('Number of databases where interaction is predicted')
print(p)
dev.off()

library(UpSetR)
myTarget=as.data.frame(targets_compared[,mget(c('MIRZA','MIRZAG','miRANDA','miRDB','targetScan'))])
myTarget=as.data.frame(apply(apply(myTarget,2,is.finite),2,ifelse,1,0))
upset(myTarget,sets=colnames(myTarget))
myTarget=as.data.table(myTarget)


###############################################################################
####        Characterize miR-gene associations  (basal state)              ####
###############################################################################

##### load CAR scores
confident_targets=list()
Enrich_bymiR=list()
all_targets=list()
all_genes=list()
GeneAnnot=fread(sprintf('%s/Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/GeneAnnot_expressed.txt',EVO_IMMUNO_POP))


for (cond in 1:5){
cat(cond,'\n')
CARSCORE_cond = fread(sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/variance_explained_by_miRNA_CARScore_withoutTranscription_withoutSVA/variance_explained_by_miRNA_CARScore_%s.tsv",EVO_IMMUNO_POP,cond))
CARSCORE_cond = CARSCORE_cond[!is.na(variance_explained_by_pop)]
#CARSCORE_cond[, eval(colnames_to_delete) := NULL]

CARSCORE_cond_melt=melt(CARSCORE_cond,measure.vars=list(miRNA=paste('miRNA',1:10,sep='_'),variance=paste('miRNA_variance',1:10,sep='_'),sign=paste('miRNA_carscore_sign',1:10,sep='_')))
CARSCORE_cond_melt[,Symbol:=GeneAnnot$'Associated Gene Name'[match(CARSCORE_cond_melt$gene,GeneAnnot$Ensembl.Gene.ID)]]
CARSCORE_cond_melt=CARSCORE_cond_melt[,code_miRNA:=paste(gene, miRNA)]
CARSCORE_cond_melt=CARSCORE_cond_melt[,code_miRNA_2:=paste(Symbol, miRNA)]
CARSCORE_cond_melt=CARSCORE_cond_melt[Symbol%in%unique(targets_compared$Symbol),]
CARSCORE_cond_melt=CARSCORE_cond_melt[miRNA%in%unique(gsub("-", "_", targets_compared$miRNA)),]
mytargets=targets_compared[Symbol%in%unique(CARSCORE_cond_melt$Symbol),]
mytargets=mytargets[gsub("-", "_", miRNA)%in%unique(CARSCORE_cond_melt$miRNA),]

# Count_BS=targets_compared[,.(NbTarget=length(unique(EnsemblGeneID))),by='miRNA']

##### add miRBS informations
CARSCORE_cond_melt_withBS=merge(CARSCORE_cond_melt,unique(mytargets,by='code'),by.x=c('code_miRNA_2'),by.y=c('code'),all.x=T,suffixes=c('','.toremove'))
CARSCORE_cond_melt_withBS[,miRNA.toremove:=NULL]
CARSCORE_cond_melt_withBS[,Symbol.toremove:=NULL]
CARSCORE_cond_melt_withBS[is.na(nbAlgo),MIRZA:=-Inf]
CARSCORE_cond_melt_withBS[is.na(nbAlgo),MIRZAG:=-Inf]
CARSCORE_cond_melt_withBS[is.na(nbAlgo),miRANDA:=-Inf]
CARSCORE_cond_melt_withBS[is.na(nbAlgo),miRDB:=-Inf]
CARSCORE_cond_melt_withBS[is.na(nbAlgo),targetScan:=-Inf]
CARSCORE_cond_melt_withBS[is.na(nbAlgo),nbAlgo:=0]

all_genes[[cond]]=CARSCORE_cond_melt_withBS

all_targets[[cond]]=CARSCORE_cond_melt_withBS[nbAlgo>0,]
confident_targets[[cond]]=CARSCORE_cond_melt_withBS[nbAlgo>0 & sign<0,]

Enrich_bymiR[[cond]]=CARSCORE_cond_melt_withBS[,.(Ntarget=sum(sign<0 & nbAlgo>0),non_target=mean(sign[nbAlgo==0]<0),target1=mean(sign[nbAlgo>0]<0),target5=mean(sign[nbAlgo==5]<0),P=suppressWarnings(as.numeric(try(fisher.test(table(nbAlgo>0,sign<0))$p,silent=TRUE)))),by=.(miRNA,condition)][Ntarget>10,]
}

confident_targets=rbindlist(confident_targets)
confident_targets[,condition:=condIndex[condition]]
confident_targets[,condition:=factor(condition,condIndex)]
confident_targets[,Algo:=paste(ifelse(is.finite(miRANDA),'miRANDA',""),ifelse(is.finite(MIRZA),'MIRZA',""),ifelse(is.finite(MIRZAG),'MIRZA-G',""),ifelse(is.finite(miRDB),'miRTargets',""),ifelse(is.finite(targetScan),'targetScan',""),sep=';')]
confident_targets[,Algo:=gsub('^;','',gsub(';$','',gsub(';+',';',Algo)))]

fwrite(confident_targets,file=sprintf('%s/19_miR_targets/SupTableX_confident_targets_raw.tsv',DATA_DIR),sep='\t')
cols=c('gene','Symbol','condition','miRNA','variance','nbAlgo','Algo')
confident_targets[,mget(cols)]
# confident_targets[,mean(variance),keyby=.(condition,nbAlgo)] # pas de relation nette
confident_targets_cast=dcast(confident_targets,miRNA+gene+Symbol+nbAlgo+Algo~condition,value.var='variance',fill=NA)
confident_targets_cast[,score:=(pmax(0,NS,na.rm=T)+pmax(0,LPS,na.rm=T)+pmax(0,PAM3CSK4,na.rm=T)+pmax(0,R848,na.rm=T)+pmax(0,IAV,na.rm=T))]
fwrite(confident_targets_cast[order(-score),],file=sprintf('%s/19_miR_targets/SupTableX_confident_targets_cast.tsv',DATA_DIR),sep='\t')

Enrich_bymiR=rbindlist(Enrich_bymiR)
Enrich_bymiR[,condition:=condIndex[condition]]
fwrite(Enrich_bymiR,file=sprintf('%s/19_miR_targets/Enrichment_negCor_bymiR.tsv',DATA_DIR))

all_genes=rbindlist(all_genes)
all_genes[,condition:=condIndex[condition]]
all_genes[,condition:=factor(condition,condIndex)]
all_genes[,Algo:=paste(ifelse(is.finite(miRANDA),'miRANDA',""),ifelse(is.finite(MIRZA),'MIRZA',""),ifelse(is.finite(MIRZAG),'MIRZA-G',""),ifelse(is.finite(miRDB),'miRTargets',""),ifelse(is.finite(targetScan),'targetScan',""),sep=';')]
all_genes[,Algo:=gsub('^;','',gsub(';$','',gsub(';+',';',Algo)))]
all_genes[,variance_signed:=variance*sign]

fwrite(all_genes,file=sprintf('%s/19_miR_targets/SupTableX_all_correlated_gene_raw.tsv',DATA_DIR),sep='\t')	
fwrite(all_genes[nbAlgo>0,],file=sprintf('%s/19_miR_targets/SupTableX_all_targets_raw.tsv',DATA_DIR),sep='\t')
cols=c('gene','Symbol','condition','miRNA','variance','nbAlgo','Algo')
# confident_targets[,mean(variance),keyby=.(condition,nbAlgo)] # pas de relation nette
all_genes_cast=dcast(all_genes,miRNA+gene+Symbol+nbAlgo+Algo~condition,value.var='variance_signed',fill=NA)
all_genes_cast[,score:=(pmin(0,NS,na.rm=T)+pmin(0,LPS,na.rm=T)+pmin(0,PAM3CSK4,na.rm=T)+pmin(0,R848,na.rm=T)+pmin(0,IAV,na.rm=T))]
all_genes_cast=all_genes_cast[order(score),]
fwrite(all_genes_cast[nbAlgo>0,],file=sprintf('%s/19_miR_targets/SupTableX_all_targets_cast.tsv',DATA_DIR),sep='\t')
fwrite(all_genes_cast,file=sprintf('%s/19_miR_targets/SupTableX_all_correlated_genes_cast.tsv',DATA_DIR),sep='\t')

###############################################################################
####		        Make figures for  (basal state)              		   ####
###############################################################################

cond=5

CARSCORE_cond = fread(sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/variance_explained_by_miRNA_CARScore_withoutTranscription_withoutSVA/variance_explained_by_miRNA_CARScore_%s.tsv",EVO_IMMUNO_POP,cond))
CARSCORE_cond = CARSCORE_cond[!is.na(variance_explained_by_pop)]
#CARSCORE_cond[, eval(colnames_to_delete) := NULL]


CARSCORE_cond_melt=melt(CARSCORE_cond,measure.vars=list(miRNA=paste('miRNA',1:10,sep='_'),variance=paste('miRNA_variance',1:10,sep='_'),sign=paste('miRNA_carscore_sign',1:10,sep='_')))
CARSCORE_cond_melt[,Symbol:=GeneAnnot$'Associated Gene Name'[match(CARSCORE_cond_melt$gene,GeneAnnot$Ensembl.Gene.ID)]]
CARSCORE_cond_melt=CARSCORE_cond_melt[,code_miRNA:=paste(gene, miRNA)]
CARSCORE_cond_melt=CARSCORE_cond_melt[,code_miRNA_2:=paste(Symbol, miRNA)]
CARSCORE_cond_melt=CARSCORE_cond_melt[Symbol%in%unique(targets_compared$Symbol),]
CARSCORE_cond_melt=CARSCORE_cond_melt[miRNA%in%unique(gsub("-", "_", targets_compared$miRNA)),]
mytargets=targets_compared[Symbol%in%unique(CARSCORE_cond_melt$Symbol),]
mytargets=mytargets[gsub("-", "_", miRNA)%in%unique(CARSCORE_cond_melt$miRNA),]

# Count_BS=targets_compared[,.(NbTarget=length(unique(EnsemblGeneID))),by='miRNA']

##### add miRBS informations
CARSCORE_cond_melt_withBS=merge(CARSCORE_cond_melt,unique(mytargets,by='code'),by.x=c('code_miRNA_2'),by.y=c('code'),all.x=T,suffixes=c('','.toremove'))
CARSCORE_cond_melt_withBS[,miRNA.toremove:=NULL]
CARSCORE_cond_melt_withBS[,Symbol.toremove:=NULL]
CARSCORE_cond_melt_withBS[is.na(nbAlgo),MIRZA:=-Inf]
CARSCORE_cond_melt_withBS[is.na(nbAlgo),MIRZAG:=-Inf]
CARSCORE_cond_melt_withBS[is.na(nbAlgo),miRANDA:=-Inf]
CARSCORE_cond_melt_withBS[is.na(nbAlgo),miRDB:=-Inf]
CARSCORE_cond_melt_withBS[is.na(nbAlgo),targetScan:=-Inf]
CARSCORE_cond_melt_withBS[is.na(nbAlgo),nbAlgo:=0]

dir.create(sprintf('%s/Correlation_VS_bindingStrength/',FIG_DIR))
pdf(sprintf('%s/Correlation_VS_bindingStrength/SignedR2_VS_numberOfMethods.pdf',FIG_DIR),width=3.5,height=2.5)
library(ggplot2)
p <- ggplot(CARSCORE_cond_melt_withBS,aes(x=as.factor(nbAlgo),y=sign*variance,fill=as.factor(nbAlgo)))
p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T)
p <- p + scale_fill_brewer('Blues')+xlab('Number of databases where interaction is predicted')+ylab('Percentage of variance explained by miRNA (signed)')
print(p)	
dev.off()

pdf(sprintf('%s/Correlation_VS_bindingStrength/AbsloluteR2_VS_numberOfMethods.pdf',FIG_DIR),width=3.5,height=2.5)
library(ggplot2)
p <- ggplot(CARSCORE_cond_melt_withBS,aes(x=as.factor(nbAlgo),y=variance,fill=as.factor(nbAlgo)))
p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T)
p <- p + scale_fill_brewer('Blues')+xlab('Number of databases where interaction is predicted')+ylab('Percentage of variance explained by miRNA')
print(p)
dev.off()

pdf(sprintf('%s/Correlation_VS_bindingStrength/SignOfCorrelation_VS_numberOfMethods.pdf',FIG_DIR),width=3.5,height=2.5)
library(ggplot2)
p <- ggplot(CARSCORE_cond_melt_withBS,aes(x=as.factor(nbAlgo),fill=as.factor(sign)))
p <- p + theme_classic()+geom_bar(stat='count',position="fill")
p <- p + scale_fill_manual(values=c("#1F78B4AA","#E31A1CAA"))+xlab('Number of databases where an interaction is predicted')+ylab('Percentage of miRNA-gene correlations')	
print(p)
dev.off()

##### increasing thresholds
CARSCORE_cond_melt_withBS=all_genes[condition=='NS',] #  3578 miR-gene pairs (genes and miRNAs that are not present in all databses were removed)

CARSCORE_NS_target=list(none=all_genes[condition=='NS' & nbAlgo==0,],
						miRANDA=all_genes[condition=='NS' & grepl('miRANDA',Algo) & nbAlgo>=1,],
						miRANDA_plus_1=all_genes[condition=='NS' & grepl('miRANDA',Algo) & nbAlgo>=2,],
						miRANDA_plus_2=all_genes[condition=='NS' & grepl('miRANDA',Algo) & nbAlgo>=3,],
						miRANDA_plus_3=all_genes[condition=='NS' & grepl('miRANDA',Algo) & nbAlgo>=4,],
						miRANDA_plus_4=all_genes[condition=='NS' & grepl('miRANDA',Algo) & nbAlgo>=5,])
for (i in 0:5){
	CARSCORE_NS_target[[i+1]][,Detection:=names(CARSCORE_NS_target)[i+1]]
}
CARSCORE_NS_target=rbindlist(CARSCORE_NS_target)
CARSCORE_NS_target[,Detection:=factor(Detection,c('none','miRANDA','miRANDA_plus_1','miRANDA_plus_2','miRANDA_plus_3','miRANDA_plus_4'))]
dir.create(sprintf('%s/Correlation_VS_bindingStrength/',FIG_DIR))
pdf(sprintf('%s/Correlation_VS_bindingStrength/AbsloluteR2_VS_numberOfMethods_cumulative.pdf',FIG_DIR),width=2.5,height=4)
library(ggplot2)
p <- ggplot(CARSCORE_NS_target,aes(x=as.factor(Detection),y=variance,fill=as.factor(Detection)))
p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T)
p <- p + scale_fill_brewer('Blues')+xlab('Databases where interaction is predicted')+ylab('Percentage of variance explained by miRNA')
p <- p + theme(legend.position="top")+ theme(legend.position="top", axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

pdf(sprintf('%s/Correlation_VS_bindingStrength/SignOfCorrelation_VS_numberOfMethods_cumulative.pdf',FIG_DIR),width=2.5,height=4)
library(ggplot2)
p <- ggplot(CARSCORE_NS_target,aes(x=as.factor(Detection),fill=as.factor(sign)))
p <- p + theme_classic()+geom_bar(stat='count',position="fill")
p <- p + scale_fill_manual(values=c("#1F78B4AA","#E31A1CAA"))+xlab('Databases where an interaction is predicted')+ylab('Percentage of miRNA-gene correlations')	
p <- p + theme(legend.position="top")+ theme(legend.position="top", axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
	dev.off()

for( cond in condIndex){
CARSCORE_cond_melt_withBS=all_genes[condition==cond,] #  3578 miR-gene pairs (genes and miRNAs that are not present in all databses were removed)

CARSCORE_cond_target=list(none=all_genes[condition==cond & nbAlgo==0,],
						miRANDA=all_genes[condition==cond & grepl('miRANDA',Algo) & nbAlgo>=1,],
						miRANDA_plus_1=all_genes[condition==cond & grepl('miRANDA',Algo) & nbAlgo>=2,],
						miRANDA_plus_2=all_genes[condition==cond & grepl('miRANDA',Algo) & nbAlgo>=3,],
						miRANDA_plus_3=all_genes[condition==cond & grepl('miRANDA',Algo) & nbAlgo>=4,],
						miRANDA_plus_4=all_genes[condition==cond & grepl('miRANDA',Algo) & nbAlgo>=5,])
for (i in 0:5){
	CARSCORE_cond_target[[i+1]][,Detection:=names(CARSCORE_cond_target)[i+1]]
}
CARSCORE_cond_target=rbindlist(CARSCORE_cond_target)
CARSCORE_cond_target[,Detection:=factor(Detection,c('none','miRANDA','miRANDA_plus_1','miRANDA_plus_2','miRANDA_plus_3','miRANDA_plus_4'))]


pdf(sprintf('%s/Correlation_VS_bindingStrength/AbsloluteR2_VS_numberOfMethods_cumulative_%s.pdf',FIG_DIR,cond),width=2.5,height=4)
library(ggplot2)
p <- ggplot(CARSCORE_cond_target,aes(x=as.factor(Detection),y=variance,fill=as.factor(Detection)))
p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T)
p <- p + scale_fill_brewer('Blues')+xlab('Databases where interaction is predicted')+ylab('Percentage of variance explained by miRNA')
p <- p + theme(legend.position="top")+ theme(legend.position="top", axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

pdf(sprintf('%s/Correlation_VS_bindingStrength/SignOfCorrelation_VS_numberOfMethods_cumulative_%s.pdf',FIG_DIR,cond),width=2.5,height=4)
library(ggplot2)
p <- ggplot(CARSCORE_cond_target,aes(x=as.factor(Detection),fill=as.factor(sign)))
p <- p + theme_classic()+geom_bar(stat='count',position="fill")
p <- p + scale_fill_manual(values=c("#1F78B4AA","#E31A1CAA"))+xlab('Databases where an interaction is predicted')+ylab('Percentage of miRNA-gene correlations')	
p <- p + theme(legend.position="top")+ theme(legend.position="top", axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()
}







theshold1=c(.2,.2,.6,.6,.6)
theshold2=c(.4,.4,.8,.8,.8)
names(theshold1)=c('MIRZA','MIRZAG','miRANDA','miRDB','targetScan')
names(theshold2)=c('MIRZA','MIRZAG','miRANDA','miRDB','targetScan')

for (method in c('MIRZA','MIRZAG','miRANDA','miRDB','targetScan')){
	pdf(sprintf('%s/Correlation_VS_bindingStrength/SignedR2_VS_%s.pdf',FIG_DIR,method),width=3.5,height=2.5)
	library(ggplot2)
	binding=CARSCORE_cond_melt_withBS[,get(method)]
	CARSCORE_cond_melt_withBS[,binding:=NULL]
	CARSCORE_cond_melt_withBS[,binding:=sapply(1 + is.finite(binding)+ (binding>theshold1[method])+ (binding>theshold2[method]),switch,'no','weak','moderate','strong')]
	CARSCORE_cond_melt_withBS[,binding:=factor(binding,levels=c('no','weak','moderate','strong'))]
	p <- ggplot(CARSCORE_cond_melt_withBS,aes(x=binding,y=sign*variance,fill=binding))
	p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
	p <- p + scale_fill_brewer('Blues')+xlab('Number of databases where interaction is predicted')+ylab('Percentage of variance explained by miRNA (signed)')
	print(p)
	dev.off()

	pdf(sprintf('%s/Correlation_VS_bindingStrength/AbsloluteR2_VS_%s.pdf',FIG_DIR,method),width=3.5,height=2.5)
	library(ggplot2)
	p <- ggplot(CARSCORE_cond_melt_withBS,aes(x=binding,y=variance,fill=binding))
	p <- p + theme_classic()+geom_violin(scale='width')+geom_boxplot(width=0.5,fill='#FFFFFF88',notch=T)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
	p <- p + scale_fill_brewer('Blues')+xlab('Number of databases where interaction is predicted')+ylab('Percentage of variance explained by miRNA')
	print(p)
	dev.off()

	pdf(sprintf('%s/Correlation_VS_bindingStrength/SignOfCorrelation_VS_%s.pdf',FIG_DIR,method),width=3.5,height=2.5)
	library(ggplot2)
	p <- ggplot(CARSCORE_cond_melt_withBS,aes(x=binding,fill=as.factor(sign)))
	p <- p + theme_classic()+geom_bar(stat='count',position="fill") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
	p <- p + scale_fill_manual(values=c("#1F78B4AA","#E31A1CAA"))+xlab('Number of databases where an interaction is predicted')+ylab('Percentage of miRNA-gene correlations')	
	print(p)
	dev.off()
}

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

###########################################################
###		identify miRNA pair that share targets			###
###########################################################


#### search for miRNA that have binding sites with 20 bp more often than expected by chance.
cols=c('transcripts','EnsemblGeneID' ,'GeneNames','miRNA','miRNABS_start_pos','miRNABS_length')
miR_pair=list()
Dist=20
i=0
for(GENE in allGenes){
	i=i+1
	if(i%%20==0){cat(i,'\n')}
	cat('.')
	miR_pair[[GENE]]=merge(target_miRANDA[EnsemblGeneID==GENE,mget(cols)],target_miRANDA[EnsemblGeneID==GENE,mget(cols)],by=c('transcripts','EnsemblGeneID' ,'GeneNames'),allow.cartesian=TRUE)
	miR_pair[[GENE]]=miR_pair[[GENE]][miRNA.x!=miRNA.y | miRNABS_start_pos.x != miRNABS_start_pos.y,]
	miR_pair[[GENE]]=miR_pair[[GENE]][(miRNABS_start_pos.x + miRNABS_length.x) < miRNABS_start_pos.y &  miRNABS_start_pos.y < (miRNABS_start_pos.x + miRNABS_length.x + Dist),]
}

miR_pair=rbindlist(miR_pair)
miR_pair[,NbColoc:=max(table(transcripts)),by=.(miRNA.x,miRNA.y,EnsemblGeneID,GeneNames)]
miR_pair_gene=miR_pair[, .(NbColoc_gene=sum(NbColoc>0)),by=.(miRNA.x,miRNA.y)]
miR_pair_gene = merge(merge(miR_pair_gene, miR_count_gene,by.x='miRNA.x',by.y='miRNA'), miR_count_gene,by.x='miRNA.y',by.y='miRNA')


# get number of targets per miR in MIRANDA
p <- ggplot(miR_count_gene,aes(x=log2(1+NbTarg))) + geom_histogram()
p <- p + scale_color_brewer(palette="Greens")+ scale_fill_brewer(palette="Greens") + theme_classic() 
p <- p + ylab('Number of microRNAs')+ xlab('log2(Number of targets +1)')
pdf(sprintf('%s/QC_miR_targets/NumberOtTargets_MIRANDA.pdf',FIG_DIR),width=3.5,height=2.5)
print(p)
dev.off()


####### obtain average size of 3'UTRs
library(seqinr)
UTRseq=read.fasta(sprintf('%s/Martin/miRNA/11.miRNA_binding_sites/temp_data/temp_transcript.fasta',EVO_IMMUNO_POP))
mean(sapply(UTRseq,length)) # 917.9 bp

miR_pair_gene[,expected_Coloc:= NbTarg.x* NbTarg.y/length(allGenes)*(Dist)/mean(sapply(UTRseq,length))]
miR_pair_gene[,P_Coloc_enrich:= ppois(NbColoc_gene,expected_Coloc,low=F)]
miR_pair_gene[,Padj_Coloc_enrich:= p.adjust(P_Coloc_enrich,'fdr')]

samp=sample(1:nrow(miR_pair_gene),10000);
dir.create(sprintf('%s/mir_pairs_colocalization',FIG_DIR))
### relation between expected and observed number of colocalized miRBS.
pdf(sprintf('%s/mir_pairs_colocalization/Enrichment_colocalized_miRNA_targets.pdf',FIG_DIR))
plot(log2(miR_pair_gene[samp, expected_Coloc]),log2(miR_pair_gene[samp, NbColoc_gene]),pch=15,col=ifelse(miR_pair_gene[samp,Padj_Coloc_enrich<.01],'#E31A1C44','#52525244'),cex=0.5,xlab='log2 of expected number of genes with colocalized targets',ylab='log2 of number of genes with colocalized targets')
abline(0,1,col='red')
dev.off()
### QQplot
pdf(sprintf('%s/mir_pairs_colocalization/Enrichment_colocalized_miRNA_targets_qqplot.pdf',FIG_DIR))
plot(sort(log2(miR_pair_gene[samp, expected_Coloc])),sort(log2(miR_pair_gene[samp, NbColoc_gene])),pch=15,col=ifelse(rev(sort(miR_pair_gene[samp,Padj_Coloc_enrich]))<0.01,'#E31A1C44','#52525244'),cex=0.5,xlab='log2 of expected number of genes with colocalized targets',ylab='log2 of number of genes with colocalized targets')
abline(0,1,col='red')
dev.off()
dir.create(sprintf('%s/19_miR_targets/mir_pairs_colocalization',DATA_DIR))
fwrite(miR_pair_gene,file=sprintf('%s/19_miR_targets/mir_pairs_colocalization/miR_pairs_genecount.txt',DATA_DIR),sep='\t')
fwrite(miR_pair,file=sprintf('%s/19_miR_targets/mir_pairs_colocalization/miR_pairs.txt',DATA_DIR),sep='\t')


