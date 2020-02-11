

#################################################################################
######################        load libraries            #########################
#################################################################################

require(ggplot2)

########################################################################################
######################             load the data               #########################
########################################################################################

eQTL_bestSNP=fread(sprintf("%s/Maxime/evo_immuno_pop_QTLs/eQTL/eQTL_assoc_bestSNP.txt",EVO_IMMUNO_POP))
ExonAnnot=fread(sprintf("%s/Maxime/Evo_Immuno_pop_data/ExonCoordinates_hg37_ens70.txt",EVO_IMMUNO_POP))
colnames(ExonAnnot)=make.names(colnames(ExonAnnot))
GeneAnnot=fread(sprintf("%s/Maxime/Evo_Immuno_pop_data/GeneAnnotation_hg37_ens70.txt",EVO_IMMUNO_POP))
colnames(GeneAnnot)=make.names(colnames(GeneAnnot))

sampleAnnot=fread(sprintf("%s/Maxime/Evo_Immuno_pop_data/SampleAnnot.txt",EVO_IMMUNO_POP))
nInd_byCond=table(sampleAnnot$condition)


miRqtl_cis_annot= fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4A_mirQTLs_Annoted.tsv",EVO_IMMUNO_POP))
miRqtl_cis_annot[,conserved_TSS:=Conservation>.2]
miRqtl_cis_annot[,R2:=pmax(R2_NS,R2_LPS,R2_PAM3CSK4,R2_R848,R2_IAV,na.rm=T)]

mir_TSS_annot = fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_TSS_annot_DeRie2017.txt',EVO_IMMUNO_POP))
mir_TSS_annot[,has_mirQTL:=hsa_ID%in%miRqtl_cis_annot$miRNA]
mir_TSS_annot[,conserved_TSS:=Conservation>.2]

#########################################################
##################### Annotate SNPS #####################
#########################################################

##--------------------------------------------##
## private functions to query mySQL databases ##
##--------------------------------------------##

getSNP<-function(rs,pop='.*'){
	SNP=getGenos(rs)[-(1:5)]
	ind=colnames(SNP)
	SNP=as.numeric(SNP)
	names(SNP)=ind
	SNP[grep(pop,names(SNP))]
    }

getGenos=function(snpList,path=EVO_IMMUNO_POP,phased=FALSE,chrlist=NULL){
	require(dplyr)
	mydb=src_sqlite(paste(path,'/Maxime/SNP_annotation/imputed_Genotypes/GenoDataBase.db',sep=''))
#	Map=tbl(mydb,'Map_imputed')
#	map=filter(Map,translate_sql(snp.name%in%escape(ident(snpList)))) %>% collect()
	if(is.null(chrlist)){
		map=tbl(mydb,sql(paste('SELECT * FROM Map_imputed WHERE "snp.name" in (',paste(paste('"',snpList,'"',sep=''),collapse=', '),')',sep=''))) %>%collect()
#		map=tbl(mydb,sql(paste('SELECT * FROM Map_imputed WHERE "snp.name" in',escape(ident(snpList)))))
		chrlist=unique(map$chromosome)
	}
	RES=NULL
	for(CHR in chrlist){
#		Geno=tbl(mydb,paste('Geno',CHR,sep='_'))
#		res=filter(Geno,snp.name%in%escape(ident(snpList))) %>% collect()
		res=tbl(mydb,sql(paste('SELECT * FROM ',ifelse(phased,'Haplo','Geno'),'_',CHR,' WHERE "snp.name" in (',paste(paste('"',snpList,'"',sep=''),collapse=', '),')',sep=''))) %>%collect()
		RES=rbind(RES,res)
		}
	RES=RES[match(snpList,RES$snp.name),]
	as.data.frame(RES)
}

getMapInfo <- function(snpList,path=EVO_IMMUNO_POP){
	require(dplyr)
	mydb=src_sqlite(paste(path,'/Maxime/SNP_annotation/imputed_Genotypes/GenoDataBase_v9.db',sep=''))
	map=tbl(mydb,sql(paste('SELECT * FROM Map_imputed WHERE "snp.name" in (',paste(paste('"',snpList,'"',sep=''),collapse=', '),')',sep=''))) %>%collect()
as.data.frame(map)
}

getMapInfo_locus<-function(chr,start, end,path=EVO_IMMUNO_POP){
	require(dplyr)
	mydb=src_sqlite(paste(path,'/Maxime/SNP_annotation/imputed_Genotypes/GenoDataBase_v9.db',sep=''))
	map=tbl(mydb,sql(paste('SELECT * FROM Map_imputed WHERE "chromosome" in ("',chr,'") AND "position" >',start,' AND "position" <',end ,sep=''))) %>%collect(n=Inf)
    as.data.frame(map)
    }

getRecombinationRate_locus=function(chr,start, end,path=EVO_IMMUNO_POP,pop=c('CEU','YRI')){
	pop=match.arg(pop)
	require(dplyr)
	mydb=src_sqlite(paste(path,'/Maxime/SNP_annotation/imputed_Genotypes/GenoDataBase_v6.db',sep=''))
	map=tbl(mydb,sql(paste('SELECT * FROM RecombRate WHERE "chromosome" in ("',chr,'") AND "population" in ("',pop,'") AND "Position.bp." >',start,' AND "Position.bp." <',end ,sep=''))) %>%collect(n=Inf)
as.data.frame(map)
}

##--------------------------------------------##
## private version (based on mysql database)  ##
##--------------------------------------------##

snps_information=getMapInfo(unique(eQTL_bestSNP$snps))[c("snp.name","chromosome","position","allele.2","allele.1","SNPfreq_AF","SNPfreq_EU","daf_char_AFB","daf_char_EUB","RegElt","TFBS", 'FST_adj','iHS_AFB','iHS_EUB','P_CLS2_EUB','P_CLS2_AFB','aSNP')]
fwrite(snps_information,file=sprintf("%s/Maxime/miRNA_V2/data/15.eQTL_comparisons/SupTableXX_eQTLs_snp_Annot.tsv",EVO_IMMUNO_POP),sep='\t')

##--------------------------------------------##
## private version (based on mysql database)  ##
##--------------------------------------------##
snps_information=fread(file=sprintf("%s/Maxime/miRNA_V2/data/15.eQTL_comparisons/SupTableXX_eQTLs_snp_Annot.tsv",EVO_IMMUNO_POP),sep='\t')

eQTL_cis_annot=merge(eQTL_bestSNP,snps_information,by.x='snps',by.y='snp.name',all.x=T)

#Exon_bestSNP=merge(eQTL_cis_annot,ExonAnnot,by.x='gene',by.y='Ensembl.Gene.ID',all.x=T)

distToRange=function(pos,start,end,strand=NULL,chr=NULL,chrRange=NULL){
        abs_distance_if_outside=pmin(abs(pos-start),abs(pos-end))
        is_before=(pos<start)
        is_after=(end<pos)
        is_inside=start<=pos & pos<=end
        Dist=ifelse( is_inside, 0, abs_distance_if_outside)
        if(!is.null(strand)){
            if(is.numeric(strand)){
                strand=ifelse(strand>0,'+','-')
            }
            Dist = Dist * ifelse(is_before & strand=='+' | is_after & strand=='-', -1,1)
        }
        if(!is.null(chr) & !is.null(chrRange)){
            Dist[chr!=chrRange]=Inf
        }
        Dist
}
#Exon_bestSNP[,CisDist_Gene:=distToRange(position, Gene.Start..bp., Gene.End..bp., Strand)/1e3 ]
#Exon_bestSNP[,CisDist_TSS:=distToRange(position, Transcript.Start..bp., Transcript.Start..bp., Strand)/1e3 ]
#Exon_bestSNP[,CisDist_TES:=distToRange(position, Transcript.End..bp., Transcript.End..bp., Strand)/1e3 ]
#Exon_bestSNP[,CisDist_TSS_gene:=distToRange(position, Gene.Start..bp., Gene.Start..bp., Strand)/1e3 ]
#Exon_bestSNP[,CisDist_TES_gene:=distToRange(position, Gene.End..bp., Gene.End..bp., Strand)/1e3 ]
#Exon_bestSNP[,CisDist_Exon:=distToRange(position, Exon.Chr.Start..bp., Exon.Chr.End..bp., Strand)/1e3 ]
#Exon_bestSNP[,CisDist_CDS:=distToRange(position, Genomic.coding.start, Genomic.coding.end, Strand)/1e3 ]
#Exon_bestSNP[,CisDist_5UTR:=distToRange(position, X5..UTR.Start, X5..UTR.End, Strand)/1e3 ]
#Exon_bestSNP[,CisDist_3UTR:=distToRange(position, X3..UTR.Start, X3..UTR.End, Strand)/1e3 ]
#
#CisDistGene=unique(Exon_bestSNP[order(gene,abs(CisDist_Gene)),.(CisDist_Gene,gene)],by="gene")
#CisDistExon=unique(Exon_bestSNP[order(gene,abs(CisDist_Exon)),.(CisDist_Exon,gene)],by="gene")
#CisDistCDS=unique(Exon_bestSNP[order(gene,abs(CisDist_CDS)),.(CisDist_CDS,gene)],by="gene")
#CisDist3UTR=unique(Exon_bestSNP[order(gene,abs(CisDist_3UTR)),.(CisDist_3UTR,gene)],by="gene")
#CisDist5UTR=unique(Exon_bestSNP[order(gene,abs(CisDist_5UTR)),.(CisDist_5UTR,gene)],by="gene")
#CisDistTSS=unique(Exon_bestSNP[order(gene,abs(CisDist_TSS)),.(CisDist_TSS,gene)],by="gene")
#CisDistTES=unique(Exon_bestSNP[order(gene,abs(CisDist_TES)),.(CisDist_TES,gene)],by="gene")
#CisDistTSS_gene=unique(Exon_bestSNP[order(gene,abs(CisDist_TSS_gene)),.(CisDist_TSS_gene,gene)],by="gene")
#CisDistTES_gene=unique(Exon_bestSNP[order(gene,abs(CisDist_TES_gene)),.(CisDist_TES_gene,gene)],by="gene")
#
#eQTL_cis_annot=merge(eQTL_cis_annot,CisDistGene,all.x=T,by='gene')
#eQTL_cis_annot=merge(eQTL_cis_annot,CisDistExon,all.x=T,by='gene')
#eQTL_cis_annot=merge(eQTL_cis_annot,CisDistCDS,all.x=T,by='gene')
#eQTL_cis_annot=merge(eQTL_cis_annot,CisDist3UTR,all.x=T,by='gene')
#eQTL_cis_annot=merge(eQTL_cis_annot,CisDist5UTR,all.x=T,by='gene')
#eQTL_cis_annot=merge(eQTL_cis_annot,CisDistTSS,all.x=T,by='gene')
#eQTL_cis_annot=merge(eQTL_cis_annot,CisDistTES,all.x=T,by='gene')
#eQTL_cis_annot=merge(eQTL_cis_annot,CisDistTSS_gene,all.x=T,by='gene')
#eQTL_cis_annot=merge(eQTL_cis_annot,CisDistTES_gene,all.x=T,by='gene')

enhancers_monocytes = fread(paste(EVO_IMMUNO_POP, "Martin/Project_Neanderthal/ChromHMMFunctional/RegionDefinitions/Regions/Enhancers/enhancers_", "E029", ".tsv", sep=""))
enhancers_monocytes_GR = makeGRangesFromDataFrame(enhancers_monocytes)
enhancers_monocytes_GR= reduce(enhancers_monocytes_GR)
seqlevelsStyle(enhancers_monocytes_GR) <- "NCBI"

promoters_monocytes = fread(paste(EVO_IMMUNO_POP, "Martin/Project_Neanderthal/ChromHMMFunctional/RegionDefinitions/Regions/Promoters/promoters_", "E029", ".tsv", sep=""))
promoters_monocytes_GR = makeGRangesFromDataFrame(promoters_monocytes)
promoters_monocytes_GR=reduce(promoters_monocytes_GR)
seqlevelsStyle(promoters_monocytes_GR) <- "NCBI"
#
#PIRs = fread(paste(EVO_IMMUNO_POP, "Martin/Ressources/Javierreetal/PCHiC_peak_matrix_cutoff5.tsv", sep=""))
#PIRs = PIRs[Mon>5]
#PIRs_GR = GRanges(paste("chr", PIRs$oeChr, sep=""),
#                  IRanges(start = PIRs$oeStart, end = PIRs$oeEnd))
#seqlevelsStyle(PIRs_GR) <- "NCBI"
#
#bait_GR = GRanges(paste("chr", PIRs$baitChr, sep=""),
#                  IRanges(start = PIRs$baitStart, end = PIRs$baitEnd))
#seqlevelsStyle(bait_GR) <- "NCBI"

#frequent_snps = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/08.snp_data/general_informations/all_snps_MAFover5.tsv", sep=""))
#frequent_snps[, snp_name := snp.name]
#frequent_snps[, snp.name := NULL]

snps_data = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/08.snp_data/general_informations/all_snps.tsv", sep=""))
snp_data=fread(paste(EVO_IMMUNO_POP, "/Maxime/Evo_Immuno_pop_data/SNP_annotations/Map_imputed_essential_informations.txt", sep=""))
snps_data=merge(snps_data,snp_data)
frequent_snps=snps_data[(MAF_EUB>.05 | MAF_AFB>.05),]
frequent_snps[, snp_name := snp.name]
frequent_snps[, snp.name := NULL]

frequent_snps_GR=makeGRangesFromDataFrame(frequent_snps,start.field='position',end.field='position')

#frequent_snps[,RegElt:='']
#oo=findOverlaps(frequent_snps_GR,promoters_monocytes_GR[width(promoters_monocytes_GR)>1000])
#frequent_snps[unique(queryHits(oo)),RegElt:='Promoter (long)']
#oo=findOverlaps(frequent_snps_GR,promoters_monocytes_GR[width(promoters_monocytes_GR)<=1000])
#frequent_snps[unique(queryHits(oo)),RegElt:='Promoter (short)']
#oo=findOverlaps(frequent_snps_GR,enhancers_monocytes_GR[width(enhancers_monocytes_GR)>1000])
#frequent_snps[unique(queryHits(oo)),RegElt:='Enhancer (long)']
#oo=findOverlaps(frequent_snps_GR,enhancers_monocytes_GR[width(enhancers_monocytes_GR)<=1000])
#frequent_snps[unique(queryHits(oo)),RegElt:='Enhancer (short)']
#
#oo=findOverlaps(frequent_snps_GR,PIRs_GR[width(enhancers_monocytes_GR)<=1000])
#frequent_snps[unique(queryHits(oo)),RegElt:=paste(RegElt, ifelse(RegElt=='','PIR',' - PIR'),sep='')]

frequent_snps[,RegElt_HMM:='']
oo=findOverlaps(frequent_snps_GR,promoters_monocytes_GR)
frequent_snps[unique(queryHits(oo)),RegElt_HMM:='Promoter']
oo=findOverlaps(frequent_snps_GR,enhancers_monocytes_GR)
frequent_snps[unique(queryHits(oo)),RegElt_HMM:='Enhancer']

frequent_snps[,is_EQTL:=(snp_name%in%eQTL_cis_annot[FDR<.05,snps])]
frequent_snps[,is_miRQTL:=(snp_name%in%miRqtl_cis_annot$snps)]


############################################
##      Enrichments Prom_Enhancer_PIR     ##
############################################
colorSet=c(brewer.pal(8,'Set2'),brewer.pal(8,'Pastel2'),brewer.pal(8,'Dark2'),brewer.pal(8,'Set1'),brewer.pal(8,'Pastel1'))
colorSet_ordered=colorSet[c(rep(c(1,0,2)*8,8)+rep(1:8,e=3),24+rep((1:0)*8,8)+rep(1:8,e=2))]]
colReg = c("#CCCCCC", "#B3E2CD", "#8DA0CB")
names(colReg)= c("other","Enhancer","Promoter")


####### make pie charts 

pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/eQTL_miRQTL_Prom_Enhancer_observed_simple.pdf',EVO_IMMUNO_POP),width=5,height=5)
oldmar=par()$mar
layout(matrix(1:4,2))
par(mar=c(0.1,2,2,.1))
tab0=table(frequent_snps[,RegElt_HMM])
pie(tab0,col=colReg,main='all SNPs, MAF>5%',labels=NA)
names(tab0)[1]='other'
tab1=table(frequent_snps[is_miRQTL==TRUE,RegElt_HMM])
pie(tab1,col=colReg,main='miRQTL',labels=NA)
names(tab1)[1]='other'
tab2=table(frequent_snps[is_EQTL==TRUE,RegElt_HMM])
pie(tab2,col=colReg,main='eQTL',labels=NA)
names(tab2)[1]='other'
plot.new()
legend('top',bty='n',fill=colReg[-1],legend=names(colReg)[-1], ncol=1)
par(mar=oldmar)
dev.off()

####### compute Odds ratios 

OR_miRQTL=do.call(rbind,lapply(names(colReg),function(x){
                        if(!x %in% names(tab1)) 
                             oframe <- data.frame(Category=x,LowerCI = NA, OR = NA, UpperCI = NA, P=NA)
                        else{
                            myTable=matrix(c(tab1[x],sum(tab1)-tab1[x],tab0[x],sum(tab0)-tab0[x]),2,byrow=T)
                            test=fisher.test(myTable)
                            oframe <- data.frame(Category=x,LowerCI = test$conf.int[1], OR = test$est, UpperCI = test$conf.int[2], P=test$p.value)
                        }
                        oframe}))

OR_eQTL=do.call(rbind,lapply(names(colReg),function(x){
                        myTable=matrix(c(tab2[x],sum(tab2)-tab2[x],tab0[x],sum(tab0)-tab0[x]),2,byrow=T)
                        test=fisher.test(myTable)
                        oframe <- data.frame(Category=x,LowerCI = test$conf.int[1], OR = test$est, UpperCI = test$conf.int[2], P=test$p.value)
                        oframe}))

fwrite(OR_eQTL,file=sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/eQTL_RegElt_Enrichment_simple.txt',EVO_IMMUNO_POP),sep='\t')
fwrite(OR_miRQTL,file=sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/miRQTL_RegElt_Enrichment_simple.txt',EVO_IMMUNO_POP),sep='\t')

OR_eQTL
#            Category    LowerCI         OR    UpperCI             P
#odds ratio     other  0.1014512  0.1087526  0.1166535  0.000000e+00
#odds ratio1 Enhancer  4.8816712  5.3087564  5.7664554 1.202183e-242
#odds ratio2 Promoter 17.9893322 19.8190921 21.8046911  0.000000e+00

OR_miRQTL
#            Category    LowerCI         OR   UpperCI            P
#odds ratio     other 0.09704344  0.1513752  0.242364 5.654486e-13
#odds ratio1 Enhancer 1.94337621  3.6157152  6.284897 5.921567e-05
#odds ratio2 Promoter 8.86902716 17.2887085 31.061450 4.562137e-12


tab=cbind(table(grepl('Enhancer',frequent_snps[is_miRQTL==TRUE,RegElt_HMM])),table(grepl('Enhancer',frequent_snps[is_EQTL==TRUE,RegElt_HMM])));odds.ratio(tab)
#             LowerCI       OR  UpperCI alpha         P
#odds ratio 0.7909155 1.415138 2.719653  0.05 0.2933827
tab=cbind(table(grepl('Prom',frequent_snps[is_miRQTL==TRUE,RegElt_HMM])),table(grepl('Prom',frequent_snps[is_EQTL==TRUE,RegElt_HMM])));odds.ratio(tab)
#             LowerCI       OR  UpperCI alpha         P
#odds ratio 0.6037184 1.122042 2.275623  0.05 0.8804294



########### Compute percentage of genes with an eQTL
load(paste(HOME,'/Annotation/GOterms/GOterms_12578_genes.Rdata',sep=''))
TF_genes=allGOterms$gene[allGOterms$go=='GO:0003700']
non_codingRNAs=GeneAnnot$Ensembl.Gene.ID[GeneAnnot$Gene.Biotype!='protein_coding' & GeneAnnot$Expressed]
all_genes=GeneAnnot$Ensembl.Gene.ID[GeneAnnot$Expressed]

GeneConservation=fread('/Volumes/@Home/Annotation/Conservation/GeneConservation_hg37_ens70_V4.txt')
LOF_intolerant=intersect(GeneConservation[pLI>.9,gene_id],all_genes)
LOF_Recessive=intersect(GeneConservation[pRec>.9,gene_id],all_genes)
LOF_Null=intersect(GeneConservation[pNull>.9,gene_id],all_genes)

percentage_all_genes_with_eQTL=replicate(1000,mean(sample(all_genes%in%eQTL_cis_annot[FDR<.05,gene],replace=T)))
percentage_all_miR_with_eQTL=replicate(1000,mean(sample(mir_TSS_annot[!duplicated(hsa_ID), has_mirQTL],replace=T)))
percentage_all_TF_with_eQTL=replicate(1000,mean(sample(TF_genes%in%eQTL_cis_annot[FDR<.05,gene],replace=T)))
percentage_all_lncRNA_with_eQTL=replicate(1000,mean(sample(non_codingRNAs%in%eQTL_cis_annot[FDR<.05,gene],replace=T)))
percentage_all_protein_coding_with_eQTL=replicate(1000,mean(sample(setdiff(all_genes,non_codingRNAs)%in%eQTL_cis_annot[FDR<.05,gene],replace=T)))
percentage_all_LOF_intolerant_with_eQTL=replicate(1000,mean(sample(LOF_intolerant%in%eQTL_cis_annot[FDR<.05,gene],replace=T)))
percentage_all_LOF_Recessive_with_eQTL=replicate(1000,mean(sample(LOF_Recessive%in%eQTL_cis_annot[FDR<.05,gene],replace=T)))
percentage_all_LOF_Null_with_eQTL=replicate(1000,mean(sample(LOF_Null%in%eQTL_cis_annot[FDR<.05,gene],replace=T)))

DT=data.table(class=rep(c('All','protein-coding','LOF Intolerant','LOF Recessive','LOF Null','lncRNA','TF','miRNA'),each=1000),Pct_with_eQTL=c(percentage_all_genes_with_eQTL,percentage_all_protein_coding_with_eQTL,percentage_all_LOF_intolerant_with_eQTL,percentage_all_LOF_Recessive_with_eQTL,percentage_all_LOF_Null_with_eQTL,percentage_all_lncRNA_with_eQTL,percentage_all_TF_with_eQTL,percentage_all_miR_with_eQTL))
fwrite(DT,file=sprintf("%s/Maxime/miRNA_V2/data/15.eQTL_comparisons/SourceData_Pct_eQTLs_byGeneClass.tsv",EVO_IMMUNO_POP),sep='\t')
#Olfactory_genes=allGOterms$gene[allGOterms$go=='GO:0004984']
#percentage_all_Olfactory=replicate(1000,mean(sample(Olfactory_genes%in%eQTL_cis_annot[FDR<.05,gene],replace=T)))

eQTL_cis_annot[, R2:=statistic^2/(nInd_byCond[condition]+statistic^2)]
eQTL_R2_all_genes = eQTL_cis_annot[FDR<.05 & gene %in% all_genes,R2]
eQTL_R2_miR = miRqtl_cis_annot$R2[miRqtl_cis_annot$chrom!=14]
eQTL_R2_TF = eQTL_cis_annot[FDR<.05 & gene %in% TF_genes,R2]
eQTL_R2_lncRNA = eQTL_cis_annot[FDR<.05 & gene %in% non_codingRNAs,R2]
eQTL_R2_protein_coding = eQTL_cis_annot[FDR<.05 & !gene %in% non_codingRNAs,R2]
eQTL_R2_LOF_intolerant = eQTL_cis_annot[FDR<.05 & gene %in% LOF_intolerant,R2]
eQTL_R2_LOF_Recessive = eQTL_cis_annot[FDR<.05 & gene %in% LOF_Recessive,R2]
eQTL_R2_LOF_Null = eQTL_cis_annot[FDR<.05 & gene %in% LOF_Null,R2]
DT=data.table(class=rep(c('All','protein-coding','LOF Intolerant','LOF Recessive','LOF Null','lncRNA','TF','miRNA'),c(length(eQTL_R2_all_genes),length(eQTL_R2_protein_coding),length(eQTL_R2_LOF_intolerant),length(eQTL_R2_LOF_Recessive),length(eQTL_R2_LOF_Null),length(eQTL_R2_lncRNA),length(eQTL_R2_TF),length(eQTL_R2_miR))),
              eQTL_R2=c(eQTL_R2_all_genes,eQTL_R2_protein_coding,eQTL_R2_LOF_intolerant,eQTL_R2_LOF_Recessive,eQTL_R2_LOF_Null,eQTL_R2_lncRNA,eQTL_R2_TF,eQTL_R2_miR))
fwrite(DT,file=sprintf("%s/Maxime/miRNA_V2/data/15.eQTL_comparisons/SourceData_eQTL_R2_byGeneClass.tsv",EVO_IMMUNO_POP),sep='\t')


#### make the plot (needs improvement)
DT=fread(sprintf("%s/Maxime/miRNA_V2/data/15.eQTL_comparisons/SourceData_Pct_eQTLs_byGeneClass.tsv",EVO_IMMUNO_POP),sep='\t')
pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/Pct_eQTLs_byGeneClass.pdf',EVO_IMMUNO_POP),width=3.6,height=4.5)
library(RColorBrewer)
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

###########################################################
##      eQTL R2 by gene class (detected eQTLS only)      ##
###########################################################
DT=fread(sprintf("%s/Maxime/miRNA_V2/data/15.eQTL_comparisons/SourceData_eQTL_R2_byGeneClass.tsv",EVO_IMMUNO_POP),sep='\t')
pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/eQTL_R2_byGeneClass.pdf',EVO_IMMUNO_POP),width=3.6,height=4.5)
library(RColorBrewer)
colorSet=c(brewer.pal(8,'Set2'),brewer.pal(8,'Pastel2'),brewer.pal(8,'Dark2'),brewer.pal(8,'Set1'),brewer.pal(8,'Pastel1'))
#colorGeneType=c('grey',colorSet[c(1,3,38:37,33,6,25)]) # color scheme1
colorGeneType=c('grey',colorSet[c(1,3,6,2,25,36,28)]) # color scheme2

names(colorGeneType)=c('All','lncRNA','protein-coding','LOF Null','LOF Recessive','LOF Intolerant','TF','miRNA')
DT[,class:=factor(class,levels=c('All','lncRNA','protein-coding','LOF Null','LOF Recessive','LOF Intolerant','TF','miRNA'))]
p <- ggplot(DT[class!='All'],aes(y=eQTL_R2,x=class,fill=class))+geom_violin(scale='width')
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
p <- p + scale_fill_manual(values=colorGeneType)
p <- p + geom_boxplot(notch=TRUE,fill='#FFFFFF88')
p <- p + xlab('') + ylab('R2 of best eQTL')
print(p)
dev.off()

#########################################################################
##      any link between miRNA eQTL presence/R2 and miRNA expression ? ##
#########################################################################
print("loading miRNA data")
corrected_count = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""), sep="\t")

mir_count_melted = melt(corrected_count, id = c("ID"),
                               variable.name = "sample",
                               value.name = "count")
mean_count=mir_count_melted[,.(count=mean(count)),by=ID]
mean_count$R2=0
mean_count$R2[match(miRqtl_cis_annot$miRNA,mean_count$ID)]=miRqtl_cis_annot$R2
#cor.test(mean_count$R2,mean_count$count,method='s')  p = 0.28 
#cor.test(as.numeric(mean_count$R2>0),mean_count$count,method='s')  p = 0.32

#### Conclusion: no significant relationship between eQTL detection, or R2, and miRNA expression
# hist(width(reduce(PIRs_GR))

############################################
##      Visualization_Prom_Enhancer_PIR   ##
############################################

pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/Visualization_Prom_Enhancer_PIR_TLR1_6_10_locus.pdf',EVO_IMMUNO_POP),width=5,height=4)
oldmar=par()$mar
par(mar=c(4,8,1,1))
plot(c(38,38.2),c(0,4),ylab='',xlab='Genomic location',axes=F,col='#00000000');axis(1)
axis(2,at=1:3,labels=c('PIRs','enhancers','promoters'),las=2)
x=mapply(function(st,end){rect(st/1e6,.7,end/1e6,1.3,border='#00000044',col='#FF000044')},PIRs[oeChr==4 & oeEnd>38e6 & oeStart<40e6,oeStart],PIRs[oeChr==4  & oeEnd>38e6 & oeStart<40e6,oeEnd])
x=mapply(function(st,end){rect(st/1e6,1.7,end/1e6,2.3,border='#00000044',col='#00FF0044')},enhancers_monocytes[seqnames=='chr4' & end>38e6 & start<40e6,start],enhancers_monocytes[seqnames=='chr4'  & end>38e6 & start<40e6,end])
x=mapply(function(st,end){rect(st/1e6,2.7,end/1e6,3.3,border='#00000044',col='#0000FF44')},promoters_monocytes[seqnames=='chr4' & end>38e6 & start<40e6,start],promoters_monocytes[seqnames=='chr4'  & end>38e6 & start<40e6,end])
par(mar=oldmar)
dev.off()

############################################
##      Enrichments Prom_Enhancer_PIR     ##
############################################
colorSet=c(brewer.pal(8,'Set2'),brewer.pal(8,'Pastel2'),brewer.pal(8,'Dark2'),brewer.pal(8,'Set1'),brewer.pal(8,'Pastel1'))
colorSet_ordered=colorSet[c(rep(c(1,0,2)*8,8)+rep(1:8,e=3),24+rep((1:0)*8,8)+rep(1:8,e=2))]]
colReg = c("#CCCCCC", "#B3E2CD", "#66C2A5", "#E6F5C9", "#A6D854", "#8DA0CB", "#F4CAE4", "#E78AC3", "#FDCDAC", "#FC8D62")
names(colReg)= c("other","Enhancer (long)", "Enhancer (long) - PIR", "Enhancer (short)", "Enhancer (short) - PIR",
                 "PIR", "Promoter (long)", "Promoter (long) - PIR", "Promoter (short)", "Promoter (short) - PIR"))


####### make pie charts 

pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/eQTL_miRQTL_Prom_Enhancer_PIR_observed.pdf',EVO_IMMUNO_POP),width=5,height=5)
oldmar=par()$mar
layout(matrix(1:4,2))
par(mar=c(0.1,2,2,.1))
tab0=table(frequent_snps[,RegElt])
pie(tab0,col=colReg,main='all SNPs, MAF>5%',labels=NA)
names(tab0)[1]='other'
tab1=table(frequent_snps[is_miRQTL==TRUE,RegElt])
pie(tab1,col=colReg,main='miRQTL',labels=NA)
names(tab1)[1]='other'
tab2=table(frequent_snps[is_EQTL==TRUE,RegElt])
pie(tab2,col=colReg,main='eQTL',labels=NA)
names(tab2)[1]='other'
plot.new()
legend('top',bty='n',fill=colReg[-1],legend=names(colReg)[-1], ncol=1)
par(mar=oldmar)
dev.off()

####### compute Odds ratios 

OR_miRQTL=do.call(rbind,lapply(names(colReg),function(x){
                        if(!x %in% names(tab1)) 
                             oframe <- data.frame(Category=x,LowerCI = NA, OR = NA, UpperCI = NA, P=NA)
                        else{
                            myTable=matrix(c(tab1[x],sum(tab1)-tab1[x],tab0[x],sum(tab0)-tab0[x]),2,byrow=T)
                            test=fisher.test(myTable)
                            oframe <- data.frame(Category=x,LowerCI = test$conf.int[1], OR = test$est, UpperCI = test$conf.int[2], P=test$p.value)
                        }
                        oframe}))

OR_eQTL=do.call(rbind,lapply(names(colReg),function(x){
                        myTable=matrix(c(tab2[x],sum(tab2)-tab2[x],tab0[x],sum(tab0)-tab0[x]),2,byrow=T)
                        test=fisher.test(myTable)
                        oframe <- data.frame(Category=x,LowerCI = test$conf.int[1], OR = test$est, UpperCI = test$conf.int[2], P=test$p.value)
                        oframe}))

fwrite(OR_eQTL,file=sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/eQTL_RegElt_Enrichment.txt',EVO_IMMUNO_POP),sep='\t')
fwrite(OR_miRQTL,file=sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/miRQTL_RegElt_Enrichment.txt',EVO_IMMUNO_POP),sep='\t')

####### make Odds ratio forest plots 
pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/eQTL_RegElt_Enrichment.pdf',EVO_IMMUNO_POP),height=5,width=1.8)
p <- ggplot(OR_eQTL,aes(Category,log2(OR),col=Category)) + geom_pointrange( ymin = log2(OR_eQTL$LowerCI),ymax= log2(OR_eQTL$UpperCI)) + ylim(c(-1,6)) + scale_color_manual(values=colReg)
p <- p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),axis.text.y = element_text(angle = 90, hjust = .5, vjust=0.5)) 
p <- p + geom_hline(aes(yintercept=0)) + xlab('')+ theme(legend.position = "none")
p <- p + theme(panel.grid.minor = element_blank())
print(p)
dev.off()
pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/miRQTL_RegElt_Enrichment.pdf',EVO_IMMUNO_POP),height=5,width=1.8)
p <- ggplot(OR_miRQTL,aes(Category,log2(OR),col=Category)) + geom_pointrange( ymin = log2(OR_miRQTL$LowerCI),ymax= log2(OR_miRQTL$UpperCI)) + ylim(c(-1,6)) + scale_color_manual(values=colReg)
p <- p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),axis.text.y = element_text(angle = 90, hjust = .5, vjust=0.5)) 
p <- p + geom_hline(aes(yintercept=0)) + xlab('')+ theme(legend.position = "none")
p <- p + theme(panel.grid.minor = element_blank())
print(p)
dev.off()

####### compute Odds ratios (grouping all categories)
colReg_keyWords = colReg[c(3,6,8)]
names(colReg_keyWords)=c('Promoter','Enhancer','PIR')
OR_eQTL_keywords=do.call(rbind,lapply(c('Promoter','Enhancer','PIR'),function(x){
                        myTable=matrix(c(sum(tab2[grep(x,names(tab2))]),sum(tab2[-grep(x,names(tab2))]),sum(tab0[grep(x,names(tab0))]),sum(tab0[-grep(x,names(tab0))])),2,byrow=T)
                        test=fisher.test(myTable)
                        oframe <- data.frame(Category=x, LowerCI = test$conf.int[1], OR = test$est, UpperCI = test$conf.int[2], P=test$p.value)
                        oframe}))

OR_miRQTL_keywords=do.call(rbind,lapply(c('Promoter','Enhancer','PIR'),function(x){
                        myTable=matrix(c(sum(tab1[grep(x,names(tab1))]),sum(tab1[-grep(x,names(tab1))]),sum(tab0[grep(x,names(tab0))]),sum(tab0[-grep(x,names(tab0))])),2,byrow=T)
                        test=fisher.test(myTable)
                        oframe <- data.frame(Category=x, LowerCI = test$conf.int[1], OR = test$est, UpperCI = test$conf.int[2], P=test$p.value)
                        oframe}))

fwrite(OR_eQTL_keywords,file=sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/eQTL_RegElt_Enrichment_keywords.txt',EVO_IMMUNO_POP),sep='\t')
fwrite(OR_miRQTL_keywords,file=sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/miRQTL_RegElt_Enrichment_keywords.txt',EVO_IMMUNO_POP),sep='\t')

####### make Odds ratio forest plots  (grouping all categories)
pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/miRQTL_RegElt_Enrichment_keywords.pdf',EVO_IMMUNO_POP),height=5,width=1.8)
p <- ggplot(OR_miRQTL_keywords,aes(Category,log2(OR),col=Category)) + geom_pointrange( ymin = log2(OR_miRQTL_keywords$LowerCI),ymax= log2(OR_miRQTL_keywords$UpperCI)) + ylim(c(-1,6)) + scale_color_manual(values=colReg_keyWords)
p <- p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),axis.text.y = element_text(angle = 90, hjust = .5, vjust=0.5)) 
p <- p + geom_hline(aes(yintercept=0)) + xlab('')+ theme(legend.position = "none")
p <- p + theme(panel.grid.minor = element_blank())
print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/eQTL_RegElt_Enrichment_keywords.pdf',EVO_IMMUNO_POP),height=5,width=1.8)
p <- ggplot(OR_eQTL_keywords,aes(Category,log2(OR),col=Category)) + geom_pointrange( ymin = log2(OR_eQTL_keywords$LowerCI),ymax= log2(OR_eQTL_keywords$UpperCI)) + ylim(c(-1,6)) + scale_color_manual(values=colReg_keyWords)
p <- p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),axis.text.y = element_text(angle = 90, hjust = .5, vjust=0.5)) 
p <- p + geom_hline(aes(yintercept=0)) + xlab('')+ theme(legend.position = "none")
p <- p + theme(panel.grid.minor = element_blank())
print(p)
dev.off()


###### No significant diference between  eQTL and miR_QTLs
tab=cbind(table(grepl('Enhancer',frequent_snps[is_miRQTL==TRUE,RegElt])),table(grepl('Enhancer',frequent_snps[is_EQTL==TRUE,RegElt])));odds.ratio(tab)
#             LowerCI       OR  UpperCI alpha         P
#odds ratio 0.7909155 1.415138 2.719653  0.05 0.2933827
tab=cbind(table(grepl('Prom',frequent_snps[is_miRQTL==TRUE,RegElt])),table(grepl('Prom',frequent_snps[is_EQTL==TRUE,RegElt])));odds.ratio(tab)
#             LowerCI       OR  UpperCI alpha         P
#odds ratio 0.6037184 1.122042 2.275623  0.05 0.8804294
tab=cbind(table(grepl('long',frequent_snps[is_miRQTL==TRUE,RegElt])),table(grepl('long',frequent_snps[is_EQTL==TRUE,RegElt])));odds.ratio(tab)
#             LowerCI       OR  UpperCI alpha         P
#odds ratio 0.7692566 1.253059 2.112331  0.05 0.4200558
tab=cbind(table(grepl('short',frequent_snps[is_miRQTL==TRUE,RegElt])),table(grepl('short',frequent_snps[is_EQTL==TRUE,RegElt])));odds.ratio(tab)
#             LowerCI       OR UpperCI alpha         P
#odds ratio 0.5295133 1.622452 8.08921  0.05 0.6293457
tab=cbind(table(grepl('PIR',frequent_snps[is_miRQTL==TRUE,RegElt])),table(grepl('PIR',frequent_snps[is_EQTL==TRUE,RegElt])));odds.ratio(tab)
#             LowerCI       OR  UpperCI alpha         P
#odds ratio 0.8598217 1.476181 2.676828  0.05 0.1810338


############################################################
##      localization of eQTLs (TSS, TES, exon, Intron     ##
############################################################

eQTL_cis_annot[,RegElt_Roadmap:=frequent_snps[match(snps,snp_name),RegElt]]

eQTL_cis_annot[,typeDetail:=ifelse(CisDist_Gene==0, 
                                ifelse(abs(CisDist_5UTR)==0,'Exonic (5UTR)',
                                ifelse(abs(CisDist_3UTR)==0,'Exonic (3UTR)',
                                ifelse(abs(CisDist_CDS)==0,'Exonic (CDS)',
                                ifelse(abs(CisDist_Exon)==0,'Exonic (non-coding transcript)',
                                ifelse(abs(CisDist_CDS)<.1 & abs(CisDist_CDS)<abs(CisDist_3UTR) & abs(CisDist_CDS)<abs(CisDist_5UTR),'intronic (CDS-flanking)',
                                ifelse(abs(CisDist_3UTR)<.1 & abs(CisDist_3UTR)<abs(CisDist_CDS) & abs(CisDist_3UTR)<abs(CisDist_TES),'intronic (3UTR-flanking)',
                                ifelse(abs(CisDist_5UTR)<.1 & abs(CisDist_5UTR)<abs(CisDist_CDS) & abs(CisDist_5UTR)<abs(CisDist_TSS),'intronic (5UTR-flanking)',
                                ifelse(abs(CisDist_Exon)<.1 & abs(CisDist_Exon)<abs(CisDist_TSS) & abs(CisDist_Exon)<abs(CisDist_TES),'intronic (non coding-flanking)',
                                ifelse(abs(CisDist_TSS)<1 & abs(CisDist_TSS)<abs(CisDist_5UTR),'intronic (secondary TSS)',
                                ifelse(abs(CisDist_TES)<1 & abs(CisDist_TSS)<abs(CisDist_3UTR),'intronic (secondary TES)',
                                    'intronic')))))))))),
                            ifelse(abs(CisDist_TSS_gene) < 1 & CisDist_Gene<0 ,'TSS-flanking (proximal)',
                            ifelse(abs(CisDist_TSS_gene) < 20 & CisDist_Gene<0 ,'TSS-flanking',
                            ifelse(abs(CisDist_TES_gene) < 1 & CisDist_Gene>0 ,'TES-flanking (proximal)',
                            ifelse(abs(CisDist_TES_gene) < 20 & CisDist_Gene>0 ,'TES-flanking',
                            'Distant')))))]

eQTL_cis_annot[,type:=ifelse(CisDist_Gene==0, 
                                ifelse(abs(CisDist_Exon)==0,'Exonic',
                                ifelse(abs(CisDist_Exon)<.1,'Intronic (Exon flanking)',
                                ifelse(abs(CisDist_TSS_gene) < 20 ,'Intronic (TSS-flanking)',
                                ifelse(abs(CisDist_TES_gene) < 20 ,'Intronic (TES-flanking)',
                                    'Intronic')))),
                            ifelse(abs(CisDist_TSS_gene) < 20 & CisDist_Gene<0 ,'TSS-flanking',
                            ifelse(abs(CisDist_TES_gene) < 20 & CisDist_Gene>0 ,'TES-flanking',
                            'Distant')))]



colType=c("#CCCCCC","#66C2A5", "#B3E2CD","#FFD92F","#E41A1C","#FBB4AE","#CBD5E8","#8DA0CB")
names(colType)=c('Distant','TSS-flanking','Intronic (TSS-flanking)','Intronic','Exonic','Intronic (Exon flanking)','Intronic (TES-flanking)','TES-flanking')


####### make pie charts 
pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/eQTL_TSS_Exon_Intron_localization_observed.pdf',EVO_IMMUNO_POP),width=5,height=5)
oldmar=par()$mar
layout(matrix(1:4,2))
# 5% FDR
par(mar=c(0.1,2,2,.1))
tab1=table(eQTL_cis_annot[FDR<0.05,type])
pie(tab1[names(colType)],col=colType,main='eQTL 5% FDR',labels=NA)

# use insignificant eQTLs as control (p > 1e-5)
tab2=table(eQTL_cis_annot[pvalue>1e-5,type])
pie(tab2[names(colType)],col=colType,main='eQTL p>1e-5',labels=NA)
plot.new()
# legend
legend('top',bty='n',fill=colType,legend=names(colType), ncol=1)
par(mar=oldmar)
dev.off()



####### compare localization (exon, TSS, TES, intron) with promoter definitions)
pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/eQTL_TSS_Exon_Intron_localization_VS_Prom_Enh_PIR_count_and_OR.pdf',EVO_IMMUNO_POP),width=5,height=5)
par(mar=c(10,14,1,1))
tab=table(eQTL_cis_annot[FDR<.05,RegElt_Roadmap], eQTL_cis_annot[FDR<.05,typeDetail])
Image(log2(tab+.1))
Image(log2(t(tab/apply(tab,1,sum))/apply(tab,2,sum)*sum(tab)+1e-1))
par(mar=oldmar)
dev.off()

####### compute enrichments
OR=do.call(rbind,lapply(names(colType),function(x){
                        myTable=matrix(c(tab1[x],sum(tab1)-tab1[x],tab2[x],sum(tab2)-tab2[x]),2,byrow=T)
                        test=fisher.test(myTable)
                        oframe <- data.frame(Category=x,LowerCI = test$conf.int[1], OR = test$est, UpperCI = test$conf.int[2], P=test$p.value)
                        oframe}))
}
OR$Category=factor(OR$Category,OR$Category)

####### Forest plots Odds ratios
pdf(sprintf('%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/eQTL_TSS_Exon_Intron_Enrichment.pdf',EVO_IMMUNO_POP),height=5,width=1.8)
p <- ggplot(OR,aes(Category,log2(OR),col=Category)) + geom_pointrange( ymin = log2(OR$LowerCI),ymax= log2(OR$UpperCI)) + ylim(c(-1,6)) + scale_color_manual(values=colType)
p <- p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),axis.text.y = element_text(angle = 90, hjust = .5, vjust=0.5)) 
p <- p + geom_hline(aes(yintercept=0)) + xlab('')+ theme(legend.position = "none")
p <- p + theme(panel.grid.minor = element_blank())
print(p)
dev.off()


