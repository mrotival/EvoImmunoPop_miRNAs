
###############################################################################
##       infos on miR cluster SNP rs12881760 that is under selection         ##
###############################################################################


##############
##Libariries##
##############
require(ggplot2)
luq = function(x){length(unique(x))}
colERC5 = c("#525252AA","#E31A1CAA","#33A02CAA","#1F78B4AA","#6A3D9AAA")
colERC = c("#969696", "#525252", "#FB9A99", "#E31A1C", "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4", "#CAB2D6", "#6A3D9A")
condIndex = c("NS","LPS","PAM3CSK4","R848","IAV" )

names(colERC5) = condIndex
names(colERC)=paste(rep(condIndex,e=2),'-',rep(c('EUB','AFB'),5))


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



getGenos_locus=function(chr, start, end,path=EVO_IMMUNO_POP,phased=FALSE){
	require(dplyr)
	mydb=src_sqlite(paste(path,'/Maxime/SNP_annotation/imputed_Genotypes/GenoDataBase.db',sep=''))
#	if(is.null(chrlist)){
#	map=tbl(mydb,sql(paste('SELECT * FROM Map_imputed WHERE "chromosome"==',chr,' AND "position" >',start,' AND "position" <',end ,sep=''))) %>%collect()
##	map=filter(Map,chromosome==chr & position > start & position < end ) %>% collect()
#	chrlist=unique(map$chromosome)
#	}
#	RES=NULL
#	for(CHR in chrlist){
#		res=tbl(mydb,sql(paste('SELECT * FROM ',ifelse(phased,'Haplo','Geno'),'_',CHR,' WHERE "position" >',start,' AND "position" <',end ,sep=''))) %>%collect()
		RES=tbl(mydb,sql(paste('SELECT * FROM ',ifelse(phased,'Haplo','Geno'),'_',chr,' WHERE "position" >',start,' AND "position" <',end ,sep=''))) %>%collect(n=Inf)
#		RES=rbind(RES,res)
#		}
#	RES
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


#########################################################
################### load count data #####################
#########################################################
corrected_count = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""), sep="\t")

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



#########################################################
##### load count data + genotypes for a given QTL #######
#########################################################


make_data_frame_to_study <- function(data.table.line,alleles=F){
 snp_name = data.table.line$snps[1]
  miRNA = data.table.line$miRNA[1]
  data_temp=mir_count_melted[ID==miRNA,.(ID,individual,condition,population,count_transformed,count)]
  map=getMapInfo(snp_name)
  genotypes=getSNP(snp_name)
  if(alleles){
    alleles=c(map[,'allele.1'],map[,'allele.2'])
    if(nchar(alleles)[1]>1){
        alleles=c('ins','-')
    }
    if(nchar(alleles)[2]>1){
        alleles=c('del','-')
    }
    geno_allele=paste(c(alleles[1],alleles[2],alleles[2]),c(alleles[1],alleles[1],alleles[2]),sep='/')
    inds=names(genotypes)
    genotypes = geno_allele[1+ genotypes]
    data_temp[, genotype := factor(genotypes[match(data_temp$individual, inds)],rev(geno_allele))]
    }else{
    genotypes=2-genotypes
    data_temp[, genotype := genotypes[match(data_temp$individual, names(genotypes))]]
    }
#    data_temp=data_temp[which(!is.na(genotypes)),]
  return(data_temp)
}


#########################################################
############@### load miR QTL informations ##############
#########################################################


miRqtl_cis_annot=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4A_mirQTLs_Annoted.tsv",EVO_IMMUNO_POP))
cols=c("snps",'miRNA', "daf_char_EUB","daf_char_AFB","FST_adj",'iHS_EUB','mat_strand','mat_start','mat_end','CisDist_miRNA','TSS','CisDist_TSS','RegElt','TFBS','typeDetail','Promoter.name')


#########################################################################
############### perform mediation analysis on all mirQTLS ###############
#########################################################################


library(mediation)
res=list()
for (line_n in 1:nrow(miRqtl_cis_annot)){
    DF=make_data_frame_to_study(miRqtl_cis_annot[line_n])
    for (cond in condIndex){
        mod_miR_pop=lm(count_transformed~population+genotype,data=DF[condition==cond])
        mod_pop=lm(genotype~population,data=DF[condition==cond])
        mod_miR=lm(count_transformed~population,data=DF[condition==cond])
	    med=mediate(mod_pop, mod_miR_pop,treat="population", mediator="genotype")
    	res[[paste(line_n,cond)]]=data.table(mirQTL_nb=line_n,condition=cond,delta_pop=med$tau.coef,Pct_mediated=med$n1*100,Pct_mediated_inf=med$n1.ci[1]*100,Pct_mediated_sup=med$n1.ci[2]*100,betaSNP=mod_miR_pop$coef[3])
		print(res[[paste(line_n,cond)]])
    }
}

res=rbindlist(res)
Pct_mediation = dcast(res,mirQTL_nb~condition,value.var=c('delta_pop','Pct_mediated','Pct_mediated_inf','Pct_mediated_sup','betaSNP'))
fwrite(Pct_mediation,file=sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/Percentage_mediation.txt",EVO_IMMUNO_POP))

# for (i in 1:12){
#     DF=make_data_frame_to_study(miRqtl_cis_annot[which(snps=='rs12881760')][i])
#     mod_miR=lm(count_transformed~population+genotype,data=DF[condition=='NS'])
#     mod_pop=lm(genotype~population,data=DF[condition=='NS'])
#     med=mediate(mod_pop, mod_miR,treat="population", mediator="genotype")
#     cat(miRqtl_cis_annot[which(snps=='rs12881760'),get('miRNA')][i],med$n1*100,'\n')
# }

########################################################
############### extract chrom 14 cluster ###############
########################################################



miRqtl_cis_annot[snps=='rs12881760',]

chr14_miRNAs=miRqtl_cis_annot[snps=='rs12881760',miRNA]

#         snps           miRNA         minP   maxBeta daf_char_EUB daf_char_AFB FST_adj  iHS_EUB mat_strand mat_start   mat_end CisDist_miRNA       TSS CisDist_TSS            RegElt             TFBS typeDetail                     Promoter.name
# 1: rs12881760  hsa-miR-127-3p 4.257985e-10 0.7199187         0.73        0.035   0.675 -3.10506          + 101349372 101349393     -0.173037        NA          NA CTCF Binding Site CTCF,RAD21,BACH1       <NA>
# 2: rs12881760  hsa-miR-134-5p 7.685317e-07 0.4193880         0.73        0.035   0.675 -3.10506          + 101521031 101521052     -0.344696        NA          NA CTCF Binding Site CTCF,RAD21,BACH1       <NA>
# 3: rs12881760  hsa-miR-136-3p 1.265542e-08 0.5329119         0.73        0.035   0.675 -3.10506          + 101351087 101351108     -0.174752        NA          NA CTCF Binding Site CTCF,RAD21,BACH1       <NA>
# 4: rs12881760  hsa-miR-370-3p 1.660737e-08 0.4730486         0.73        0.035   0.675 -3.10506          + 101377523 101377544     -0.201188 101364123   -0.187788 CTCF Binding Site CTCF,RAD21,BACH1    Distant	p@chr14:101364120..101364126,+
# 5: rs12881760  hsa-miR-381-3p 1.114029e-08 0.5588978         0.73        0.035   0.675 -3.10506          + 101512305 101512326     -0.335970        NA          NA CTCF Binding Site CTCF,RAD21,BACH1       <NA>
# 6: rs12881760  hsa-miR-409-3p 1.708642e-10 0.6739905         0.73        0.035   0.675 -3.10506          + 101531683 101531704     -0.355348        NA          NA CTCF Binding Site CTCF,RAD21,BACH1       <NA>
# 7: rs12881760  hsa-miR-410-3p 1.531626e-10 0.6856698         0.73        0.035   0.675 -3.10506          + 101532298 101532318     -0.355963        NA          NA CTCF Binding Site CTCF,RAD21,BACH1       <NA>
# 8: rs12881760  hsa-miR-411-5p 6.187000e-09 0.5915436         0.73        0.035   0.675 -3.10506          + 101489677 101489697     -0.313342        NA          NA CTCF Binding Site CTCF,RAD21,BACH1       <NA>
# 9: rs12881760  hsa-miR-431-5p 6.297919e-09 0.5099542         0.73        0.035   0.675 -3.10506          + 101347363 101347383     -0.171028        NA          NA CTCF Binding Site CTCF,RAD21,BACH1       <NA>
#10: rs12881760  hsa-miR-432-5p 8.932399e-10 0.5247004         0.73        0.035   0.675 -3.10506          + 101350833 101350855     -0.174498        NA          NA CTCF Binding Site CTCF,RAD21,BACH1       <NA>
#11: rs12881760 hsa-miR-487b-3p 2.306316e-07 0.4713975         0.73        0.035   0.675 -3.10506          + 101512842 101512863     -0.336507        NA          NA CTCF Binding Site CTCF,RAD21,BACH1       <NA>
#12: rs12881760  hsa-miR-654-3p 2.267917e-09 0.5890732         0.73        0.035   0.675 -3.10506          + 101506606 101506627     -0.330271        NA          NA CTCF Binding Site CTCF,RAD21,BACH1       <NA>



#### how many miRNA in the chr14 miRNA cluster
miRNA_coordinate = fread(paste(EVO_IMMUNO_POP, "ERCPilot_SharedDBs/mirbase20/miRNA_mature_coordinates_strandinfo.bed", sep=""))
names(miRNA_coordinate) = c("chromosome", "start", "end", "miRNA_name", "V5", "V6", "V7", "V8")
miRNA_coordinate[chromosome=='chr14' & start>101e6 & end<102e6]
# 97

######################################################################
############### plot selection in the chrom 14 cluster ###############
######################################################################

pos=getMapInfo('rs12881760')$position
Map_locus_chr14=getMapInfo_locus(14,100.9e6,101.7e6)
RR=getRecombinationRate_locus(14,100.9e6,101.7e6)

pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/chr14_rs12881760_locus_selection.pdf",EVO_IMMUNO_POP),height=5,width=4.5)
layout(1:2)
par(mar=c(4,4,.2,4))
w=which(Map_locus_chr14$maf_EUB>.05)
plot(Map_locus_chr14$position/1e6, Map_locus_chr14$FST_adj,cex=0.5+Map_locus_chr14$FST_adj,pch=16,col=colERC5[1],ylim=c(0,1),ylab='FST',xlab='',las=1)
library(zoo)
w=which(Map_locus_chr14$SNPfreq>.1)
lines(RR$Position.bp./1e6,RR$Rate.cM.Mb./max(RR$Rate.cM.Mb.),col='lightblue')
lines(rollmean(Map_locus_chr14$position[w]/1e6,100),rollmean(Map_locus_chr14$FST_adj[w]>.4,100),col='red',lwd=3)
abline(v=pos/1e6,col='blue',lty=3)
abline(h=.4,lty=3,col='lightgrey')
axis(4,at=seq(0,1,l=5),labels=seq(0,100,l=5),las=2)

Map_locus_chr14$iHS_EUB_rank=frank(-abs(Map_locus_chr14$iHS_EUB),na.last=FALSE, ties.method='random')/nrow(Map_locus_chr14)
w=which(Map_locus_chr14$maf_EUB>.05 & !is.na(Map_locus_chr14$iHS_EUB))
plot(Map_locus_chr14$position[w]/1e6, abs(Map_locus_chr14$iHS_EUB[w]),cex=1.5-Map_locus_chr14$iHS_EUB_rank[w],pch=16,col=ifelse(Map_locus_chr14$iHS_EUB[w]>0,colERC5[2],colERC5[4]),ylim=c(0,6),ylab='|iHS|',xlab='Genomic position (Mb)',las=1)
library(zoo)
lines(RR$Position.bp./1e6,RR$Rate.cM.Mb./max(RR$Rate.cM.Mb.)*6,col='lightblue')
lines(rollmean(Map_locus_chr14$position[w]/1e6,100),rollmean(abs(Map_locus_chr14$iHS_EUB)[w]>2.5,100)*6,col='black',lwd=3)
abline(v=pos/1e6,col='blue',lty=3)
abline(h=2.5,lty=3,col='lightgrey')
axis(4,at=seq(0,6,l=5),labels=seq(0,100,l=5),las=2)
dev.off()

########################################################################################
############### Compare selection statistics to GenomeWide distributions ###############
########################################################################################

library(zoo)
Map_imputed=fread('/Volumes/evo_immuno_pop/Maxime/SNP_annotation/Map_imputed_essential_informations.txt')
Map_imputed[,daf_EU:=as.numeric(gsub('\\((.*)\\)','\\1',daf_char_EUB))]
Map_imputed[,daf_AF:=as.numeric(gsub('\\((.*)\\)','\\1',daf_char_AFB))]
Map_imputed[,maf_EUB:=pmin(daf_EU,1-daf_EU)]
Map_imputed[,maf_AFB:=pmin(daf_AF,1-daf_AF)]
Map_imputed[,sum(maf_EUB>=0.05 | maf_AFB>=0.05)]
#10261270

w=which(Map_imputed$maf_EUB>.05 & !is.na(Map_imputed$iHS_EUB))
window_length=(rollmax(Map_imputed$position[w],200)+rollmax(-Map_imputed$position[w],200))/1e3
mean(window_length[window_length<10000])


### compute % of outliers in a 100 SNP (~50 kb) window around the SNP 
Pct_outlier=rollmean(abs(Map_imputed$iHS_EUB[w])>2.5,101,na.rm=T)
quantile(Pct_outlier,.99) #
Map_imputed$Pct_outlier=NA
Map_imputed$Pct_outlier[w[51:(length(w)-50)]]=Pct_outlier
Map_imputed[snp.name=='rs12881760',]

####### plot the mirQTLs

data_temp=mir_count_melted[ID%in%chr14_miRNAs,.(ID,individual,condition,population,count_transformed,count)]
data_temp[,SNP:=2-getSNP('rs12881760')[individual]]
data_temp[,cond_pop:=paste(condition,population,sep=' - ')]
#data_temp=dcast(data_temp,SNP+individual+condition+population~ID,value.var='count_transformed')
pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/chr_14_clusterQTL.pdf",EVO_IMMUNO_POP),height=2.5,width=12)
data_temp=data_temp[,.(population=unique(population), condition=unique(condition) , SNP=unique(SNP) , cond_pop=unique(cond_pop) ,count=mean(count_transformed)),by= individual]
p <- ggplot(data_temp[!is.na(SNP)],aes(x=factor(SNP),y=count,fill=cond_pop))+geom_violin(scale='width')+facet_grid(~factor(cond_pop,names(colERC))) + scale_fill_manual(values=colERC)
    p <- p + geom_boxplot(fill="#FFFFFF88", outlier.size=0.5, notch=TRUE,width=0.4)
    p <- p + theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank() , legend.position = "none")
    print(p)
dev.off()

###############################################################################
############### look for a secondary mirQTL near the chr14 cluster ############
###############################################################################

Genos_chr14=getGenos_locus(14,101176335-1e6,101176335+1e6)
data_temp=mir_count_melted[ID%in%chr14_miRNAs,.(ID,individual,condition,population,count_transformed,count)]
data_temp=dcast(data_temp,individual+condition+population~ID,value.var='count_transformed')
mat_temp=as.matrix(data_temp[,-(1:3)])
annot=data_temp[,mget(c('individual','condition','population'))]
library(MatrixEQTL)

for(cond in condIndex){}
cond='NS'
RPM=SlicedData$new()
RPM$CreateFromMatrix(t(mat_temp[annot$cond==cond,]))    
colnames(RPM)=annot$individu[annot$cond==cond]

Geno=SlicedData$new()
Geno$CreateFromMatrix(as.matrix(Genos_chr14[,annot$individu[annot$cond==cond]]))    
rownames(Geno)=Genos_chr14$snp.name

Covar=SlicedData$new()
mySNP=2-as.numeric(unlist(Genos_chr14[Genos_chr14$snp.name=='rs12881760',annot$individu[annot$cond==cond]]))
myPop=as.numeric(annot[annot$cond==cond,population]=='AFB')
Covar$CreateFromMatrix(t(as.matrix(cbind(myPop,mySNP))))

    
# no other eQTL around

#############################################################
############### miR-targets of the chr14 cluster ############
#############################################################
# load GO annotation data
source(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/scripts/00a.GOannotation/01.Load_and_format_databases.R", sep=""))
GO_informations = load_BHF_UCL_miRNA()

# load targets data
targets=fread(sprintf('%s/Martin/miRNA/11.miRNA_binding_sites/results/miRNA_binding_sites_on_protein_coding_3primeUTR_simplified.tsv',EVO_IMMUNO_POP))
targets$per1=as.numeric(gsub('%','',targets$per1))
targets$per2=as.numeric(gsub('%','',targets$per1))
chr14_cluster_targets=unique(targets[which(miRNA%in%chr14_miRNAs & per1>80 & per2>80),],by=c('miRNA','EnsemblGeneID','GeneNames'))
GOSeq(names(table(chr14_cluster_targets$EnsemblGeneID))[table(chr14_cluster_targets$EnsemblGeneID)>3],unique(targets$EnsemblGeneID))

chr14_miRNAs_extended=miRNA_coordinate[chromosome=='chr14' & start>101e6 & end<102e6,miRNA_name]
chr14_cluster_targets_extended=unique(targets[which(miRNA%in%chr14_miRNAs_extended ),],by=c('miRNA','EnsemblGeneID','GeneNames'))

targets_cor=fread(sprintf('%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/glmnet_correlation_miRNA_geneRPKM_transcriptionAdjusted_svaAdjusted/gene_expression_glmnet_alpha0.5_allConds.tsv',EVO_IMMUNO_POP))
allGenes=unique(targets_cor$gene)
targets_cor=targets_cor[grep('hsa',selected_variables),]
# the SNP on chr 14 is associated with Semeno-gellin - 2
for (cond in 1:5){
    tab=table(targets_cor[condition==cond,selected_variables])
    for (miR in names(rev(sort(tab)))[1:7]){
        cat(cond,miR,'\n')
        target_miR=merge(targets_cor[condition==cond & selected_variables==miR,],targets[miRNA==miR,],by.x='gene',by.y='EnsemblGeneID',all.x=T,all.y=T)
        print(odds.ratio(table((allGenes %in% target_miR[!is.na(selected_variables),gene]),(allGenes %in% target_miR[!is.na(miRNA),gene]))))
    }
}




################################################################################################################
############### Compare selection statistics to GenomeWide FST European-African & European-Asian ###############
################################################################################################################

Map_1kg=list()
SELINK='/pasteur/projets/policy01/selink'
for(pop in c('CEU','YRI','CHB')){
	Map_1kg[[pop]]=NULL
	for(chr in 22:1){
		cat(chr,'')
		Selink2=fread(paste(SELINK,'/1000G_Ph1_ch/Selink_inputs/100k_allPOP/JOB_',chr,'/normalized_OUT_Selink2_chr',chr,'_',pop,'.out',sep=''))
		colnames(Selink2)=c('SNP','start','stop','core_id','core_pos','minor_allele','ancestral','Theta_Pi','Theta_S','Tajima_D','Fay&WuH','maf','daf','iHHa','iHHd','ihs','Deltat_iHH','piA','piD','piA/piD','iHS_norm','piA/piD_norm','Deltat_iHH_norm')
		Selink2=Selink2[-1,]
		Selink2=Selink2[Selink2$ancestral!='.',]
		Selink2=Selink2[Selink2$maf>0.05,]
		Selink2$daf=as.numeric(Selink2$daf)
		Selink2=Selink2[!is.na(Selink2$daf),]
		Selink2$posID=paste(chr,Selink2$core_pos,sep=':')
		cat(nrow(Selink2),'\n')
		Map_1kg[[pop]][[chr]]=Selink2
	}
	Map_1kg[[pop]]=rbindlist(Map_1kg[[pop]])
}
Map_1kg[['interpop']]=NULL
for(chr in 22:1){
    Selink2_interpop=fread(sprintf('%s/1000G_Ph1_ch/Selink_inputs/100k_allPOP/JOB_%s/OUT_Selink2_chr%s.interpop',SELINK,chr,chr))
	Map_1kg[['interpop']][[chr]]=Selink2_interpop
	}
Map_1kg[['interpop']]=rbindlist(Map_1kg[['interpop']])


plot(Selink2[core_pos>100.8e6 & core_pos<101.8e6,core_pos],Selink2[core_pos>100.8e6 & core_pos<101.8e6,FST_CEU_CHB])
points(Selink2[core_pos>100.8e6 & core_pos<101.8e6,core_pos],Selink2[core_pos>100.8e6 & core_pos<101.8e6,FST_CEU_YRI],col='blue')
 points(Selink2[core_pos>100.8e6 & core_pos<101.8e6,core_pos],Selink2[core_pos>100.8e6 & core_pos<101.8e6,FST_YRI_CHB],col='red')

cor_pos_rs12881760='14_101176335'

# obtain FSTs
mean(cor_pos_rs12881760)
FST_AFR=Map_1kg[['interpop']][core_id=='rs12881760',FST_CEU_YRI] 
FST_ASI=Map_1kg[['interpop']][core_id=='rs12881760',FST_CEU_CHB] 
DAF_CEU=Map_1kg[['CEU']][core_id=='rs12881760',daf] 

Map_merged = merge(Map_1kg[['interpop']],Map_1kg[['CEU']],by='core_id',all.x=T)

mean(Map_merged[(daf-DAF_CEU)<0.01,FST_CEU_YRI]>=FST_AFR,na.rm=T)
mean(Map_merged[(daf-DAF_CEU)<0.01,FST_CEU_CHB]>=FST_ASI,na.rm=T)
mean(Map_1kg[['interpop']][,FST_CEU_CHB]>=FST_ASI)



Map_imputed=fread('/Volumes/evo_immuno_pop/Maxime/SNP_annotation/Map_imputed_essential_informations.txt')
Map_imputed[,daf_EU:=as.numeric(gsub('\\((.*)\\)','\\1',daf_char_EUB))]
Map_imputed[,daf_AF:=as.numeric(gsub('\\((.*)\\)','\\1',daf_char_AFB))]
Map_imputed[,maf_EUB:=pmin(daf_EU,1-daf_EU)]
Map_imputed[,maf_AFB:=pmin(daf_AF,1-daf_AF)]
Map_imputed[,sum(maf_EUB>=0.05 | maf_AFB>=0.05)]
#10261270

w=which(Map_imputed$maf_EUB>.05 & !is.na(Map_imputed$iHS_EUB))
window_length=(rollmax(Map_imputed$position[w],200)+rollmax(-Map_imputed$position[w],200))/1e3
mean(window_length[window_length<10000])


### compute % of outliers in a 100 SNP (~50 kb) window around the SNP 
Pct_outlier=rollmean(abs(Map_imputed$iHS_EUB[w])>2.5,101,na.rm=T)
quantile(Pct_outlier,.99) #
Map_imputed$Pct_outlier=NA
Map_imputed$Pct_outlier[w[51:(length(w)-50)]]=Pct_outlier
Map_imputed[snp.name=='rs12881760',]

####### plot the mirQTLs

data_temp=mir_count_melted[ID%in%chr14_miRNAs,.(ID,individual,condition,population,count_transformed,count)]
data_temp[,SNP:=2-getSNP('rs12881760')[individual]]
data_temp[,cond_pop:=paste(condition,population,sep=' - ')]
#data_temp=dcast(data_temp,SNP+individual+condition+population~ID,value.var='count_transformed')
pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/chr_14_clusterQTL.pdf",EVO_IMMUNO_POP),height=2.5,width=12)
data_temp=data_temp[,.(population=unique(population), condition=unique(condition) , SNP=unique(SNP) , cond_pop=unique(cond_pop) ,count=mean(count_transformed)),by= individual]
p <- ggplot(data_temp[!is.na(SNP)],aes(x=factor(SNP),y=count,fill=cond_pop))+geom_violin(scale='width')+facet_grid(~factor(cond_pop,names(colERC))) + scale_fill_manual(values=colERC)
    p <- p + geom_boxplot(fill="#FFFFFF88", outlier.size=0.5, notch=TRUE,width=0.4)
    p <- p + theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank() , legend.position = "none")
    print(p)
dev.off()

