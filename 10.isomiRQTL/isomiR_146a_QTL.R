###############################################################################
##       infos on isomiR SNP rs290164 that controls isomiR of mir-146a 3p    ##
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

isomiR_QTL=fread(sprintf('%s/Maxime/miRNA_V2/data/10.isomiRQTL/cis_isomiRQTLs_with_FDR_filtered_best_isomiRNA_snp_association.tsv',EVO_IMMUNO_POP))



#################################
##  load isomiR ratio Data     ##
#################################

isomir_ratio_cast=fread(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_ratios_aggregated_nosubs.GCRL_Batch_lane_corrected.tsv",EVO_IMMUNO_POP))
#isomir_ratio_cast=isomir_ratio_cast[-grep('other',isomir_ID)]
#dim(isomir_ratio_cast)
isomir_ratio_melted = melt(isomir_ratio_cast, id = c("mirID", "isomir_ID"),
                               variable.name = "sample",
                               value.name = "ratio")
isomir_ratio_melted[, individual := substr(sample, 1,6)]
isomir_ratio_melted[, condition := condIndex[as.numeric(substr(sample, 8,8))]]
isomir_ratio_melted[, population := substr(sample, 1,3)]
isomir_ratio_melted[,hsa_ID := gsub('(hsa.*)_(MIMAT.*)_(MI[0-9]+)','\\1',mirID)]
inverseNormalRankTransform=function(x){n=length(x);
									qnorm(frank(x,na.last=FALSE,ties.method='random')/(n+1),mean(x,na.rm=T),sd(x,na.rm=T))}
isomir_ratio_melted[, ratio_transformed := inverseNormalRankTransform(ratio),by=.(mirID,isomir_ID)]

# get meanRatio NS
meanRatio=isomir_ratio_melted[condition==1,.(baseRatio=mean(ratio,na.rm=T)),by=isomir_ID]

#################################
##  load isomiR count Data     ##
#################################

isomir_count_cast=fread(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_counts_aggregated_nosubs.log2RPM.GCRL_Batch_lane_corrected.tsv",EVO_IMMUNO_POP))
isomir_count_melted = melt(isomir_count_cast, id = c("V1"),
                               variable.name = "sample",
                               value.name = "count")
isomir_count_melted[, individual := substr(sample, 1,6)]
isomir_count_melted[, condition := condIndex[as.numeric(substr(sample, 8,8))]]
isomir_count_melted[, population := substr(sample, 1,3)]
isomir_count_melted[, hsa_ID := gsub('.*(hsa.*)_(MIMAT.*)_(MI[0-9]+).*','\\1',V1)]
isomir_count_melted[, mirID := gsub('.*(hsa.*)_(MIMAT.*)_(MI[0-9]+).*','\\1_\\2_\\3',V1)]
isomir_count_melted[, isomir_ID := V1]
inverseNormalRankTransform=function(x){n=length(x);
									qnorm(frank(x,na.last=FALSE,ties.method='random')/(n+1),mean(x,na.rm=T),sd(x,na.rm=T))}
isomir_count_melted[, count_transformed := inverseNormalRankTransform(count),by=.(mirID,isomir_ID)]


#########################################################
##### load count data + genotypes for a given QTL #######
#########################################################


make_data_frame_to_study <- function(data.table.line,alleles=F){
 snp_name = data.table.line$snps[1]
  isomiR = data.table.line$gene[1]
  data_temp=isomir_ratio_melted[isomir_ID==isomiR,.(mirID,isomir_ID,individual,condition,population,ratio_transformed,ratio)]
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



#########################################################################
############### perform mediation analysis on all mirQTLS ###############
#########################################################################

w=which(grepl('hsa-miR-146a-3p',isomiR_QTL$gene) & isomiR_QTL$snps=='rs2910164')[1:2]

library(mediation)
res=list()
for (line_n in w){
    DF=make_data_frame_to_study(isomiR_QTL[line_n])
    for (cond in condIndex){
        mod_miR_pop=lm(ratio_transformed~population+genotype,data=DF[condition==cond])
        mod_pop=lm(genotype~population,data=DF[condition==cond])
        mod_miR=lm(ratio_transformed~population,data=DF[condition==cond])
	    med=mediate(mod_pop, mod_miR_pop,treat="population", mediator="genotype")
    	res[[paste(line_n,cond)]]=data.table(mirQTL_nb=line_n,condition=cond,delta_pop=med$tau.coef,Pct_mediated=med$n1*100,Pct_mediated_inf=med$n1.ci[1]*100,Pct_mediated_sup=med$n1.ci[2]*100,betaSNP=mod_miR_pop$coef[3])
		print(res[[paste(line_n,cond)]])
    }
}

res=rbindlist(res)
mean(res$Pct_mediated)
Pct_mediation = dcast(res,mirQTL_nb~condition,value.var=c('delta_pop','Pct_mediated','Pct_mediated_inf','Pct_mediated_sup','betaSNP'))
fwrite(Pct_mediation,file=sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/Percentage_mediation.txt",EVO_IMMUNO_POP))



