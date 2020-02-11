###############
##Librairires##
###############

sample_names=fread(sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/sample_names_977_highQuality.tsv",EVO_IMMUNO_POP))[[1]]
df_cond=table(substr(sample_names,8,8))

require(ggplot2)
luq = function(x){length(unique(x))}
colERC5 = c("#525252AA","#E31A1CAA","#33A02CAA","#1F78B4AA","#6A3D9AAA")
colERC = c("#969696", "#525252", "#FB9A99", "#E31A1C", "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4", "#CAB2D6", "#6A3D9A")
condIndex = c("NS","LPS","PAM3CSK4","R848","IAV" )

names(colERC5) = condIndex
names(colERC)=paste(rep(condIndex,e=2),'-',rep(c('EUB','AFB'),5))

#############
##Load data##
#############

isomiRqtl_cis = fread(sprintf("%s/Maxime/miRNA_V2/data/10.isomirQTL/cis_isomiRQTLs_with_FDR_filtered_best_isomiRNA_snp_association.tsv",EVO_IMMUNO_POP))
isomiRqtl_cis[,isomiR:=gene]
isomiRqtl_cis[,minP:=min(pvalue),by=isomiR]
isomiRqtl_cis[,FDR:=min(FDR_Perm),by=isomiR]
isomiRqtl_cis[,maxBeta:=max(abs(beta)),by=isomiR]
isomiRqtl_cis[,df_cond:=df_cond[condition]-3]
isomiRqtl_cis[,R2:=statistic^2/(statistic^2+df_cond)]
isomiRqtl_cis[,condition:=factor(condIndex[condition],levels=condIndex)]
isomiRqtl_cis[, miRNA := gsub('.*(hsa.*)_(MIMAT.*)_(MI[0-9]+).*','\\1',isomiR)]
isomiRqtl_cis[, miR_ID := gsub('.*(hsa.*)_(MIMAT.*)_(MI[0-9]+).*','\\1_\\2_\\3',isomiR)]

cols=c("snps", "isomiR","miRNA", "FDR", "beta_NS", "pvalue_NS", "R2_NS",
                          "beta_LPS", "pvalue_LPS", "R2_LPS", 
                          "beta_PAM3CSK4", "pvalue_PAM3CSK4", "R2_PAM3CSK4",
                          "beta_R848", "pvalue_R848", "R2_R848", 
                          "beta_IAV", "pvalue_IAV", "R2_IAV")

isomiRqtl_cis_table=dcast(isomiRqtl_cis,snps + isomiR + miRNA + minP + maxBeta + FDR ~ condition, value.var=c('beta', 'pvalue', 'R2'))
fwrite(isomiRqtl_cis_table[order(minP),mget(cols)],file=sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/SupTable4C_isomirQTLs_raw.tsv",EVO_IMMUNO_POP),sep='\t')
# isomiRqtl_cis_table = isomiRqtl_cis_table[order(minP),mget(cols)]

length(isomiRqtl_cis_table[,unique(miRNA)]) # 13
length(isomiRqtl_cis_table[,unique(isomiR)]) # 25

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

snps_information=getMapInfo(unique(isomiRqtl_cis_table$snps))[c("snp.name","chromosome","position","allele.2","allele.1","SNPfreq_AF","SNPfreq_EU","daf_char_AFB","daf_char_EUB","RegElt","TFBS", 'FST_adj','iHS_AFB','iHS_EUB','P_CLS2_EUB','P_CLS2_AFB','aSNP')]
fwrite(snps_information,file=sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/SupTable4C_isomirQTLs_snp_Annot.tsv",EVO_IMMUNO_POP),sep='\t')



mir_TSS_annot=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_TSS_annot_DeRie2017.txt',EVO_IMMUNO_POP))
##--------------------------------------------##
## private version (based on mysql database)  ##
##--------------------------------------------##


isomiRqtl_cis_annot=merge(isomiRqtl_cis_table,snps_information,by.x='snps',by.y='snp.name',all.x=T)
isomiRqtl_cis_annot=merge(isomiRqtl_cis_annot,mir_TSS_annot[Nb_TSS<=1, ],by.x='miRNA',by.y='hsa_ID',all.x=T)

distToRange=function(pos,start,end,strand=NULL,chr=NULL,chrRange=NULL){
        abs_distance_if_outside=pmin(abs(pos-start),abs(pos-end))
        is_before=(pos<start)
        is_after=(end<pos)
        is_inside=start<=pos & pos<=end
        Dist=ifelse( is_inside, 0, abs_distance_if_outside)
        if(!is.null(strand)){
            Dist = Dist * ifelse(is_before & strand=='+' | is_after & strand=='-', -1,1)
        }
        if(!is.null(chr) & !is.null(chrRange)){
            Dist[chr!=chrRange]=Inf
        }
        Dist
}
isomiRqtl_cis_annot[,CisDist_TSS:=distToRange(position, TSS, TSS, mat_strand)/1e6 ]
isomiRqtl_cis_annot[,CisDist_miRNA:=distToRange(position, mat_start, mat_end, mat_strand)/1e6 ]
isomiRqtl_cis_annot[,CisDist_hairpin:=distToRange(position, hairpin_start, hairpin_end, mat_strand)/1e6 ]
isomiRqtl_cis_annot[,Dist_TSS_hairpin:=distToRange(TSS, hairpin_start, hairpin_end, mat_strand)/1e6 ]

isomiRqtl_cis_annot[,type:=ifelse(abs(CisDist_TSS)<.5*abs(Dist_TSS_hairpin),'TSS',
                        ifelse(abs(CisDist_hairpin)<.5*abs(Dist_TSS_hairpin),'Hairpin',
                        ifelse(abs(CisDist_TSS) > .5*abs(Dist_TSS_hairpin) & abs(CisDist_hairpin) > pmax(.5*abs(Dist_TSS_hairpin),.01),'Distant','Undecided')))]


#Namely, mirQTL were first classified as 'mirNA-altering' or 'hairpin-altering' if the overlapped the mature sequence of their associated miRNA or its hairpin.
#Then mirQTL that were less than 10kb from the hairpin or TSS were annotated as 'hairpin-' or 'TSS-flanking' according to which of these feature was the closest.
#Finally, mirQTLs that were located >10kb form both TSS and hairpin were classified as 'Distant'

#
#miRqtl_cis_annot[,typeDetail:=ifelse(abs(CisDist_miRNA)==0,'miRNA-altering',
#                        ifelse(abs(CisDist_hairpin)==0,'Hairpin-altering',
#                        ifelse(abs(CisDist_TSS)<.01 & abs(CisDist_TSS) < abs(CisDist_hairpin),'TSS-flanking',
#                        ifelse(abs(CisDist_hairpin)<.01 & abs(CisDist_hairpin) < abs(CisDist_TSS),'Hairpin-flanking',
#                        ifelse(abs(CisDist_TSS) > .01 & abs(CisDist_hairpin) > .01,'Distant','Undecided')))))]


myDist=.02
isomiRqtl_cis_annot[,typeDetail:=ifelse(abs(CisDist_miRNA)==0,'miRNA-altering',
                        ifelse(abs(CisDist_hairpin)==0,'Hairpin-altering',
                        ifelse(abs(CisDist_TSS) < myDist & abs(CisDist_TSS) < abs(CisDist_hairpin),'TSS-flanking',
                        ifelse(abs(CisDist_hairpin) < myDist & abs(CisDist_hairpin) < abs(CisDist_TSS),'Hairpin-flanking',
                        ifelse(abs(CisDist_TSS) > myDist & abs(CisDist_hairpin) > myDist,'Distant','Undecided')))))]

isomiRqtl_cis_annot[,conserved_TSS:=(Conservation>.2)]


cols=c("snps","chromosome","position", "SNPfreq_AF","SNPfreq_EU", #"ancestral_allele","allele.1","allele.2","daf_char_EUB","daf_char_AFB","FST_adj","iHS_AFB",'iHS_EUB',
            "isomiR","miRNA","assigned_arm",'mat_strand','mat_start','mat_end','CisDist_miRNA',
            "hairpin_name", 'hairpin_start','hairpin_end','CisDist_hairpin',
            "FDR", "beta_NS", "pvalue_NS", "R2_NS",
            "beta_LPS", "pvalue_LPS", "R2_LPS", 
            "beta_PAM3CSK4", "pvalue_PAM3CSK4", "R2_PAM3CSK4", 
            "beta_R848", "pvalue_R848", "R2_R848", 
            "beta_IAV", "pvalue_IAV", "R2_IAV",
            'Promoter.name','TSS','CisDist_TSS','Conservation','RegElt','TFBS','typeDetail')
            
fwrite(isomiRqtl_cis_annot[order(minP),mget(cols)],file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4C_isomirQTLs_Annoted.tsv",EVO_IMMUNO_POP),sep='\t')
isomiRqtl_cis_annot=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable3C_isomirQTLs_Annoted.tsv",EVO_IMMUNO_POP))

#################################
##  load isomiR annotation     ##
#################################

isomir_annot=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_nosubs_FULL.tsv',EVO_IMMUNO_POP))



#######################################
##  compmlete annotation isomiR QTLs ##
#######################################

Map=getMapInfo(isomiRqtl_cis_annot$snps)
isomiRqtl_cis_annot=cbind(isomiRqtl_cis_annot,
                            Map[match(isomiRqtl_cis_annot$snps,Map$snp.name),c('daf_char_AFB','daf_char_EUB'),],
                            isomir_annot[match(isomiRqtl_cis_annot$isomiR,ID),mget(c('hsa_ID','MIMAT','MI_ID','is_cannonical','shift_5p','shift_3p','isomiR_subs','isomir_sequence'))])

enhancers_monocytes = fread(paste(EVO_IMMUNO_POP, "Martin/Project_Neanderthal/ChromHMMFunctional/RegionDefinitions/Regions/Enhancers/enhancers_", "E029", ".tsv", sep=""))
enhancers_monocytes_GR = makeGRangesFromDataFrame(enhancers_monocytes)
enhancers_monocytes_GR= reduce(enhancers_monocytes_GR)
seqlevelsStyle(enhancers_monocytes_GR) <- "NCBI"

promoters_monocytes = fread(paste(EVO_IMMUNO_POP, "Martin/Project_Neanderthal/ChromHMMFunctional/RegionDefinitions/Regions/Promoters/promoters_", "E029", ".tsv", sep=""))
promoters_monocytes_GR = makeGRangesFromDataFrame(promoters_monocytes)
promoters_monocytes_GR=reduce(promoters_monocytes_GR)
seqlevelsStyle(promoters_monocytes_GR) <- "NCBI"

isomiRQTL_GR=makeGRangesFromDataFrame(isomiRqtl_cis_annot[,.(position,chromosome)],start.field='position',end.field='position',seqnames='chromosome')
isomiRqtl_cis_annot[,RegElt_HMM:='']
oo=findOverlaps(isomiRQTL_GR,promoters_monocytes_GR)
isomiRqtl_cis_annot[unique(queryHits(oo)),RegElt_HMM:='Promoter']
oo=findOverlaps(isomiRQTL_GR,enhancers_monocytes_GR)
isomiRqtl_cis_annot[unique(queryHits(oo)),RegElt_HMM:='Enhancer']

cols=c("snps","chromosome","position", 'daf_char_AFB','daf_char_EUB', #"ancestral_allele","allele.1","allele.2","daf_char_EUB","daf_char_AFB","FST_adj","iHS_AFB",'iHS_EUB',
            'hsa_ID','MIMAT','MI_ID','is_cannonical','shift_5p','shift_3p','isomiR_subs','isomir_sequence',
            "FDR", "beta_NS", "pvalue_NS", "R2_NS",
            "beta_LPS", "pvalue_LPS", "R2_LPS", 
            "beta_PAM3CSK4", "pvalue_PAM3CSK4", "R2_PAM3CSK4", 
            "beta_R848", "pvalue_R848", "R2_R848", 
            "beta_IAV", "pvalue_IAV", "R2_IAV","RegElt_HMM",'typeDetail')
isomiRqtl_cis_annot[,isomir_sequence:=gsub('T','U',isomir_sequence)]
isomiRqtl_cis_annot[,isomiR_subs:=gsub('T','U',isomiR_subs)]
fwrite(isomiRqtl_cis_annot[,mget(cols)],file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable3C_isomirQTLs_Annoted_V2.tsv",EVO_IMMUNO_POP),sep='\t')


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




#########################################
##  test condition x snp interactions  ##
#########################################

make_data_frame_to_study <- function(data.table.line){
 snp_name = data.table.line$snps[1]
  miRNA = data.table.line$miRNA[1]
  isomiR = data.table.line$isomiR[1]
  data_temp=isomir_ratio_melted[isomir_ID==isomiR,.(sample,individual,condition,population,ratio_transformed)]
  genotypes = getSNP(snp_name)
  data_temp[, genotype := genotypes[match(data_temp$individual, names(genotypes))]]
  return(data_temp)
}

compute_interaction_effect <- function(data.table.line, cond_to_test){
  data_table_to_study = make_data_frame_to_study(data.table.line)
  snp_name = data.table.line$snps[1]
  miRNA = data.table.line$miRNA[1]
  isomiR = data.table.line$isomiR[1]
  
  data_to_test = data_table_to_study[condition %in% c("NS", cond_to_test)]
  data_to_test[, condition := ifelse(condition == "NS", 0, 1)]
  res = summary(glm(data_to_test, formula = ratio_transformed ~ genotype + condition + population + genotype*condition, family = 'gaussian'))$coefficients
  estimate_interaction_term = res["genotype:condition", "Estimate"]
  pvalue_interaction_term = res["genotype:condition", "Pr(>|t|)"]
  return(data.table(snp_name, miRNA, isomiR, condition = cond_to_test, estimate_interaction_term, pvalue_interaction_term))
}


####################
##Main computation##
####################

results = list()
for (line_n in 1:nrow(isomiRqtl_cis_table)){
  for (cd in c("LPS", "PAM3CSK4", "R848", "IAV")){
    print(paste(line_n, cd))
    results[[paste(line_n, cd)]] = compute_interaction_effect(isomiRqtl_cis_table[line_n], cd)
  }
}
results = rbindlist(results)
results[, pvalue_bonferroni_corrected := p.adjust(pvalue_interaction_term, method = "bonferroni")]
results[, fdr := p.adjust(pvalue_interaction_term, method = "fdr")]

# test = data.table(gene = "hsa-miR-4746-5p", snps = "rs932276", chromosome = "19")
# for (cd in c("LPS", "PAM3", "R848", "IAV")){
#   print(compute_interaction_effect(test, cd))
# }
results= results[order(pvalue_interaction_term),]
results[,minP:=min(pvalue_interaction_term),by=isomiR]
results[,FDR:=min(fdr),by=isomiR]
write.table(results, sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/interaction_genetic_stimulation_on_isomirs.tsv",EVO_IMMUNO_POP), quote = F, row.names = F, sep="\t")
# results=fread(sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/interaction_genetic_stimulation_on_miRNA.tsv",EVO_IMMUNO_POP))
#miRqtl_cis_annot=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4A_mirQTLs_Annoted.tsv",EVO_IMMUNO_POP))

response_isomirQTLs=dcast(results, snp_name + isomiR + miRNA + minP + FDR ~ condition, value.var=c('estimate_interaction_term','pvalue_interaction_term'))
response_cols=c("Beta_response_IAV","Beta_response_LPS","Beta_response_PAM3CSK4","Beta_response_R848","pvalue_response_IAV","pvalue_response_LPS","pvalue_response_PAM3CSK4","pvalue_response_R848")
colnames(response_isomirQTLs)=c("snps","isomiR","miRNA",'minP','FDR',response_cols)

cols=c("snps","chromosome","position", "SNPfreq_AF","SNPfreq_EU", #"ancestral_allele","allele.1","allele.2","daf_char_EUB","daf_char_AFB","FST_adj","iHS_AFB",'iHS_EUB',
            "isomiR","miRNA","assigned_arm",'mat_strand','mat_start','mat_end','CisDist_miRNA',
            "hairpin_name", 'hairpin_start','hairpin_end','CisDist_hairpin',
            'Promoter.name','TSS','CisDist_TSS','Conservation','RegElt','TFBS','typeDetail')

response_isomirQTLs=merge(response_isomirQTLs,isomiRqtl_cis_annot[,mget(cols)],all.x=T,by=c('snps','miRNA','isomiR'))

fwrite(response_isomirQTLs[order(minP),],file=sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/SupTable4D_resp-isomirQTLs_Annoted_raw.tsv",EVO_IMMUNO_POP),sep='\t')










####################################################
##     find best ISOMIR QTL sharing model         ##
####################################################

make_modelisation_joint_isomir <-function(data.table.line){
	temp = make_data_frame_to_study(data.table.line)
	temp[, condition := factor(condition,levels=condIndex)]

	allModels_STIM=cbind(NS=rep(0:1,16),
				LPS=rep(rep(0:1,e=16),1),
				PAM3CSK4=rep(rep(0:1,e=8),2),
				R848=rep(rep(0:1,e=4),4),
				IAV=rep(rep(0:1,e=2),8))
	res=list()
	for (num_model in 1:32){
		model_to_test=allModels_STIM[num_model,]
		temp[, condition_snp:= (genotype)*model_to_test[as.numeric(condition)]]
		snp_name = data.table.line$snps
		miRNA = data.table.line$miRNA
		isomiR = data.table.line$isomiR
		if ( !all(temp$ratio == 1)){#We test if the miRNA is truly expressed in that condition
    		mod = lm(temp, formula = ratio_transformed ~ condition + population + condition_snp)
		    Likelihood = logLik(mod)
		    if( sum(model_to_test)>0 ){
        		Coeffs=summary(mod)$coeff
		        estimate=Coeffs["condition_snp","Estimate"]
        		pvalue=Coeffs["condition_snp","Pr(>|t|)"]
		    }else{
        		estimate=NA
		        pvalue=NA
        	}
		    res[[num_model]] = data.table(miRNA = miRNA, isomiR = isomiR, condition_test = paste(model_to_test,collapse=''), logLik = as.numeric(Likelihood), df = attr(Likelihood,'df'), Beta = estimate, pval = pvalue)
    	}else{
		    res[[num_model]] = data.table(miRNA = character(0), isomiR=character(0), condition_test = character(0), logLik = numeric(0), df = numeric(0), Beta=numeric(0), pval=numeric(0))
   		}
   	}
   	res=rbindlist(res)
   	res[,ProbModel:=exp(logLik-min(logLik))/sum(exp(logLik-min(logLik)))]
   	res
}

make_data_frame_to_study <- function(data.table.line,alleles=F,count=F){
 snp_name = data.table.line$snps[1]
  miRNA = data.table.line$miRNA[1]
  isomiR = data.table.line$isomiR[1]
	if(!count){
	  data_temp=isomir_ratio_melted[isomir_ID==isomiR,.(sample,individual,condition,population,ratio_transformed,ratio)]
	  data_temp[,ratio:=100*ratio]
	  }else{
	  data_temp=isomir_count_melted[isomir_ID==isomiR,.(sample,individual,condition,population,count_transformed,count)]	  
	  }
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

######################
## Main computation ##
######################

results = list()

for (line_n in 1:nrow(isomiRqtl_cis_table)){
    print(paste(line_n))
    results[[paste(line_n)]] = make_modelisation_joint_isomir(isomiRqtl_cis_table[line_n])
  }
results = rbindlist(results)

fwrite(results, sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/sharing_isomirQTL_conditions_likelihoods.tsv",EVO_IMMUNO_POP), sep="\t")
# results=fread(sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/sharing_isomirQTL_conditions_likelihoods.tsv",EVO_IMMUNO_POP),colClass=c('character','character','numeric','numeric','numeric','numeric','numeric'))

bestModel = unique(results[order(isomiR, -logLik)],by= 'isomiR')
bestModel[,bestModel:=condition_test]

response_cols=c("snps","chromosome","position", "SNPfreq_AF","SNPfreq_EU", #"ancestral_allele","allele.1","allele.2","daf_char_EUB","daf_char_AFB","FST_adj","iHS_AFB",'iHS_EUB',
            "isomiR","miRNA","assigned_arm",'RegElt','TFBS','typeDetail',"FDR",
            "Beta_response_LPS","pvalue_response_LPS",
            "Beta_response_PAM3CSK4","pvalue_response_PAM3CSK4",
            "Beta_response_R848","pvalue_response_R848",
            "Beta_response_IAV","pvalue_response_IAV")

response_isomirQTLs_withBM=merge(response_isomirQTLs,bestModel[,mget(c('miRNA','isomiR','bestModel','ProbModel'))],by=c('miRNA','isomiR'))

fwrite(response_isomirQTLs_withBM[order(minP),mget(c(response_cols,'bestModel','ProbModel'))],file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4D_resp-isomirQTLs_Annoted.tsv",EVO_IMMUNO_POP),sep='\t')

response_isomirQTLs_withBM=response_isomirQTLs_withBM[order(minP),mget(c(response_cols,'bestModel','ProbModel'))]

tab_NbCond_isomirQTL=table(sapply(strsplit(response_isomirQTLs_withBM$bestModel,''),function(x){sum(x=='1')}))
tab_NbCond_isomirQTL=tab_NbCond_isomirQTL[as.character(1:5)]
names(tab_NbCond_isomirQTL)=1:5
tab_NbCond_isomirQTL[is.na(tab_NbCond_isomirQTL)]=0
###### piechart of isomirQTL sharing
pieCol=c("#41B6C4","#A1DAB4","#FFFFB2","#FECC5C","#E31A1C")
par(mar=c(3,3,3,3))
tabpct=paste(round(100*tab_NbCond_isomirQTL/sum(tab_NbCond_isomirQTL), 1),'%')
pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/pie_NbCond_isomiR_qtl.pdf",EVO_IMMUNO_POP),width=4,height=4)
pie(tab_NbCond_isomirQTL,col=pieCol,init.angle=90,labels=tabpct) # 5.5 x 5.5 inches
pie(tab_NbCond_isomirQTL,col=pieCol,init.angle=90,labels=rep(' ',5)) # 4 x 4 inches
dev.off()




############################################################
##       plot isomirQTLs across 5 conditions, 1 isomiR    ##
############################################################

plot_isomiR_snp=function(data.table.line,cond=1:5,add.violin=T,count=F,Range=c(0,100)){
    snp_name = data.table.line$snps[1]
    miRNA = data.table.line$miRNA[1]
    isomiR = data.table.line$isomiR[1]
    temp = make_data_frame_to_study(data.table.line,alleles=T,count=count)
    temp = temp[which(condition %in% condIndex[cond]),]
    temp = temp[which(!is.na(genotype)),]
    temp[,condition := factor(condition,levels=condIndex[cond],ordered=TRUE),]
#   temp[,cond_snp := factor(paste(condition,'-',genotype),levels=paste(rep(condIndex,e=2),'-',rep(c('AFB','EUB'),5)))]
    # plot
	if(count){
    	p <- ggplot(temp,aes(x=genotype,y=count,fill=condition))+theme_bw()+facet_grid(~condition)
    }else{
	    p <- ggplot(temp,aes(x=genotype,y=ratio,fill=condition))+theme_bw()+facet_grid(~condition)    
    }
    p <- p + scale_fill_manual(values=colERC5) 
  if(add.violin){
        p <- p + geom_violin(scale='width') + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot(fill="#FFFFFF88", outlier.size=0, notch=TRUE,width=0.4)
    }else{
        p <- p +  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot( outlier.size=0, notch=TRUE)
    }
	  if(!count){
	  p <- p + ylim(Range[1],Range[2]) + xlab('') + ylab('Percentage of miRNA reads')
	  }else{
	  p <- p + ylab('log2(RPM)') + xlab('')
	  }
    p <- p + geom_jitter(width=0.2,size=0.5,colour="#52525266")
    print(p)
}

#################################
##     plot isomiR QTL 1 by 1     ##
#################################
isomiRqtl_cis_table=isomiRqtl_cis_table[order(minP)]

regexpr_isomiR='chr([0-9XYMT]+)_([-0-9]+)_([-0-9]+)_([+-])_([ATGC]+)_[0-2];(.*)_(hsa.*)_(MIMAT.*)_(MI[0-9]+)'
regexpr_other='(hsa.*)_(MIMAT.*)_(MI[0-9]+)_other'

for (line_n in 1:nrow(isomiRqtl_cis_table)){
    print(paste(line_n))
    mymiR=isomiRqtl_cis_table[line_n,miRNA]
    myisomiR=isomiRqtl_cis_table[line_n,isomiR]
    mysnp=isomiRqtl_cis_table[line_n,snps]
	isomir_ID_simple=ifelse(!grepl('other',myisomiR),gsub(regexpr_isomiR,'\\7_\\2;\\3;\\6',myisomiR),gsub(regexpr_other,'\\1_other',myisomiR))

	pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/exemples/isomiRs/%s_%s_%s_%s-ALL.pdf",EVO_IMMUNO_POP,line_n,mymiR,isomir_ID_simple,mysnp),height=2.8,width=7)
	plot_isomiR_snp(isomiRqtl_cis_table[line_n],1:5)
	dev.off()
	pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/exemples/isomiRs/%s_%s_%s_%s-ALL_count.pdf",EVO_IMMUNO_POP,line_n,mymiR,isomir_ID_simple,mysnp),height=2.8,width=7)
	plot_isomiR_snp(isomiRqtl_cis_table[line_n],1:5,count=T)
	dev.off()
}

line_n=1
    print(paste(line_n))
    mymiR=response_isomirQTLs_withBM[line_n,miRNA]
    myisomiR=response_isomirQTLs_withBM[line_n,isomiR]
    mysnp=response_isomirQTLs_withBM[line_n,snps]
	isomir_ID_simple=ifelse(!grepl('other',myisomiR),gsub(regexpr_isomiR,'\\7_\\2;\\3;\\6',myisomiR),gsub(regexpr_other,'\\1_other',myisomiR))

	pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/exemples/isomiRs/%s_%s_%s_%s-ALL_rescaled.pdf",EVO_IMMUNO_POP,line_n,mymiR,isomir_ID_simple,mysnp),height=2.8,width=7)
	plot_isomiR_snp(response_isomirQTLs_withBM[line_n],1:5,Range=c(6,24))
	dev.off()


############################################################
##     plot isomirQTLs across 5 conditions, all isomiRs   ##
############################################################

make_data_frame_to_study_allisomiRs <- function(data.table.line,alleles=F,count=F){
  snp_name = data.table.line$snps[1]
  miRNA = data.table.line$miRNA[1]
	if(!count){
	  data_temp=isomir_ratio_melted[hsa_ID==miRNA,.(sample,isomir_ID,individual,condition,population,ratio_transformed,ratio)]
	  data_temp[,ratio:=100*ratio]
	  }else{
	  data_temp=isomir_count_melted[hsa_ID==miRNA,.(sample,isomir_ID,individual,condition,population,count_transformed,count)]	  
	  }
    regexpr_isomiR='chr([0-9XYMT]+)_([-0-9]+)_([-0-9]+)_([+-])_([ATGC]+)_[0-2];(.*)_(hsa.*)_(MIMAT.*)_(MI[0-9]+)'
    regexpr_other='(hsa.*)_(MIMAT.*)_(MI[0-9]+)_other'
    data_temp[,isomir_ID_simple:=ifelse(!grepl('other',isomir_ID),gsub(regexpr_isomiR,'\\7_\\2;\\3;\\6',isomir_ID),gsub(regexpr_other,'\\1_other',isomir_ID))]

  map=getMapInfo(snp_name)
  genotypes=getSNP(snp_name)
  genoNum=2-genotypes
  data_temp[, genoNum := genoNum[match(data_temp$individual, names(genotypes))]]

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
    }
    data_temp
}


plot_all_isomiRs_snp=function(Data,cond,add.violin=T,count=F){
    temp = Data
    temp = temp[which(condition %in% condIndex[cond]),]
    temp = temp[which(!is.na(genotype)),]
    temp[,condition := factor(condition,levels=condIndex[cond],ordered=TRUE),]
    # plot
	if(count){
    	p <- ggplot(temp,aes(x=genotype,y=count,fill=condition))+theme_bw()+facet_grid(condition~factor(isomir_ID_simple))
    }else{
	    p <- ggplot(temp,aes(x=genotype,y=ratio,fill=condition))+theme_bw()+facet_grid(condition~factor(isomir_ID_simple))    
    }
    p <- p + scale_fill_manual(values=colERC5) 
  if(add.violin){
        p <- p + geom_violin(scale='width') + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot(fill="#FFFFFF88", outlier.size=0, notch=TRUE,width=0.4)
    }else{
        p <- p +  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot( outlier.size=0, notch=TRUE)
    }
	  if(!count){
	  p <- p + ylim(0,100) + xlab('') + ylab('Percentage of miRNA reads')
	  }else{
	  p <- p + ylab('log2(RPM)') + xlab('')
	  }
    p <- p + geom_jitter(width=0.2,size=0.5,colour="#52525266")
    print(p)
}




for (line_n in 1:nrow(unique(isomiRqtl_cis_table,by='miRNA'))){
	myline = unique(isomiRqtl_cis_table,by='miRNA')[line_n]
	mymiR=myline[,miRNA]
    mysnp=myline[,snps]
    Data=make_data_frame_to_study_allisomiRs(myline, allele=T)
    Data_count=make_data_frame_to_study_allisomiRs(myline, allele=T,count=T)
    for (cond in 1:5){
		pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/exemples/all_isomiRs/%s_%s_%s-cond%s.pdf",EVO_IMMUNO_POP,line_n,mymiR,mysnp,cond),height=2.8,width=7)
		plot_all_isomiRs_snp(Data, cond)
		dev.off()
		pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/exemples/all_isomiRs/%s_%s_%s-cond%s_count.pdf",EVO_IMMUNO_POP,line_n,mymiR,mysnp,cond),height=2.8,width=7)
		plot_all_isomiRs_snp(Data_count,cond,count=T)
		dev.off()
	}
	pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/exemples/all_isomiRs/%s_%s_%s-ALL_count.pdf",EVO_IMMUNO_POP,line_n,mymiR,mysnp),height=2.8,width=7)
	plot_all_isomiRs_snp(Data_count,1:5,count=T)
	dev.off()
}


###############################################################################
## 			      infos on isomiR QTLs of miR-146a-3p/5p  					 ##
###############################################################################
# (and their targets)

