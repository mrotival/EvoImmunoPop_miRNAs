
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
miRqtl_cis = fread(sprintf("%s/Maxime/miRNA_V2/data/09.mirQTL/cis_mirQTLs_with_FDR_filtered_best_miRNA_snp_association.tsv",EVO_IMMUNO_POP))
miRqtl_cis[,miRNA:=gene]
miRqtl_cis[,minP:=min(pvalue),by=miRNA]
miRqtl_cis[,FDR:=min(FDR_Perm),by=miRNA]
miRqtl_cis[,maxBeta:=max(abs(beta)),by=miRNA]
miRqtl_cis[,df_cond:=df_cond[condition]-3]
miRqtl_cis[,R2:=statistic^2/(statistic^2+df_cond)]
miRqtl_cis[,condition:=factor(condIndex[condition],levels=condIndex)]
cols=c("snps", "miRNA", "FDR", "beta_NS", "pvalue_NS", "R2_NS",
                          "beta_LPS", "pvalue_LPS", "R2_LPS", 
                          "beta_PAM3CSK4", "pvalue_PAM3CSK4", "R2_PAM3CSK4", 
                          "beta_R848", "pvalue_R848", "R2_R848", 
                          "beta_IAV", "pvalue_IAV", "R2_IAV")

miRqtl_cis_table=dcast(miRqtl_cis,snps + miRNA + minP + maxBeta + FDR ~ condition, value.var=c('beta', 'pvalue', 'R2'))
fwrite(miRqtl_cis_table[order(minP),mget(cols)],file=sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/SupTable4A_mirQTLs_raw.tsv",EVO_IMMUNO_POP),sep='\t')
#miRqtl_cis_table=fread(sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/SupTable4A_mirQTLs_raw.tsv",EVO_IMMUNO_POP))
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

snps_information=getMapInfo(unique(miRqtl_cis_table$snps))[c("snp.name","chromosome","position","allele.2","allele.1","SNPfreq_AF","SNPfreq_EU","daf_char_AFB","daf_char_EUB","RegElt","TFBS", 'FST_adj','iHS_AFB','iHS_EUB','P_CLS2_EUB','P_CLS2_AFB','aSNP')]
fwrite(snps_information,file=sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/SupTable4A_mirQTLs_snp_Annot.tsv",EVO_IMMUNO_POP),sep='\t')
snps_information=fread(sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/SupTable4A_mirQTLs_snp_Annot.tsv",EVO_IMMUNO_POP))

mir_TSS_annot=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_TSS_annot_DeRie2017.txt',EVO_IMMUNO_POP))
##--------------------------------------------##
## private version (based on mysql database)  ##
##--------------------------------------------##


miRqtl_cis_annot=merge(miRqtl_cis_table,snps_information,by.x='snps',by.y='snp.name',all.x=T)
miRqtl_cis_annot=merge(miRqtl_cis_annot,mir_TSS_annot[Nb_TSS<=1, ],by.x='miRNA',by.y='hsa_ID',all.x=T)

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
miRqtl_cis_annot[,CisDist_TSS:=distToRange(position, TSS, TSS, mat_strand)/1e6 ]
miRqtl_cis_annot[,CisDist_miRNA:=distToRange(position, mat_start, mat_end, mat_strand)/1e6 ]
miRqtl_cis_annot[,CisDist_hairpin:=distToRange(position, hairpin_start, hairpin_end, mat_strand)/1e6 ]
miRqtl_cis_annot[,Dist_TSS_hairpin:=distToRange(TSS, hairpin_start, hairpin_end, mat_strand)/1e6 ]

miRqtl_cis_annot[,type:=ifelse(abs(CisDist_TSS)<.5*abs(Dist_TSS_hairpin),'TSS',
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
miRqtl_cis_annot[,typeDetail:=ifelse(abs(CisDist_miRNA)==0,'miRNA-altering',
                        ifelse(abs(CisDist_hairpin)==0,'Hairpin-altering',
                        ifelse(abs(CisDist_TSS) < myDist & abs(CisDist_TSS) < abs(CisDist_hairpin),'TSS-flanking',
                        ifelse(abs(CisDist_hairpin) < myDist & abs(CisDist_hairpin) < abs(CisDist_TSS),'Hairpin-flanking',
                        ifelse(abs(CisDist_TSS) > myDist & abs(CisDist_hairpin) > myDist,'Distant','Undecided')))))]

# mirQTL for which the TSS of the pri-miRNA were annotated based on their location relative to the miRNA.
# Namely, mirQTL were first classified as 'mirNA-altering' or 'hairpin-altering' if the overlapped the sequence of their associated (mature) miRNA or its hairpin.
# Then for each mirQTL we computed the distance between the SNP and both the hairpin and the TSS of the pri-miRNA and annotated 
# mirQTL for which the that were less than 10kb from the hairpin or TSS were annotated as 'hairpin-' or 'TSS-flanking' according to which of these feature was the closest.
#Finally, mirQTLs that were located >10kb form both TSS and hairpin were classified as 'Distant'

#miRqtl_cis_annot[,typeDetail_2:=ifelse(abs(CisDist_TSS)<.5*abs(Dist_TSS_hairpin),'TSS',
#                       ifelse(abs(CisDist_miRNA)==0,'miRNA',
#                       ifelse(abs(CisDist_hairpin)==0,'Hairpin',
#                       ifelse(abs(CisDist_hairpin)<.5*abs(Dist_TSS_hairpin),'Hairpin flanks',
#                       ifelse(abs(CisDist_TSS) > .5*abs(Dist_TSS_hairpin) & abs(CisDist_hairpin) > pmax(.5*abs(Dist_TSS_hairpin),.01),'Distant','Undecided')))))]
#

# mirQTL for which the distance between the were less than 10kb from the hairpin or TSS were annotated 
miRqtl_cis_annot[,conserved_TSS:=(Conservation>.2)]


cols=c("snps","chromosome","position", "SNPfreq_AF","SNPfreq_EU", #"ancestral_allele","allele.1","allele.2","daf_char_EUB","daf_char_AFB","FST_adj","iHS_AFB",'iHS_EUB',
            "miRNA","assigned_arm",'mat_strand','mat_start','mat_end','CisDist_miRNA',
            "hairpin_name", 'hairpin_start','hairpin_end','CisDist_hairpin',
            "FDR", "beta_NS", "pvalue_NS", "R2_NS",
            "beta_LPS", "pvalue_LPS", "R2_LPS", 
            "beta_PAM3CSK4", "pvalue_PAM3CSK4", "R2_PAM3CSK4", 
            "beta_R848", "pvalue_R848", "R2_R848", 
            "beta_IAV", "pvalue_IAV", "R2_IAV",
            'Promoter.name','TSS','CisDist_TSS','Conservation','RegElt','TFBS','typeDetail')

fwrite(miRqtl_cis_annot[order(minP),mget(cols)],file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4A_mirQTLs_Annoted.tsv",EVO_IMMUNO_POP),sep='\t')
# miRqtl_cis_annot=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4A_mirQTLs_Annoted.tsv",EVO_IMMUNO_POP))

#########################################################
###     add infos on DAF, HMM and detail miR IDs      ###
#########################################################

cols=c("snps","chromosome","position", "daf_char_EUB","daf_char_AFB", "miRNA","MIMAT","MI_ID",
            "FDR", "beta_NS", "pvalue_NS", "R2_NS",
            "beta_LPS", "pvalue_LPS", "R2_LPS", 
            "beta_PAM3CSK4", "pvalue_PAM3CSK4", "R2_PAM3CSK4", 
            "beta_R848", "pvalue_R848", "R2_R848", 
            "beta_IAV", "pvalue_IAV", "R2_IAV",
            'typeDetail','RegElt_HMM')
            

enhancers_monocytes = fread(paste(EVO_IMMUNO_POP, "Martin/Project_Neanderthal/ChromHMMFunctional/RegionDefinitions/Regions/Enhancers/enhancers_", "E029", ".tsv", sep=""))
enhancers_monocytes_GR = makeGRangesFromDataFrame(enhancers_monocytes)
enhancers_monocytes_GR= reduce(enhancers_monocytes_GR)
seqlevelsStyle(enhancers_monocytes_GR) <- "NCBI"

promoters_monocytes = fread(paste(EVO_IMMUNO_POP, "Martin/Project_Neanderthal/ChromHMMFunctional/RegionDefinitions/Regions/Promoters/promoters_", "E029", ".tsv", sep=""))
promoters_monocytes_GR = makeGRangesFromDataFrame(promoters_monocytes)
promoters_monocytes_GR=reduce(promoters_monocytes_GR)
seqlevelsStyle(promoters_monocytes_GR) <- "NCBI"

miRQTL_GR=makeGRangesFromDataFrame(miRqtl_cis_annot[,.(position,chromosome)],start.field='position',end.field='position',seqnames='chromosome')
miRqtl_cis_annot[,RegElt_HMM:='']
oo=findOverlaps(miRQTL_GR,promoters_monocytes_GR)
miRqtl_cis_annot[unique(queryHits(oo)),RegElt_HMM:='Promoter']
oo=findOverlaps(miRQTL_GR,enhancers_monocytes_GR)
miRqtl_cis_annot[unique(queryHits(oo)),RegElt_HMM:='Enhancer']
fwrite(miRqtl_cis_annot[order(minP),mget(cols)],file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable3A_mirQTLs_Annoted_V2.tsv",EVO_IMMUNO_POP),sep='\t')

#######################
##Necessary functions##
#######################

# We need, given a SNP and a miRNA to get
# The genotypes of each individual
# The covariates of each individual (population)
# The expression of the miRNA for each individual
# We want to order this so that each sample (individual*condition) is one line

#################################
##    load expression Data     ##
#################################

print("loading miRNA data")
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


make_data_frame_to_study <- function(data.table.line){
 snp_name = data.table.line$snps[1]
  miRNA = data.table.line$miRNA[1]
  data_temp=mir_count_melted[ID==miRNA,.(ID,individual,condition,population,count_transformed)]
  genotypes = getSNP(snp_name)
  data_temp[, genotype := genotypes[match(data_temp$individual, names(genotypes))]]
  return(data_temp)
}

compute_interaction_effect <- function(data.table.line, cond_to_test){
  data_table_to_study = make_data_frame_to_study(data.table.line)
  snp_name = data.table.line$snps[1]
  miRNA = data.table.line$miRNA[1]
  data_to_test = data_table_to_study[condition %in% c("NS", cond_to_test)]
  data_to_test[, condition := ifelse(condition == "NS", 0, 1)]
  res = summary(glm(data_to_test, formula = count_transformed ~ genotype + condition + population + genotype*condition, family = 'gaussian'))$coefficients
  estimate_interaction_term = res["genotype:condition", "Estimate"]
  pvalue_interaction_term = res["genotype:condition", "Pr(>|t|)"]
  return(data.table(snp_name, miRNA, condition = cond_to_test, estimate_interaction_term, pvalue_interaction_term))
}

####################
##Main computation##
####################

results = list()
for (line_n in 1:nrow(miRqtl_cis_table)){
  for (cd in c("LPS", "PAM3CSK4", "R848", "IAV")){
    print(paste(line_n, cd))
    results[[paste(line_n, cd)]] = compute_interaction_effect(miRqtl_cis_table[line_n], cd)
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
results[,minP:=min(pvalue_interaction_term),by=miRNA]
results[,FDR:=min(fdr),by=miRNA]
write.table(results, sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/interaction_genetic_stimulation_on_miRNA.tsv",EVO_IMMUNO_POP), quote = F, row.names = F, sep="\t")
# results=fread(sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/interaction_genetic_stimulation_on_miRNA.tsv",EVO_IMMUNO_POP))
#miRqtl_cis_annot=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4A_mirQTLs_Annoted.tsv",EVO_IMMUNO_POP))

response_mirQTLs=dcast(results, snp_name + miRNA + minP + FDR ~ condition, value.var=c('estimate_interaction_term','pvalue_interaction_term'))
response_cols=c("Beta_response_IAV","Beta_response_LPS","Beta_response_PAM3CSK4","Beta_response_R848","pvalue_response_IAV","pvalue_response_LPS","pvalue_response_PAM3CSK4","pvalue_response_R848")
colnames(response_mirQTLs)=c("snps","miRNA",'minP','FDR',response_cols)

cols=c("snps","chromosome","position", "SNPfreq_AF","SNPfreq_EU", #"ancestral_allele","allele.1","allele.2","daf_char_EUB","daf_char_AFB","FST_adj","iHS_AFB",'iHS_EUB',
            "miRNA","assigned_arm",'mat_strand','mat_start','mat_end','CisDist_miRNA',
            "hairpin_name", 'hairpin_start','hairpin_end','CisDist_hairpin',
            'Promoter.name','TSS','CisDist_TSS','Conservation','RegElt','TFBS','typeDetail')

response_mirQTLs=merge(response_mirQTLs,miRqtl_cis_annot[,mget(cols)],all.x=T,by=c('snps','miRNA'))

fwrite(response_mirQTLs[order(minP),],file=sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/SupTable4B_resp-mirQTLs_Annoted_raw.tsv",EVO_IMMUNO_POP),sep='\t')



####################################################
##	       find best QTL sharing model         	  ##
####################################################

make_modelisation_joint_mir <-function(data.table.line){
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
		if ( !all(temp$count == 0)){#We test if the miRNA is truly expressed in that condition
    		mod = lm(temp, formula = count_transformed ~ condition + population + condition_snp)
		    Likelihood = logLik(mod)
		    if( sum(model_to_test)>0 ){
        		Coeffs=summary(mod)$coeff
		        estimate=Coeffs["condition_snp","Estimate"]
        		pvalue=Coeffs["condition_snp","Pr(>|t|)"]
		    }else{
        		estimate=NA
		        pvalue=NA
        	}
		    res[[num_model]] = data.table(miRNA = miRNA, condition_test = paste(model_to_test,collapse=''), logLik = as.numeric(Likelihood), df = attr(Likelihood,'df'), Beta = estimate, pval = pvalue)
    	}else{
		    res[[num_model]] = data.table(miRNA = character(0), condition_test = character(0), logLik = numeric(0), df = numeric(0), Beta=numeric(0), pval=numeric(0))
   		}
   	}
   	res=rbindlist(res)
   	res[,ProbModel:=exp(logLik-min(logLik))/sum(exp(logLik-min(logLik)))]
   	res
}

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

######################
## Main computation ##
######################

results = list()

for (line_n in 1:nrow(miRqtl_cis_table)){
    print(paste(line_n))
    results[[paste(line_n)]] = make_modelisation_joint_mir(miRqtl_cis_table[line_n])
  }
results = rbindlist(results)

fwrite(results, sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/sharing_mirQTL_conditions_likelihoods.tsv",EVO_IMMUNO_POP), sep="\t")
# results=fread(sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/sharing_mirQTL_conditions_likelihoods.tsv",EVO_IMMUNO_POP),colClass=c('character','character','numeric','numeric','numeric','numeric','numeric'))

bestModel = unique(results[order(miRNA, -logLik)],by= 'miRNA')
bestModel[,bestModel:=condition_test]

response_cols=c("snps","chromosome","position", "SNPfreq_AF","SNPfreq_EU", #"ancestral_allele","allele.1","allele.2","daf_char_EUB","daf_char_AFB","FST_adj","iHS_AFB",'iHS_EUB',
            "miRNA","assigned_arm",'RegElt','TFBS','typeDetail',"FDR",
            "Beta_response_LPS","pvalue_response_LPS",
            "Beta_response_PAM3CSK4","pvalue_response_PAM3CSK4",
            "Beta_response_R848","pvalue_response_R848",
            "Beta_response_IAV","pvalue_response_IAV")

response_mirQTLs_withBM=merge(response_mirQTLs,bestModel[,mget(c('miRNA','bestModel','ProbModel'))],by='miRNA')

fwrite(response_mirQTLs_withBM[order(minP),mget(c(response_cols,'bestModel','ProbModel'))],file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4B_resp-mirQTLs_Annoted.tsv",EVO_IMMUNO_POP),sep='\t')

response_mirQTLs_withBM=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable3B_resp-mirQTLs_Annoted.tsv",EVO_IMMUNO_POP))
response_mirQTLs_withBM=cbind(response_mirQTLs_withBM,miRqtl_cis_annot[match(response_mirQTLs_withBM$snps,snps),mget(c("daf_char_EUB","daf_char_AFB", "miRNA",'arm',"MIMAT","MI_ID",'typeDetail','RegElt_HMM'))])

response_cols=c("snps","chromosome","position", "daf_char_EUB","daf_char_AFB", 
            "miRNA","miRNA","MIMAT","MI_ID","FDR",
            "Beta_response_LPS","pvalue_response_LPS",
            "Beta_response_PAM3CSK4","pvalue_response_PAM3CSK4",
            "Beta_response_R848","pvalue_response_R848",
            "Beta_response_IAV","pvalue_response_IAV",'bestModel','ProbModel')
fwrite(response_mirQTLs_withBM[,mget(response_cols)],file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable3B_resp-mirQTLs_Annoted_V2.tsv",EVO_IMMUNO_POP))

tab_NbCond_mirQTL=table(sapply(strsplit(response_mirQTLs_withBM$bestModel,''),function(x){sum(x=='1')}))
tab_NbCond_mirQTL=tab_NbCond_mirQTL[as.character(1:5)]
names(tab_NbCond_mirQTL)=1:5
tab_NbCond_mirQTL[is.na(tab_NbCond_mirQTL)]=0
###### piechart of mirQTL sharing
pieCol=c("#41B6C4","#A1DAB4","#FFFFB2","#FECC5C","#E31A1C")
par(mar=c(3,3,3,3))
tabpct=paste(round(100*tab_NbCond_mirQTL/sum(tab_NbCond_mirQTL), 1),'%')
pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/pie_NbCond_miR_qtl.pdf",EVO_IMMUNO_POP),width=4,height=4)
pie(tab_NbCond_mirQTL,col=pieCol,init.angle=90,labels=tabpct) # 5.5 x 5.5 inches
pie(tab_NbCond_mirQTL,col=pieCol,init.angle=90,labels=rep(' ',5)) # 4 x 4 inches
dev.off()

####################################################
##       plot mirQTLs across 5 conditions         ##
####################################################

plot_miRNA_snp=function(data.table.line,cond=1:5,add.violin=T,transformed=FALSE,pop='.*'){
    snp_name = data.table.line$snps[1]
    miRNA = data.table.line$miRNA[1]
    temp = make_data_frame_to_study(data.table.line,alleles=T)
    temp = temp[which(condition %in% condIndex[cond]),]
    temp = temp[which(!is.na(genotype)),]
    temp = temp[grep(pop,population)]
    temp[,condition := factor(condition,levels=condIndex[cond],ordered=TRUE),]
#   temp[,cond_snp := factor(paste(condition,'-',genotype),levels=paste(rep(condIndex,e=2),'-',rep(c('AFB','EUB'),5)))]
    # plot
    if(!transformed){
    p <- ggplot(temp,aes(x=genotype,y=count,fill=condition))+theme_bw()+facet_grid(~condition)
    }else{
    p <- ggplot(temp,aes(x=genotype,y=count_transformed,fill=condition))+theme_bw()+facet_grid(~condition)
    }
    p <- p + scale_fill_manual(values=colERC5) 
  if(add.violin){
        p <- p + geom_violin(scale='width') + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot(fill="#FFFFFF88", outlier.size=0, notch=TRUE,width=0.4)
    }else{
        p <- p +  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot( outlier.size=0, notch=TRUE)
    }
    p <- p + geom_jitter(width=0.2,size=0.5,colour="#52525266")
    print(p)
}

#################################
##     plot miR QTL 1 by 1     ##
#################################
miRqtl_cis_table=miRqtl_cis_table[order(minP)]
for (line_n in 1:nrow(miRqtl_cis_table)){
    print(paste(line_n))
    mymiR=miRqtl_cis_table[line_n,miRNA]
    mysnp=miRqtl_cis_table[line_n,snps]
	pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/exemples/%s_%s_%s-ALL.pdf",EVO_IMMUNO_POP,line_n,mymiR,mysnp),height=2.8,width=7)
	plot_miRNA_snp(miRqtl_cis_table[line_n],1:5)
	dev.off()
}


#####################################
##     plot miR QTL pop by pop     ##
#####################################

plot_miRNA_snp_pop=function(data.table.line,cond=1:5,add.violin=T,transformed=FALSE){
    snp_name = data.table.line$snps[1]
    miRNA = data.table.line$miRNA[1]
    temp = make_data_frame_to_study(data.table.line,alleles=T)
    temp = temp[which(condition %in% condIndex[cond]),]
    temp = temp[which(!is.na(genotype)),]
    temp[,condition := factor(condition,levels=condIndex[cond],ordered=TRUE),]
#   temp[,cond_snp := factor(paste(condition,'-',genotype),levels=paste(rep(condIndex,e=2),'-',rep(c('AFB','EUB'),5)))]
    # plot
    if(!transformed){
    p <- ggplot(temp,aes(x=genotype,y=count,fill=condition))+theme_bw()+facet_grid(condition~population)
    }else{
    p <- ggplot(temp,aes(x=genotype,y=count_transformed,fill=condition))+theme_bw()+facet_grid(condition~population)
    }
    p <- p + scale_fill_manual(values=colERC5) 
  if(add.violin){
        p <- p + geom_violin(scale='width') + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot(fill="#FFFFFF88", outlier.size=0, notch=TRUE,width=0.4)
    }else{
        p <- p +  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot( outlier.size=0, notch=TRUE)
    }
    p <- p + geom_jitter(width=0.2,size=0.5,colour="#52525266")
    print(p)
}


#####################################
##     plot miR QTL rs127881760     ##
#####################################
line_n=which(miRqtl_cis_table[,snps]=='rs12881760')[1]
   mymiR=miRqtl_cis_table[line_n,miRNA]
   mysnp=miRqtl_cis_table[line_n,snps]

	pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/exemples/%s_%s_%s-NS_both.pdf",EVO_IMMUNO_POP,line_n,mymiR,mysnp),height=2.8,width=4)
	plot_miRNA_snp_pop(miRqtl_cis_table[line_n],1)
	dev.off()

#####################################
##     plot miR QTL rs5743168     ##
#####################################


###############################################
##     plot miR QTL miR-146a-3p     ##
###############################################

miRqtl_cis_table=miRqtl_cis_table[order(minP)]
line_n=122
    print(paste(line_n))
    mymiR=miRqtl_cis_table[line_n,miRNA]
    mysnp=miRqtl_cis_table[line_n,snps]
	pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/response_miRQTL_miR-146a-5p-AFB.pdf",EVO_IMMUNO_POP,line_n,mymiR,mysnp),height=2.8,width=7)
	plot_miRNA_snp(miRqtl_cis_table[line_n],1:5,transformed=FALSE,pop='AFB')
	dev.off()
	pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/response_miRQTL_miR-146a-5p-EUB.pdf",EVO_IMMUNO_POP,line_n,mymiR,mysnp),height=2.8,width=7)
	plot_miRNA_snp(miRqtl_cis_table[line_n],1:5,transformed=FALSE,pop='EUB')
	dev.off()

################################
##       DATA LOADING         ##
##  version for external use  ##
################################

#snps_informations = list()
#for (c in 1:22){
# print(c)
# snps_informations[[c]] = fread(paste("../../06.snp_data/data/general_informations/snps_chr", c, ".tsv", sep=""))
#}
#rm(c)
#snps_informations =rbindlist(snps_informations)

#get_covariates <- function(){
#  cov_inf = fread(paste(EVO_IMMUNO_POP, "/Maxime/Evo_Immuno_pop_data/GenotypePCs/EvoImmunoPop_Behar_GabonDiv_ordered_ADMIXTURE.tsv", sep=""), select = c("ID", "K=2"))
#  cov_inf = cov_inf[, 1:2, with = F]
#  names(cov_inf) = c("individual", "covar")
#  cov_inf = cov_inf[grepl("AFB", individual) | grepl("EUB", individual)]
#  return(cov_inf)
#}

#covariates = get_covariates()

#get_genotypes <- function(chr, snp_name){
#  all_genos = fread(paste("../../06.snp_data/data/genotypes/genotype_matrix_chr", chr, ".tsv", sep=""))
#  interesting_genos = melt(all_genos[snp_id == snp_name], id = "snp_id", variable.name = "individual", value.name = "genotypes")
#  return(interesting_genos)
#}

#
#get_expression <- function(condition, miRNA){
#  expression = fread(paste("../../07.eQTL/FilesForComputations/expression1_", condition, ".tsv", sep=""))
#  expression = melt(expression[gene_id == miRNA],id = c("gene_id"), variable.name = "individual", value.name ="RPM"  )
#  expression = expression[, condition := condition]
#  return(expression)
#}
#a = get_expression(5, "hsa-miR-4632-3p")
#
#make_data_frame_to_study <-function(data.table.line){
#  condition_name = c("1"="NS", "2"="LPS", "3" = "PAM3", "4" ="R848", "5" = "IAV")
#  snp_name = data.table.line$snps[1]
#  miRNA = data.table.line$gene[1]
#  chromosome = data.table.line$chromosome[1]
#  data_temp = rbindlist(lapply(1:5, get_expression, miRNA = miRNA))
#  genotypes = get_genotypes(chromosome, snp_name)
#  data_temp[, genotype := genotypes[match(data_temp$individual, individual), genotypes]]
#  covariates_db = get_covariates()
#  data_temp[, covariate := covariates_db[match(data_temp$individual, individual), covar]]
#  data_temp[, condition := condition_name[as.character(condition)]]
#  data_temp[, population := substr(individual, 1,3)]
#  print(data_temp)
#
#  return(data_temp)
#}
