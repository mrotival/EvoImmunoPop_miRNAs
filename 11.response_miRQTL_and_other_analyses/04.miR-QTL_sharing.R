


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
		gene = data.table.line$gene
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



response_mirQTLs_withBM=response_mirQTLs_withBM[order(minP),mget(c(response_cols,'bestModel','ProbModel'))]
response_mirQTLs_withBM=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4B_resp-mirQTLs_Annoted.tsv",EVO_IMMUNO_POP))


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



