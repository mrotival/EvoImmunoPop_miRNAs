##################################################################################
##Aim of the script : investigate which isomirs change the *ratio* of expression##
##Upon stimulation                                                              ##
##################################################################################

##############
##Librairies##
##############

condIndex=c("NS","LPS","PAM3CSK4","R848","IAV")
luq = function(x){length(unique(x))}

###############
## Load data ##
###############

isomir_ratio_cast=fread(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_ratios_aggregated_nosubs.GCRL_Batch_lane_corrected.tsv",EVO_IMMUNO_POP))
isomir_ratio_cast=isomir_ratio_cast[-grep('other',isomir_ID)]
dim(isomir_ratio_cast)
isomir_ratio_melted = melt(isomir_ratio_cast, id = c("mirID", "isomir_ID"),
                               variable.name = "sample",
                               value.name = "ratio")
isomir_ratio_melted[, individual := substr(sample, 1,6)]
isomir_ratio_melted[, condition := condIndex[as.numeric(substr(sample, 8,8))]]
isomir_ratio_melted[, population := substr(sample, 1,3)]
inverseNormalRankTransform=function(x){n=length(x);
									qnorm(frank(x,na.last=FALSE,ties.method='random')/(n+1),mean(x,na.rm=T),sd(x,na.rm=T))}
isomir_ratio_melted[, ratio_transformed := inverseNormalRankTransform(ratio),by=.(mirID,isomir_ID)]

######################################
##    Convert into response data    ##
######################################

isomir_response = dcast(isomir_ratio_melted, mirID + isomir_ID + individual + population ~ condition, value.var=list('ratio_transformed','ratio'))
isomir_response_melted=melt(isomir_response,id.vars=c('mirID','isomir_ID','individual','population','ratio_NS','ratio_transformed_NS'), measure.var=list(c('ratio_transformed_LPS','ratio_transformed_PAM3CSK4','ratio_transformed_R848','ratio_transformed_IAV'),c('ratio_LPS','ratio_PAM3CSK4','ratio_R848','ratio_IAV')),value.name=c('ratio_transformed','ratio'),variable.name='stimulus')
isomir_response_melted[,condition := condIndex[as.numeric(stimulus)+1]]
isomir_response_melted[,response := ratio-ratio_NS]
isomir_response_melted[,response_transformed := ratio_transformed-ratio_transformed_NS]


#####################################################
##Let's remove miRNA with several possible hairpins##
#####################################################
# extract infos from hairpin ID

# miRNA_with_multiple_hairpins = isomir_counts_melted[, .(numbermir2 = length(unique(mir2)), numberhairpin = length(unique(hairpin))), by = mir1]
# miRNA_with_multiple_hairpins = miRNA_with_multiple_hairpins[numberhairpin !=1]
# isomir_counts_melted = isomir_counts_melted[!(mir1 %in% miRNA_with_multiple_hairpins[, mir1])]

##################################
##Transform count in proportions##
##################################
# isomir_counts_melted[, ratio_expression := count/sum(count), by = c("mir1", "condition", "individual")]
# write.table(isomir_counts_melted, "../../09.isomirs/data/isomirs_ratios.tsv", sep="\t", quote = F, row.names = F)
## We will probably have to group on more precise cases of isomirs to do this analysis

##################################################################
##     Test for differential expression upon stimulation        ##
##################################################################

condIndex=c("NS","LPS","PAM3CSK4","R848","IAV")

library(lme4)
library(lmerTest)




# function 
make_modelisation <-function(isomir_ID_to_test, condition_to_test){
  temp = isomir_ratio_melted[(isomir_ID == isomir_ID_to_test) & (condition %in% c(condition_to_test))]
  MeanAFB = mean(temp[population=='AFB',ratio],na.rm=T)
  MeanEUB = mean(temp[population=='EUB',ratio],na.rm=T)
  temp[, population := factor(population,levels=c("AFB","EUB"))]
  if ( !all(temp$ratio == 1, na.rm=T)){#We test if there is several isomirs for a miRNA
      res = summary(lm(temp, formula = ratio_transformed ~ population))
      estimate = res$coefficients["populationEUB", "Estimate"]
      pvalue = res$coefficients["populationEUB", "Pr(>|t|)"]
    return(data.table(isomir = isomir_ID_to_test, miRNA = temp$mirID[1], condition_test = condition_to_test, MeanAFB, MeanEUB, log2FC_EminusA = estimate, pvalue ))
  }else{
    return(data.table(isomir = character(0), miRNA = character(0), condition_test = numeric(0),  MeanAFB = numeric(0), MeanEUB = numeric(0), log2FC_EminusA = numeric(0), pvalue = numeric(0)))
  }
}



# Run for every condition
results = list()
for(c in condIndex){
  j=0
  for (isomir_to_test in isomir_ratio_melted[, unique(isomir_ID)]){
    j=j+1
    print(paste(j,isomir_to_test, c, 'pop'))
    results[[paste(isomir_to_test, c)]] = make_modelisation(isomir_to_test, c)
  }
}

# combine results and write
results = rbindlist(results)
results[, fdr := p.adjust(pvalue, method = "fdr")]
results[, pbonf := p.adjust(pvalue)]
results = results[order(pvalue)]

write.table(results, sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/isomirs_differentially_expressed_between_populations.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)


#####################################################
##    test differential response between pops      ##
#####################################################

# function 
make_modelisation_response <-function(isomir_ID_to_test, condition_to_test){
  temp = isomir_response_melted[(isomir_ID == isomir_ID_to_test) & (condition %in% c(condition_to_test))]
  MeanAFB =mean(temp[population=='AFB',response], na.rm=T)
  MeanEUB =mean(temp[population=='EUB',response], na.rm=T)
  temp[, population := factor(population,levels=c("AFB","EUB"))] 
  if ( !all(temp$ratio == 1, na.rm=T)){#We test if the miRNA is truly expressed in that condition
    res = summary(lm(temp, formula = response_transformed ~ population))
    estimate = res$coefficients["populationEUB", "Estimate"]
    pvalue = res$coefficients["populationEUB", "Pr(>|t|)"]
    return(data.table(isomir = isomir_ID_to_test, miRNA = temp$mirID[1], condition_test = condition_to_test, MeanAFB, MeanEUB, log2FC_EminusA = estimate, pvalue))
  }else{
    return(data.table(isomir = character(0), miRNA = character(0), condition_test = numeric(0), MeanAFB = numeric(0), MeanEUB = numeric(0), log2FC_EminusA = numeric(0), pvalue = numeric(0)))
  }
}

# Run for every condition
results = list()
for(c in condIndex[2:5]){
  j=0
  for (isomir_to_test in isomir_ratio_melted[, unique(isomir_ID)]){
    j=j+1
    print(paste(j,isomir_to_test, c,'pop response'))
    results[[paste(isomir_to_test, c)]] = make_modelisation_response(isomir_to_test, c)
  }
}

# combine results and write
results = rbindlist(results)
results[, fdr := p.adjust(pvalue, method = "fdr")]
results[, pbonf := p.adjust(pvalue)]
results = results[order(pvalue)]

write.table(results, sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/isomirs_differential_response_between_populations.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)



########################################################################################################
###            Bayesian Modelling of isomiR ratios to identify sharing across conditions             ###
########################################################################################################

				
allModels_STIM=cbind(NS=rep(0:1,16),
				LPS=rep(rep(0:1,e=16),1),
				PAM3CSK4=rep(rep(0:1,e=8),2),
				R848=rep(rep(0:1,e=4),4),
				IAV=rep(rep(0:1,e=2),8))
rownames(allModels_STIM)=apply(allModels_STIM,1,paste,collapse='')

				

make_modelisation_joint_isomir <-function(model_to_test,isomir_ID_to_test){
  temp = isomir_ratio_melted[(isomir_ID == isomir_ID_to_test)]
  temp[, population := factor(population,levels=c("AFB","EUB"))]
  temp[, condition := factor(condition,levels=condIndex)]
  temp[, condition_population := (population=='EUB')*model_to_test[as.numeric(condition)]]

  if (any(temp$ratio != 1,na.rm=T)){#We test if there are several isomirs for a miRNA
	# res = summary(glm(temp, formula = ratio_transformed ~ condition + population, family = "gaussian"))
    mod = lm(temp, formula = ratio_transformed ~ condition + condition_population)
    Likelihood = logLik(mod)
    if( sum(model_to_test)>0 ){   
        Coeffs=summary(mod)$coeff
        estimate=Coeffs["condition_population","Estimate"]
        pvalue=Coeffs["condition_population","Pr(>|t|)"]  
    }else{
        estimate=NA
        pvalue=NA
        }
    return(data.table(isomir = isomir_ID_to_test, miRNA = temp$mirID[1], condition_test = paste(model_to_test,collapse=''), logLik = as.numeric(Likelihood), df = attr(Likelihood,'df'), log2FC_EminusA = estimate, pval = pvalue))
  }else{
    return(data.table(isomir = character(0), miRNA = character(0), condition_test = numeric(0),  logLik  = numeric(0), df = numeric(0), log2FC_EminusA = numeric(0), pval = numeric(0)))
  }
}



results = list()
j=0
for (isomir_to_test in isomir_ratio_melted[, unique(isomir_ID)]){
    j=j+1
    print(paste(j,isomir_to_test,'pop, joint'))
    results[[paste(isomir_to_test)]] = rbindlist(apply(allModels_STIM,1,make_modelisation_joint_isomir,isomir_to_test))
  }
results = rbindlist(results)
results[,ProbModel:=exp(logLik-min(logLik))/sum(exp(logLik-min(logLik))),by=.(miRNA,isomir)]
write.table(results, sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/isomirs_differentially_expressed_between_populations_Likelihoods.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)
