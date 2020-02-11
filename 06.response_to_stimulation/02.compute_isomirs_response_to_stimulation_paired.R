##################################################################################
##Aim of the script : investigate which isomirs change the *ratio* of expression##
##Upon stimulation                                                              ##
##################################################################################

##############
##Librairies##
##############


#############
##Load data##
#############

isomir_ratio_cast=fread(sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_ratios_aggregated_nosubs.GCRL_Batch_lane_corrected.tsv",EVO_IMMUNO_POP))
isomir_ratio_cast=isomir_ratio_cast[-grep('other',isomir_ID)]
dim(isomir_ratio_cast)
isomir_ratio_melted = melt(isomir_ratio_cast, id = c("mirID", "isomir_ID"),
                               variable.name = "sample",
                               value.name = "ratio")
isomir_ratio_melted[, individual := substr(sample, 1,6)]
isomir_ratio_melted[, condition := substr(sample, 8,8)]
isomir_ratio_melted[, population := substr(sample, 1,3)]
inverseNormalRankTransform=function(x){n=length(x);
									qnorm(frank(x,na.last=FALSE,ties.method='random')/(n+1),mean(x,na.rm=T),sd(x,na.rm=T))}
isomir_ratio_melted[, ratio_transformed := inverseNormalRankTransform(ratio),by=.(mirID,isomir_ID)]

# isomir_annot=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_FULL_V2.0.tsv',EVO_IMMUNO_POP))
# isomir_ratio_melted



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

make_modelisation <-function(isomir_ID_to_test, condition_to_test){
  temp = isomir_ratio_melted[(isomir_ID == isomir_ID_to_test) & (condition %in% c(1,condition_to_test))]
  temp[, condition := factor(condIndex[as.numeric(condition)],levels=condIndex)] #to be sure it doesn't use the different values
  if (FALSE %in% (temp$ratio == 1)){#We test if there is several isomirs for a miRNA
	# res = summary(glm(temp, formula = ratio_transformed ~ condition + population, family = "gaussian"))
    res = summary(lmer(temp, formula = ratio_transformed ~ condition + (1|individual)))
    estimate = res$coefficients[paste("condition", condIndex[condition_to_test], sep=""), "Estimate"]
    pvalue = res$coefficients[paste("condition", condIndex[condition_to_test], sep=""), "Pr(>|t|)"]
    return(data.table(isomir = isomir_ID_to_test, miRNA = temp$mirID[1], condition_test = condIndex[condition_to_test], beta = estimate, pvalue ))
  }else{
    return(data.table(isomir = character(0), miRNA = character(0), condition_test = numeric(0), beta = numeric(0), pvalue = numeric(0)))
  }
}

results = list()
for(c in 2:5){
  j=0
  for (isomir_to_test in isomir_ratio_melted[, unique(isomir_ID)]){
    j=j+1
    print(paste(j,isomir_to_test, c))
    results[[paste(isomir_to_test, c)]] = make_modelisation(isomir_to_test, c)
  }
}

results = rbindlist(results)
results[, fdr := p.adjust(pvalue, method = "fdr")]
results[, pbonf := p.adjust(pvalue)]
results = results[order(pvalue)]

write.table(results, sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/isomirs_differentially_expressed_in_conditions.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)


########################################################################################################
###            Bayesian Modelling of isomiR ratios to identify sharing across conditions             ###
########################################################################################################

				
allModels_STIM=cbind(NS=rep(0,8),
				LPS=rep(rep(0:1,e=8),1),
				PAM3CSK4=rep(rep(0:1,e=4),2),
				R848=rep(rep(0:1,e=2),4),
				IAV=rep(rep(0:1,e=1),8))
				

make_modelisation_joint <-function(model_to_test,isomir_ID_to_test){
  temp = isomir_ratio_melted[(isomir_ID == isomir_ID_to_test)]
  temp[, condition := model_to_test[as.numeric(condition)]]

  if (any(temp$ratio != 1,na.rm=T)){#We test if there is several isomirs for a miRNA
	# res = summary(glm(temp, formula = ratio_transformed ~ condition + population, family = "gaussian"))
    res = logLik(lmer(temp, formula = ratio_transformed ~ condition + (1|individual)))
    return(data.table(isomir = isomir_ID_to_test, miRNA = temp$mirID[1], condition_test = paste(model_to_test,collapse=''), logLik = as.numeric(res), df = attr(res,'df')))
  }else{
    return(data.table(isomir = character(0), miRNA = character(0), condition_test = numeric(0),  logLik  = numeric(0), df = numeric(0)))
  }
}

results = list()
j=0
for (isomir_to_test in isomir_ratio_melted[, unique(isomir_ID)]){
    j=j+1
    print(paste(j,isomir_to_test))
    results[[paste(isomir_to_test)]] = rbindlist(apply(allModels_STIM,1,make_modelisation_joint,isomir_to_test))
  }
results = rbindlist(results)
results[,ProbModel:=exp(logLik-min(logLik))/sum(exp(logLik-min(logLik))),by=.(miRNA,isomir)]
write.table(results, sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/isomirs_differentially_expressed_in_conditions_Likelihoods.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)


test_model_joint <-function(model_to_test,isomir_ID_to_test){
    model_to_test=as.numeric(model_to_test)
    modelID=as.numeric(apply(allModels_STIM,1,paste,collapse=''))
    modelVect=allModels_STIM[match(model_to_test,modelID),]
  temp = isomir_ratio_melted[(isomir_ID == isomir_ID_to_test)]
  temp[, condition := modelVect[as.numeric(condition)]]

  if (any(temp$ratio != 1,na.rm=T)){#We test if there is several isomirs for a miRNA
	# res = summary(glm(temp, formula = ratio_transformed ~ condition + population, family = "gaussian"))
    res = summary(lmer(temp, formula = ratio_transformed ~ condition + (1|individual)))$coeff[2,c("Estimate","Pr(>|t|)")]
    return(data.table(isomir = isomir_ID_to_test, miRNA = temp$mirID[1], condition_test = paste(model_to_test,collapse=''), beta = as.numeric(res[1]), pval =as.numeric(res[2])))
  }else{
    return(data.table(isomir = character(0), miRNA = character(0), condition_test = numeric(0),  beta  = numeric(0), pval=0))
  }
}

# input = fread(sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/isomirs_differentially_expressed_in_conditions_Likelihoods.tsv",EVO_IMMUNO_POP), colClasses=c("character","character","character","numeric","integer","numeric"))
results=list()
for (i in which(input$condition_test>0)){
    cat(i,'')
    results[[i]]=test_model_joint(input$condition_test[i],input$isomir[i])
}
results = rbindlist(results)
results_2 = merge(results[,mget(c('isomir','miRNA','condition_test','beta','pval'))],input,by=c('isomir','miRNA','condition_test'))
results_2[,ProbModel:=exp(logLik-min(logLik))/sum(exp(logLik-min(logLik))),by=.(miRNA,isomir)]
write.table(results_2, sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/isomirs_differentially_expressed_in_conditions_Likelihoods_withBeta.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)

##############################################################################################
##     Test for differential expression upon stimulation  (adjusted on end site shift )    ###
##############################################################################################

condIndex=c("NS","LPS","PAM3CSK4","R848","IAV")

library(lme4)
library(lmerTest)

Pct_modif_perSample=fread(file=sprintf('%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/Pct_modifs_perSample.txt',EVO_IMMUNO_POP))


make_modelisation_adj <-function(isomir_ID_to_test, condition_to_test){
  temp = isomir_ratio_melted[(isomir_ID == isomir_ID_to_test) & (condition %in% c(1,condition_to_test))]
  temp[, condition := factor(condIndex[as.numeric(condition)],levels=condIndex)] #to be sure it doesn't use the different values
  temp[,shift_3p := Pct_modif_perSample$shift_3p[match(sample,Pct_modif_perSample$ID)]]
  if (FALSE %in% (temp$ratio == 1)){#We test if there is several isomirs for a miRNA
	# res = summary(glm(temp, formula = ratio_transformed ~ condition + population, family = "gaussian"))
    res = summary(lmer(temp, formula = ratio_transformed ~ condition + (1|individual) + shift_3p))
    estimate = res$coefficients[paste("condition", condIndex[condition_to_test], sep=""), "Estimate"]
    pvalue = res$coefficients[paste("condition", condIndex[condition_to_test], sep=""), "Pr(>|t|)"]
    return(data.table(isomir = isomir_ID_to_test, miRNA = temp$mirID[1], condition_test = condIndex[condition_to_test], beta = estimate, pvalue ))
  }else{
    return(data.table(isomir = character(0), miRNA = character(0), condition_test = numeric(0), beta = numeric(0), pvalue = numeric(0)))
  }
}

results = list()
for(c in 2:5){
  j=0
  for (isomir_to_test in isomir_ratio_melted[, unique(isomir_ID)]){
    j=j+1
    print(paste(j,isomir_to_test, c))
    results[[paste(isomir_to_test, c)]] = make_modelisation_adj(isomir_to_test, c)
  }
}

results = rbindlist(results)
results[, fdr := p.adjust(pvalue, method = "fdr")]
results[, pbonf := p.adjust(pvalue)]
results = results[order(pvalue)]

write.table(results, sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/isomirs_differentially_expressed_in_conditions_3pshift.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)




#######################################################################################################################
###     Bayesian Modelling of isomiR ratios to identify sharing across conditions  (adjusted on end site shift )    ###
#######################################################################################################################

allModels_STIM=cbind(NS=rep(0,8),
				LPS=rep(rep(0:1,e=8),1),
				PAM3CSK4=rep(rep(0:1,e=4),2),
				R848=rep(rep(0:1,e=2),4),
				IAV=rep(rep(0:1,e=1),8))
				
Pct_modif_perSample=fread(file=sprintf('%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/Pct_modifs_perSample.txt',EVO_IMMUNO_POP))

make_modelisation_joint_adj3pShift <-function(model_to_test,isomir_ID_to_test){
  temp = isomir_ratio_melted[(isomir_ID == isomir_ID_to_test)]
  temp[, condition := model_to_test[as.numeric(condition)]]
  temp[,shift_3p := Pct_modif_perSample$shift_3p[match(sample,Pct_modif_perSample$ID)]]
  if (any(temp$ratio != 1,na.rm=T)){#We test if there is several isomirs for a miRNA
	# res = summary(glm(temp, formula = ratio_transformed ~ condition + population, family = "gaussian"))
    res = logLik(lmer(temp, formula = ratio_transformed ~ condition + (1|individual) + shift_3p))
    return(data.table(isomir = isomir_ID_to_test, miRNA = temp$mirID[1], condition_test = paste(model_to_test,collapse=''), logLik = as.numeric(res), df = attr(res,'df')))
  }else{
    return(data.table(isomir = character(0), miRNA = character(0), condition_test = numeric(0),  logLik  = numeric(0), df = numeric(0)))
  }
}

results = list()
j=0
for (isomir_to_test in isomir_ratio_melted[, unique(isomir_ID)]){
    j=j+1
    print(paste(j,isomir_to_test))
    results[[paste(isomir_to_test)]] = rbindlist(apply(allModels_STIM,1,make_modelisation_joint_adj3pShift,isomir_to_test))
  }
results = rbindlist(results)
results[,ProbModel:=exp(logLik-min(logLik))/sum(exp(logLik-min(logLik))),by=.(miRNA,isomir)]
write.table(results, sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/isomirs_differentially_expressed_in_conditions_Likelihoods_3pshift.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)





