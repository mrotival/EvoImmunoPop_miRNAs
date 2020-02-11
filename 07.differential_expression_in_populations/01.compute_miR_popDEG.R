##############
##librairies##
##############
suppressMessages(require(DESeq2))

# load GO annotation data
source(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/scripts/00a.GOannotation/01.Load_and_format_databases.R", sep=""))
GO_informations = load_BHF_UCL_miRNA()

######################
##Handling arguments##
######################

##################################################################
##     Test for differential expression upon stimulation        ##
##################################################################

condIndex=c("NS","LPS","PAM3CSK4","R848","IAV")
luq = function(x){length(unique(x))}

library(lme4)
library(lmerTest)


######################################
##     Load read count data         ##
######################################

print("loading data")
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
######################################
##    Convert into response data    ##
######################################

mir_response = dcast(mir_count_melted,ID + individual + population ~ condition, value.var=list('count_transformed','count'))
mir_response_melted=melt(mir_response,id.vars=c('ID','individual','population','count_transformed_NS','count_NS'), measure.var=list(c('count_transformed_LPS','count_transformed_PAM3CSK4','count_transformed_R848','count_transformed_IAV'),c('count_LPS','count_PAM3CSK4','count_R848','count_IAV')),value.name=c('count_transformed','count'),variable.name='stimulus')
mir_response_melted[,condition := condIndex[as.numeric(stimulus)+1]]
mir_response_melted[,response := count-count_NS]
mir_response_melted[,response_transformed := count_transformed-count_transformed_NS]


#####################################################
##    test differential expression between pops    ##
#####################################################

# function 
make_modelisation <-function(mir_ID_to_test, condition_to_test){
  temp = mir_count_melted[(ID == mir_ID_to_test) & (condition %in% c(condition_to_test))]
  MeanAFB =mean(temp[population=='AFB',count])
  MeanEUB =mean(temp[population=='EUB',count])
  temp[, population := factor(population,levels=c("AFB","EUB"))]

  if ( !all(temp$count == 0)){#We test if the miRNA is truly expressed in that condition
    res = summary(lm(temp, formula = count_transformed ~ population))
    estimate = res$coefficients["populationEUB", "Estimate"]
    pvalue = res$coefficients["populationEUB", "Pr(>|t|)"]
    return(data.table(miRNA = mir_ID_to_test, condition_test = condition_to_test, MeanAFB, MeanEUB, log2FC_EminusA = estimate, pvalue))
  }else{
    return(data.table(miRNA = character(0), condition_test = numeric(0), MeanAFB = numeric(0), MeanEUB = numeric(0), log2FC_EminusA = numeric(0), pvalue = numeric(0)))
  }
}

# Run for every condition
results = list()
for(c in condIndex){
  j=0
  for (mir_to_test in mir_response_melted[, unique(ID)]){
    j=j+1
    print(paste(j,mir_to_test, c, 'pop'))
    results[[paste(mir_to_test, c)]] = make_modelisation(mir_to_test, c)
  }
}

# combine results and write
results = rbindlist(results)
results[, fdr := p.adjust(pvalue, method = "fdr")]
results[, pbonf := p.adjust(pvalue)]
results = results[order(pvalue)]

write.table(results, sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/miRNA_differentially_expressed_between_populations.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)


#####################################################
##    test differential response between pops      ##
#####################################################

# function 
make_modelisation_response <-function(mir_ID_to_test, condition_to_test){
  temp = mir_response_melted[(ID == mir_ID_to_test) & (condition %in% c(condition_to_test))]
  MeanAFB =mean(temp[population=='AFB',response],na.rm=T)
  MeanEUB =mean(temp[population=='EUB',response],na.rm=T)
  temp[, population := factor(population,levels=c("AFB","EUB"))] 
  if ( !all(temp$count == 0)){#We test if the miRNA is truly expressed in that condition
    res = summary(lm(temp, formula = response_transformed ~ population))
    estimate = res$coefficients["populationEUB", "Estimate"]
    pvalue = res$coefficients["populationEUB", "Pr(>|t|)"]
    return(data.table(miRNA = mir_ID_to_test, condition_test = condition_to_test, MeanAFB, MeanEUB, log2FC_EminusA = estimate, pvalue))
  }else{
    return(data.table(miRNA = character(0), condition_test = numeric(0), MeanAFB = numeric(0), MeanEUB = numeric(0), log2FC_EminusA = numeric(0), pvalue = numeric(0)))
  }
}

# Run for every condition
results = list()
for(c in condIndex[2:5]){
  j=0
  for (mir_to_test in mir_count_melted[, unique(ID)]){
    j=j+1
    print(paste(j,mir_to_test, c,'pop response'))
    results[[paste(mir_to_test, c)]] = make_modelisation_response(mir_to_test, c)
  }
}

# combine results and write
results = rbindlist(results)
results[, fdr := p.adjust(pvalue, method = "fdr")]
results[, pbonf := p.adjust(pvalue)]
results = results[order(pvalue)]

write.table(results, sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/miRNA_differential_response_between_populations.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)

################################################################
####         create joint model of miRNA expression         ####
################################################################

allModels_STIM=cbind(NS=rep(0:1,16),
				LPS=rep(rep(0:1,e=16),1),
				PAM3CSK4=rep(rep(0:1,e=8),2),
				R848=rep(rep(0:1,e=4),4),
				IAV=rep(rep(0:1,e=2),8))

make_modelisation_joint_mir <-function(model_to_test,mir_ID_to_test){
  temp = mir_count_melted[(ID == mir_ID_to_test)]
  temp[, population := factor(population,levels=c("AFB","EUB"))]
  temp[, condition := factor(condition,levels=condIndex)]
  temp[, condition_population := (population=='EUB')*model_to_test[as.numeric(condition)]]
  if ( !all(temp$count == 0)){#We test if the miRNA is truly expressed in that condition
    mod = lm(temp, formula = count_transformed ~ condition + condition_population)
    Likelihood = logLik(mod)
    if( sum(model_to_test)>0 ){
        Coeffs=summary(mod)$coeff
        estimate=Coeffs["condition_population","Estimate"]
        pvalue=Coeffs["condition_population","Pr(>|t|)"]
    }else{
        estimate=NA
        pvalue=NA
        }
    return(data.table(miRNA = mir_ID_to_test, condition_test = paste(model_to_test,collapse=''), logLik = as.numeric(Likelihood), df = attr(Likelihood,'df'), log2FC_EminusA = estimate, pval = pvalue))
    }else{
    return(data.table(miRNA = character(0), condition_test = character(0), logLik = numeric(0), df = numeric(0), log2FC_EminusA=numeric(0), pval=numeric(0)))
   }
}

results = list()
j=0
for (mir_ID_to_test in mir_count_melted[, unique(ID)]){
    j=j+1
    print(paste(j,mir_ID_to_test,'pop, joint'))
    results[[paste(mir_ID_to_test)]] = rbindlist(apply(allModels_STIM,1,make_modelisation_joint_mir,mir_ID_to_test))
  }
results = rbindlist(results)
results[,ProbModel:=exp(logLik-min(logLik))/sum(exp(logLik-min(logLik))),by=.(miRNA)]
write.table(results, sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/miRNA_differentially_expressed_between_populations_Likelihoods.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)
