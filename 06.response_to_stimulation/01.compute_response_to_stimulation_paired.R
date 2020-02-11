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

library(lme4)
library(lmerTest)

print("loading data")
corrected_count = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""), sep="\t")

mir_count_melted = melt(corrected_count, id = c("ID"),
                               variable.name = "sample",
                               value.name = "count")
mir_count_melted[, individual := substr(sample, 1,6)]
mir_count_melted[, condition := substr(sample, 8,8)]
mir_count_melted[, population := substr(sample, 1,3)]
inverseNormalRankTransform=function(x){n=length(x);
									qnorm(frank(x,na.last=FALSE,ties.method='random')/(n+1),mean(x,na.rm=T),sd(x,na.rm=T))}
mir_count_melted[, count_transformed := inverseNormalRankTransform(count),by=.(ID)]


make_modelisation <-function(mir_ID_to_test, condition_to_test){
  temp = mir_count_melted[(ID == mir_ID_to_test) & (condition %in% c(1,condition_to_test))]
  baseMean =mean(temp[condition==1,count])
  temp[, condition := factor(condIndex[as.numeric(condition)],levels=condIndex)] #to be sure it doesn't use the different values
  if ( !all(temp$count == 0)){#We test if there is several isomirs for a miRNA
	# res = summary(glm(temp, formula = ratio_transformed ~ condition + population, family = "gaussian"))
    res = summary(lmer(temp, formula = count_transformed ~ condition + (1|individual)))
    estimate = res$coefficients[paste("condition", condIndex[condition_to_test], sep=""), "Estimate"]
    pvalue = res$coefficients[paste("condition", condIndex[condition_to_test], sep=""), "Pr(>|t|)"]
    return(data.table(miRNA = mir_ID_to_test, condition_test = condIndex[condition_to_test], log2FC = estimate, pvalue , baseMean))
  }else{
    return(data.table(miRNA = character(0), condition_test = numeric(0), log2FC = numeric(0), pvalue = numeric(0), baseMean=numeric(0)))
  }
}

results = list()
for(c in 2:5){
  j=0
  for (mir_to_test in mir_count_melted[, unique(ID)]){
    j=j+1
    print(paste(j,mir_to_test, c))
    results[[paste(mir_to_test, c)]] = make_modelisation(mir_to_test, c)
  }
}

results = rbindlist(results)
results[, fdr := p.adjust(pvalue, method = "fdr")]
results[, pbonf := p.adjust(pvalue)]
results = results[order(pvalue)]

write.table(results, sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/miRNA_differentially_expressed_in_conditions.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)

################################################################
####         create joint model of miRNA expression         ####
################################################################

allModels_STIM=cbind(NS=rep(0,8),
				LPS=rep(rep(0:1,e=8),1),
				PAM3CSK4=rep(rep(0:1,e=4),2),
				R848=rep(rep(0:1,e=2),4),
				IAV=rep(rep(0:1,e=1),8))
				

make_modelisation_joint_mir <-function(model_to_test,mir_ID_to_test){
  temp = mir_count_melted[(ID == mir_ID_to_test)]
  temp[, condition := model_to_test[as.numeric(condition)]]
    res = logLik(lmer(temp, formula = count_transformed ~ condition + (1|individual) ))
    return(data.table(miRNA = mir_ID_to_test, condition_test = paste(model_to_test,collapse=''), logLik = as.numeric(res), df = attr(res,'df')))
  }

results = list()
j=0
for (mir_ID_to_test in mir_count_melted[, unique(ID)]){
    j=j+1
    print(paste(j,mir_ID_to_test))
    results[[paste(mir_ID_to_test)]] = rbindlist(apply(allModels_STIM,1,make_modelisation_joint_mir,mir_ID_to_test))
  }
results = rbindlist(results)
results[,ProbModel:=exp(logLik-min(logLik))/sum(exp(logLik-min(logLik))),by=.(miRNA)]
write.table(results, sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/miRNA_differentially_expressed_in_conditions_Likelihoods.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)

# for (population in c("None", "AFB", "EUB")){
#     for (cond in 2:5){
#     
#     print(paste(cond, population))
#     ################
#     ##Loading data##
#     ################
#     print("loading data")
#     corrected_count = read.delim(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""), sep="\t")
#     colnames(corrected_count) = gsub("\\.", "-", colnames(corrected_count))
# 
#     ##################
#     ##Filtering data##
#     ##################
#     print("filtering")
#     for_filtering = data.table(samples = colnames(corrected_count))
#     for_filtering[, pop := substr(samples, 1, 3)]
#     for_filtering[, condition := substr(samples, 8,8)]
#     if (population!="None"){
#       for_filtering = for_filtering[pop == population]
#     }
#     for_filtering = for_filtering[(condition == 1) | (condition == cond)]
# 
#     corrected_count = corrected_count[, for_filtering$samples]
# 
# 
#     print("creating metadata")
#     colData = data.table(samples = colnames(corrected_count))
#     colData[, pop := factor(substr(samples, 1, 3))]
#     colData[, condition := factor(substr(samples, 8,8), levels = c(1,cond))]
#     colData[, individual :=  factor(substr(samples, 1,6))]
#     rownames(colData) = colData$samples
# 
#     print("Launching DESEQ")
#     design <- ~ individual + condition
#     dds = DESeqDataSetFromMatrix(countData = round(corrected_count), colData = colData, design = design)
#     sizeFactors(dds)<-1
# 
#     dds = estimateDispersions(dds) ##Carefull, it's pretty long
#     dds = nbinomWaldTest(dds)
#     results = as.data.table(results(dds))
#     results[, miRNA := rownames(corrected_count)]
#     write.table(results, paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/data/06.response_to_stimulation/paired/response_to_stimulation_basalVScond", cond,ifelse(population=="None", "", paste("_",population, sep = "")),".tsv",sep=""), sep="\t", quote = F, row.names = F)
#   }
# }
# 
# 
# results=list()
# for(i in 2:5){
#     results[[i-1]]=fread(sprintf('%s/Maxime/miRNA_V2/data/06.response_to_stimulation/paired/response_to_stimulation_basalVScond%s.tsv',EVO_IMMUNO_POP,i))
#     results[[i-1]]$cond=condIndex[i]
# }
# results=rbindlist(results)
# fwrite(results,sprintf('%s/Maxime/miRNA_V2/data/06.response_to_stimulation/paired/response_to_stimulation_all.tsv',EVO_IMMUNO_POP),sep='\t')

# 
# print("loading data")
# corrected_count = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""), sep="\t")
# 
# 
# mir_count_melted = melt(corrected_count, id = c("ID"),
#                                variable.name = "sample",
#                                value.name = "count")
# mir_count_melted[, individual := substr(sample, 1,6)]
# mir_count_melted[, condition := substr(sample, 8,8)]
# mir_count_melted[, population := substr(sample, 1,3)]
# inverseNormalRankTransform=function(x){n=length(x);
# 									qnorm(frank(x,na.last=FALSE,ties.method='random')/(n+1),mean(x,na.rm=T),sd(x,na.rm=T))}
# mir_count_melted[, count_transformed := inverseNormalRankTransform(count),by=.(ID)]
