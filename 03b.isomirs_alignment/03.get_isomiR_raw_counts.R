################################################################
##This script aims to extract the raw miRNA transcripts counts##
################################################################

###########
##Imports##
###########
suppressMessages(require("reshape"))
suppressMessages(require("DESeq2"))
suppressMessages(require(ggplot2))
suppressMessages(require("sva"))

###############################
##Functions of multiple usage##
###############################
plot_pca <-function(data_table, samples_to_get){
  pca<-prcomp(t(data_table[, mget(samples_to_get)]),center=T,scale=F)
  pcs.prop.var<-round(summary(pca)$importance[2,]*100,1)
  to_plot = as.data.table(pca$x[,1:10])
  to_plot[, library_ID := samples_to_get]
  to_plot[, color := paste(substr(library_ID,8,8), substr(library_ID, 1,3), sep="_")]
  p <- ggplot(to_plot, aes(x=PC1, y = PC2, color = color))
  p <- p + geom_point()
  p <- p + scale_color_manual(values = colERC)
  p <- p + xlab(paste("PC1 : ", pcs.prop.var[1], "%", sep=""))
  p <- p + ylab(paste("PC2 : ", pcs.prop.var[2], "%", sep=""))
  p <- p + theme_bw()
  return(p)
}


remove_confounder <- function(data,variables, samples_to_get){
  temp_data = data[, mget(samples_to_get)]
  temp_data = t(temp_data)
  for(i in 1:ncol(temp_data)){
    fit = lm(temp_data[, i]~.,data=as.data.frame(variables))
    temp_data[, i] = mean(temp_data[, i])+fit$residuals
  }
  pouloup = as.data.table(t(temp_data))
  pouloup[, isomir_ID := data$isomir_ID]
  test = pouloup[1:.N]
  return(test)
}


######################################
##Getting reads in a workable format##
######################################

raw_isomiR_counts_melted = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03b.isomirs_alignment/miRNA_transcripts_raw_counts_melded.tsv", sep=""))
raw_isomiR_counts_melted[, isomir_ID := paste(chromosome , start , end , strand , sequence , subsitutions , mir1 , mir2 , hairpin, sep="_")]
#we keep the notation the most diverse possible
first <- function(x){
  if (is.na(x[1])){
    return(NA)
  }else{
    return(x[1])
  }
}
raw_isomiR_counts = dcast(raw_isomiR_counts_melted,chromosome + start + end + strand + sequence + subsitutions + mir1 + mir2 + hairpin + isomir_ID ~ library_ID, value.var = "count", first, fill=0 )


######################
##Loading covariates##
######################
all_test_var_cat = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/covariates/all_categorical_variables.tsv", sep=""))
names(all_test_var_cat)[1] = "library_ID"
all_test_var_cont<-fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/covariates/all_continuous_variables_V2.0_MR.tsv", sep=""))
names(all_test_var_cont)[1] = "library_ID"

######################
##Color for plotting##
######################
colERC=c("#969696", "#525252", "#FB9A99", "#E31A1C", "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4", "#CAB2D6", "#6A3D9A")
names(colERC)<-sort(unique(paste(all_test_var_cat$condition,all_test_var_cat$pop,sep="_")))

#####################################
##We remove rows with to few counts##
#####################################
samples_names = sort(names(raw_isomiR_counts)[c(grep("EUB", names(raw_isomiR_counts)), grep("AFB", names(raw_isomiR_counts)))])
write.table(samples_names,file=sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/sample_names_1001_libraries.tsv",EVO_IMMUNO_POP),sep='\t',quote=F,row.names=F)

mean_expression_by_isomir = raw_isomiR_counts_melted[, .(mean_expresion = sum(count)/length(samples_names)) , by = c("isomir_ID")]
isomirs_passing_expression_threshold = mean_expression_by_isomir[mean_expresion>1, isomir_ID]
raw_isomiR_counts = raw_isomiR_counts[isomir_ID %in% isomirs_passing_expression_threshold]
raw_isomiR_counts_melted = raw_isomiR_counts_melted[isomir_ID %in% isomirs_passing_expression_threshold]

###################
## Transform RPM ##
###################
isomirs_rpm_count_rounded = raw_isomiR_counts[1:.N]
for (samp in samples_names){
  print(samp)
  isomirs_rpm_count_rounded[, eval(samp) := round(1e6*get(eval(samp))/sum(get(eval(samp))))]
}

################
##Scale counts##
################
colData = all_test_var_cat[, .(pop = as.factor(pop), condition = as.factor(condition))]
dds<-DESeqDataSetFromMatrix(countData=isomirs_rpm_count_rounded[, mget(samples_names)], colData=colData, design= ~ pop + condition)
dds<-estimateSizeFactors(dds)
isomir_counts_file_norm<-as.data.table(counts(dds,normalized=TRUE))
isomir_counts_file_norm[, isomir_ID := isomirs_rpm_count_rounded$isomir_ID]

#############################
##Remove doublet of samples##
#############################
##We keep the same samples as for miRNA analysis
miRNA_expression = read.delim(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""), sep="\t")
colnames(miRNA_expression) = gsub("\\.", "-", colnames(miRNA_expression))
samples_names = colnames(miRNA_expression)
miRNA_toKeep=rownames(miRNA_expression)
rm(miRNA_expression)

all_test_var_cat = all_test_var_cat[match(samples_names, library_ID)]
all_test_var_cont = all_test_var_cont[match(samples_names, library_ID)]
isomir_counts_file_norm = isomir_counts_file_norm[isomirs_rpm_count_rounded$mir1%in%miRNA_toKeep & !grepl('N',isomirs_rpm_count_rounded$seq), mget(c("isomir_ID", samples_names))]

fwrite(isomir_counts_file_norm, file=sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/isomiR_counts.RPM.raw_%s_V2.0_MR.tsv",EVO_IMMUNO_POP,nrow(isomir_counts_file_norm)), sep="\t")
write.table(samples_names,file=sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/sample_names_977_highQuality.tsv",EVO_IMMUNO_POP),sep='\t',quote=F,row.names=F)
write.table(miRNA_toKeep,file=sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/mirID_658_expressed.tsv",EVO_IMMUNO_POP),sep='\t',quote=F,row.names=F)

