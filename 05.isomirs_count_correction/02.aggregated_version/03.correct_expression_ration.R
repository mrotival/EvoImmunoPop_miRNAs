############################################################################
##This script aims to correct and analise the raw miRNA transcripts counts##
############################################################################

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


remove_confounder <- function(data,variable, samples_to_get){
  temp_data = data[, mget(samples_to_get)]
  temp_data = t(temp_data)
  for(i in 1:ncol(temp_data)){
    fit = lm(temp_data[, i]~variable)
    temp_data[, i] = temp_data[, i]-(fit$coefficients[2]*variable)
  }
  pouloup = as.data.table(t(temp_data))
  pouloup[, isomir_ID := data$isomir_ID]
  test = pouloup[1:.N]
  return(test)
}

raw_isomiR_counts=fread(file=paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/05.isomirs_count_correction/miRNA_isomiR_raw_counts_aggregated_nosubs_cast.tsv", sep=""),sep='\t')

######################
##Loading covariates##
######################
samples_names=unlist(fread(sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/sample_names_977_highQuality.tsv",EVO_IMMUNO_POP)))
all_test_var_cat = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/covariates/all_categorical_variables.tsv", sep=""))
names(all_test_var_cat)[1] = "library_ID"
all_test_var_cat = all_test_var_cat[match(samples_names, library_ID)]

all_test_var_cont<-fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/covariates/all_continuous_variables_V2.0_MR.tsv", sep=""))
names(all_test_var_cont)[1] = "library_ID"
all_test_var_cont = all_test_var_cont[match(samples_names, library_ID)]

######################
##Color for plotting##
######################
colERC=c("#969696", "#525252", "#FB9A99", "#E31A1C", "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4", "#CAB2D6", "#6A3D9A")
names(colERC)<-sort(unique(paste(all_test_var_cat$condition,all_test_var_cat$pop,sep="_")))

#####################################
##We remove rows with to few counts##
#####################################
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
miRNA_expression = read.delim(paste(EVO_IMMUNO_POP, "Martin/miRNA/03.alignment_and_count_correction/data/miRNA_counts.log2RPM.GC_Batch_corrected.tsv", sep=""), sep="\t")
colnames(miRNA_expression) = gsub("\\.", "-", colnames(miRNA_expression))
samples_names = colnames(miRNA_expression)
isomir_counts_file_norm = isomir_counts_file_norm[, mget(c("isomir_ID", samples_names))]
all_test_var_cat = all_test_var_cat[match(samples_names, library_ID)]
all_test_var_cont = all_test_var_cont[match(samples_names, library_ID)]
rm(miRNA_expression)

################################################################################################
##Filter lowly expressed isomirs (so that we can compute GC correlation, also confuses COMBAT)##
################################################################################################
mean_count = isomir_counts_file_norm[, rowMeans(.SD[, mget(samples_names)])]
isomir_counts_file_norm = isomir_counts_file_norm[mean_count > 1]

#############################
##Transform count (log2RPM)##
#############################
isomir_counts_file_norm_transformed = cbind( isomir_counts_file_norm[, .(isomir_ID)], log2(isomir_counts_file_norm[, mget(samples_names)] + 1))
p <- plot_pca(isomir_counts_file_norm_transformed, samples_names)


##########################
##Correct for GC content##
##########################
GC_prime<-all_test_var_cont[, (GC-mean(GC))/sd(GC)]
isomir_counts_file_norm_transformed_GC_corrected = remove_confounder(isomir_counts_file_norm_transformed,  GC_prime, samples_names)
isomir_counts_file_norm_transformed_GC_corrected = isomir_counts_file_norm_transformed_GC_corrected[, mget(c("isomir_ID",samples_names) )]
p <- plot_pca(isomir_counts_file_norm_transformed_GC_corrected, samples_names)


###################################
##Remove batch effect with combat##
###################################
isomir_counts_file_norm_log_GCCorrected_df = as.data.frame(isomir_counts_file_norm_transformed_GC_corrected)
rownames(isomir_counts_file_norm_log_GCCorrected_df) = isomir_counts_file_norm_transformed_GC_corrected$isomir_ID
isomir_counts_file_norm_log_GCCorrected_df$isomir_ID = NULL

#library_prep_both #
batch<-all_test_var_cat$library_prep_both
mod<-model.matrix(~as.factor(condition) * as.factor(pop),data=all_test_var_cat)
iso_mir_combat_counts_file_gc_norm<-ComBat(dat=isomir_counts_file_norm_log_GCCorrected_df,batch=batch,mod=mod) # do on log2 values

# date_manip #
batch<-all_test_var_cat$date_manip
mod<-model.matrix(~as.factor(condition) * as.factor(pop),data=all_test_var_cat)
iso_mir_combat_counts_file_gc_norm_2<-ComBat(dat=iso_mir_combat_counts_file_gc_norm,batch=batch,mod=mod) # do on log2 values


########################
## Untransform counts ##
########################
iso_NT_combat_counts_file_gc_norm_2<-(2^iso_mir_combat_counts_file_gc_norm_2)-1
iso_NT_combat_counts_file_gc_norm_2[iso_NT_combat_counts_file_gc_norm_2<0]<-0
iso_NT_combat_counts_file_gc_norm_2_dt = as.data.table(iso_NT_combat_counts_file_gc_norm_2)
iso_NT_combat_counts_file_gc_norm_2_dt[, isomir_ID := isomir_counts_file_norm_transformed_GC_corrected$isomir_ID]


iso_NT_combat_counts_file_gc_norm_2_logged = log2(iso_NT_combat_counts_file_gc_norm_2+1)
iso_NT_combat_counts_file_gc_norm_2 = cbind(isomirs_informations[match(rownames(iso_NT_combat_counts_file_gc_norm_2), isomir_ID)], iso_NT_combat_counts_file_gc_norm_2)
iso_NT_combat_counts_file_gc_norm_2_logged = cbind(isomirs_informations[match(rownames(iso_NT_combat_counts_file_gc_norm_2_logged), isomir_ID)], iso_NT_combat_counts_file_gc_norm_2_logged)

write.table(iso_NT_combat_counts_file_gc_norm_2, "../../data/isomiR_counts.RPM.GC_Batch_corrected.aggregated.tsv", sep="\t", quote = F, row.names = F)
write.table(iso_NT_combat_counts_file_gc_norm_2_logged, "../../data/miRNA_counts.log2RPM.GC_Batch_corrected.aggregated.tsv", sep="\t", quote = F, row.names = F)
