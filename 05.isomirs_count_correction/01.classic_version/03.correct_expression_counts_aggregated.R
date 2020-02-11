#######################################################################################################################
## This script aims to correct the raw miRNA transcripts after aggregation of minor isomiRs, with all substitutions  ##
#######################################################################################################################

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

#############################
## get Raw Read Count data ##
#############################
raw_isomiR_counts=fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/05.isomirs_count_correction/miRNA_isomiR_raw_counts_aggregated_cast.tsv", sep=""),sep='\t')

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
##Transform count (log2RPM)##
#############################
isomir_counts_file_norm_transformed = cbind( isomir_counts_file_norm[, .(isomir_ID)], log2(isomir_counts_file_norm[, mget(samples_names)] + 1))
p <- plot_pca(isomir_counts_file_norm_transformed, samples_names)

##########################
##Correct for GC content##
##########################
GC_prime<-all_test_var_cont[, (GC_pct_miRNAaligned-mean(GC_pct_miRNAaligned))/sd(GC_pct_miRNAaligned)]
ReadLength_prime<-all_test_var_cont[, (MeanReadLength-mean(MeanReadLength))/sd(MeanReadLength)]

isomir_counts_file_norm_transformed_GC_corrected = remove_confounder(isomir_counts_file_norm_transformed,  data.frame(GC_prime,ReadLength_prime), samples_names)
isomir_counts_file_norm_transformed_GC_corrected = isomir_counts_file_norm_transformed_GC_corrected[, mget(c("isomir_ID",samples_names) )]
p <- plot_pca(isomir_counts_file_norm_transformed_GC_corrected, samples_names)

###################################
##Remove batch effect with combat##
###################################
isomir_counts_file_norm_log_GCCorrected_df = as.data.frame(isomir_counts_file_norm_transformed_GC_corrected)
rownames(isomir_counts_file_norm_log_GCCorrected_df) = isomir_counts_file_norm_transformed_GC_corrected$isomir_ID
isomir_counts_file_norm_log_GCCorrected_df$isomir_ID = NULL

# date_manip #
batch<-all_test_var_cat$date_manip
mod<-model.matrix(~as.factor(condition) * as.factor(pop),data=all_test_var_cat)
iso_mir_combat_counts_file_gc_norm<-ComBat(dat=as.matrix(isomir_counts_file_norm_log_GCCorrected_df),batch=batch,mod=mod) # do on log2 values

#library_prep_both #
batch<-all_test_var_cat$library_prep_both
mod<-model.matrix(~as.factor(condition) * as.factor(pop),data=all_test_var_cat)
iso_mir_combat_counts_file_gc_norm_2<-ComBat(dat=as.matrix(iso_mir_combat_counts_file_gc_norm),batch=batch,mod=mod) # do on log2 values

#sequencing lane
batch<-paste(all_test_var_cat$machine, all_test_var_cat$lane)
mod<-model.matrix(~as.factor(condition) * as.factor(pop),data=all_test_var_cat)
iso_mir_combat_counts_file_gc_norm_3<-ComBat(dat=as.matrix(iso_mir_combat_counts_file_gc_norm_2),batch=batch,mod=mod) # do on log2 values

########################
## Untransform counts ##
########################
iso_NT_combat_counts_file_gc_norm_3<-(2^iso_mir_combat_counts_file_gc_norm_3)-1
iso_NT_combat_counts_file_gc_norm_3[iso_NT_combat_counts_file_gc_norm_3<0]<-0
iso_NT_combat_counts_file_gc_norm_3_dt = as.data.table(iso_NT_combat_counts_file_gc_norm_3)
iso_NT_combat_counts_file_gc_norm_3_dt[, isomir_ID := isomir_counts_file_norm_transformed_GC_corrected$isomir_ID]
iso_NT_combat_counts_file_gc_norm_3_dt[, mir_ID := isomir_counts_file_norm_transformed_GC_corrected$mir_ID]

iso_NT_combat_counts_file_gc_norm_3_logged = log2(iso_NT_combat_counts_file_gc_norm_3+1)

write.table(iso_NT_combat_counts_file_gc_norm_3, sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_counts_aggregated.RPM.GCRL_Batch_lane_corrected.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = T)
write.table(iso_NT_combat_counts_file_gc_norm_3_logged, sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_counts_aggregated.log2RPM.GCRL_Batch_lane_corrected.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = T)


########################
##   Compute ratios   ##
########################

raw_isomiR_counts[,mirID:=paste(hsa_ID,MIMAT,MI_ID,sep='_')]
all(iso_NT_combat_counts_file_gc_norm_3_dt$isomir_ID== raw_isomiR_counts$isomir_ID)
iso_NT_combat_counts_file_gc_norm_3_dt[, mirID := raw_isomiR_counts$mirID]
isomir_counts_melted=melt(iso_NT_combat_counts_file_gc_norm_3_dt, id.vars=c("isomir_ID","mirID"), measure.vars=samples_names,  variable.name = "library_ID", value.name = "count")

mir_grouped_counts_melted=isomir_counts_melted[,.(total_miR_counts=sum(count)),by=.(mirID,library_ID)]
isomir_counts_melted=merge(isomir_counts_melted,mir_grouped_counts_melted,by=c('mirID','library_ID'))
isomir_counts_melted[,ratio:=count/total_miR_counts]

isomir_ratio_cast=dcast(isomir_counts_melted,mirID+isomir_ID ~ library_ID, value.var = "ratio", first, fill=NA)

fwrite(isomir_ratio_cast,file=sprintf("%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_ratios_aggregated.GCRL_Batch_lane_corrected.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = T)


