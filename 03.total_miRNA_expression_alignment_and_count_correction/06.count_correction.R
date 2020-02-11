################################################################
##This script aims to correct and analise the raw miRNA counts##
################################################################

###########
##Imports##
###########
suppressMessages(require("reshape"))
suppressMessages(require("DESeq2"))
suppressMessages(require(ggplot2))
suppressMessages(require("sva"))
suppressMessages(require("pcaMethods"))

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
  pouloup[, miRNA_name := data$miRNA_name]
  test = pouloup[1:.N]
  return(test)
}

######################################
##Getting reads in a workable format##
######################################
raw_miRNA_counts = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_raw_counts.tsv", sep=""))
raw_miRNA_counts = raw_miRNA_counts[, .(mimat_ID = paste(unique(miRNA_name2), collapse="//"),
                                        count = sum(raw_counts)), by = c("library_ID", "miRNA_name")]

libraries_inspected = raw_miRNA_counts[, unique(library_ID)]

raw_miRNA_counts = dcast(raw_miRNA_counts, miRNA_name ~library_ID, value.var = "count", fill=0)


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

#################
##Transform_rpm##
#################
rpm_count_rounded = raw_miRNA_counts[1:.N]
for (lib in libraries_inspected){
  rpm_count_rounded[, eval(lib) := round(1e6*get(eval(lib))/sum(get(eval(lib))))]
}

################
##Scale counts##
################
colData = all_test_var_cat[, .(pop = as.factor(pop), condition = as.factor(condition))]
dds<-DESeqDataSetFromMatrix(countData=rpm_count_rounded[, mget(libraries_inspected)], colData=colData, design= ~ pop + condition)
dds<-estimateSizeFactors(dds)
counts_file_norm<-as.data.table(counts(dds,normalized=TRUE))
counts_file_norm[, miRNA_name := rpm_count_rounded$miRNA_name]

#####################
## Detect outliers ##
#####################
# per sample correlation #
expo <- 0.01 # value retaken from Katie who calculated it using findP $maxIQR of package ops
adjusted_miRNA_count = raw_miRNA_counts[1:.N]
for (lib in libraries_inspected){
  adjusted_miRNA_count[, eval(lib) := get(eval(lib))**expo]
}
cor_matrix = cor(adjusted_miRNA_count[, mget(libraries_inspected)], adjusted_miRNA_count[, mget(libraries_inspected)])
median_corelation = apply(cor_matrix,1,median)
outlier_samples = names(median_corelation)[round(median_corelation,2)<=0.6]


########################
## Detect outliers V2 ##            DEFUNCT
########################
#median_count=apply(counts_file_norm[, mget(libraries_inspected)],1,median)
#median_cor=cor(counts_file_norm,median_count)


###############################
##Detect replicates to remove##
###############################
replicate_detection_df = data.table(samples = libraries_inspected)
replicate_detection_df[, condition := substr(samples,8,8)]
replicate_detection_df[, individual := substr(samples,1,6)]
replicate_detection_df = replicate_detection_df[, .N, by=c("condition", "individual")]
replicates_to_remove = replicate_detection_df[N>1, paste(individual, "-", condition, "_", "1", sep="")]


###########################
##Remove unwanted samples##
###########################
counts_file_norm[, eval(replicates_to_remove) := NULL]
counts_file_norm[, eval(outlier_samples):=NULL]

all_test_var_cat = all_test_var_cat[! (library_ID %in% c(replicates_to_remove, outlier_samples))]
all_test_var_cont = all_test_var_cont[! (library_ID %in% c(replicates_to_remove, outlier_samples))]

samples = all_test_var_cont[, library_ID]

##########################
##Filter low miRNA count##
##########################
mean_count<-rowMeans(counts_file_norm[, mget(samples)])
counts_file_norm<-counts_file_norm[mean_count>1,] # around 5 reads

######################
## Transform counts ##
######################
counts_file_norm_log = counts_file_norm[1:.N]
for (s in samples){
  counts_file_norm_log[, eval(s) := log2(get(eval(s)) + 1)]
}
p = plot_pca(counts_file_norm_log, samples)

############################
## Correct for GC content ##
############################
GC_prime<-all_test_var_cont[, (GC_pct_miRNAaligned-mean(GC_pct_miRNAaligned))/sd(GC_pct_miRNAaligned)]
ReadLength_prime<-all_test_var_cont[, (MeanReadLength-mean(MeanReadLength))/sd(MeanReadLength)]

counts_file_norm_log_GCCorrected = remove_confounder(counts_file_norm_log,  data.frame(GC_prime,ReadLength_prime), samples)
q = plot_pca(counts_file_norm_log_GCCorrected, samples)

######################################
## Remove batch effects with combat ##
######################################
counts_file_norm_log_GCCorrected_dt = as.data.frame(counts_file_norm_log_GCCorrected)
rownames(counts_file_norm_log_GCCorrected_dt) = counts_file_norm_log_GCCorrected_dt$miRNA_name
counts_file_norm_log_GCCorrected_dt$miRNA_name = NULL

# # date_manip #
batch<-all_test_var_cat$date_manip
mod<-model.matrix(~as.factor(condition) * as.factor(pop),data=all_test_var_cat)
combat_counts_file_gc_norm<-ComBat(dat=as.matrix(counts_file_norm_log_GCCorrected_dt),batch=batch,mod=mod) # do on log2 values

#library_prep_both #
batch<-all_test_var_cat$library_prep_both
mod<-model.matrix(~as.factor(condition) * as.factor(pop),data=all_test_var_cat)
combat_counts_file_gc_norm_2<-ComBat(dat=combat_counts_file_gc_norm,batch=batch,mod=mod) # do on log2 values

########################
## Untransform counts ##
########################
NT_combat_counts_file_gc_norm_2<-(2^combat_counts_file_gc_norm_2)-1
NT_combat_counts_file_gc_norm_2[NT_combat_counts_file_gc_norm_2<0]<-0
NT_combat_counts_file_gc_norm_2_dt = as.data.table(NT_combat_counts_file_gc_norm_2)
NT_combat_counts_file_gc_norm_2_dt[, miRNA_name := rownames(combat_counts_file_gc_norm)]

s = plot_pca(NT_combat_counts_file_gc_norm_2_dt, samples)
test = as.data.table(combat_counts_file_gc_norm_2)
d = plot_pca(test, samples)

NT_combat_counts_file_gc_norm_2_logged = log2(NT_combat_counts_file_gc_norm_2+1)
write.table(NT_combat_counts_file_gc_norm_2, paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""), sep="\t", quote = F, row.names = T)
write.table(NT_combat_counts_file_gc_norm_2_logged, paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""), sep="\t", quote = F, row.names = T)


# NT_combat_counts_file_gc_norm_2=fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""))
# NT_combat_counts_file_gc_norm_2_logged=fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""))
# samples=unlist(fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03b.isomirs_alignment/sample_names_977_highQuality.tsv", sep="")))

################################################
###### RUN PCA on Read Length adjusted data   ##
################################################
pdf(sprintf('%s/Maxime/miRNA_V2/figures/03.total_miRNA_expression_alignment_and_count_correction/PCA_ajusted_GCmiRNA_readlength.pdf',EVO_IMMUNO_POP),width=5,height=4)
test = as.data.table(NT_combat_counts_file_gc_norm_2_logged)
d = plot_pca(test, samples)
d
dev.off()

############ check the most contributing miRNAS

library(pcaMethods)
rownames(NT_combat_counts_file_gc_norm_2_logged)=NT_combat_counts_file_gc_norm_2_logged$ID
NT_combat_counts_file_gc_norm_2_logged$ID=NULL
PCs_new=pca(t(as.matrix(NT_combat_counts_file_gc_norm_2_logged)))

pop=substr(samples,1,3)
cond=substr(samples,8,8)
apply(PCs_new@scores,2,function(x){summary(lm(x~pop+cond))$coeff[2,4]})
#          PC1          PC2 
# 1.002197e-79 1.181126e-11 
apply(PCs_new@scores,2,function(x){kruskal.test(x,as.factor(pop))$p.value})
#          PC1          PC2 
# 3.229389e-13 5.154721e-08 


PCs_new@loadings[order(PCs_new@loadings[,1])[1:10],1,drop=F]
#hsa-miR-155-5p hsa-miR-146a-3p hsa-miR-125a-3p hsa-miR-146a-5p  hsa-miR-222-5p hsa-miR-146b-3p    hsa-miR-9-5p  hsa-miR-155-3p    hsa-miR-9-3p   hsa-let-7e-3p 
#     -0.2621451      -0.2446855      -0.2214307      -0.2178937      -0.2171443      -0.2089318      -0.1861198      -0.1813685      -0.1791713      -0.1454359 
PCs_new@loadings[order(-PCs_new@loadings[,1])[1:10],1,drop=F]
# hsa-miR-215-5p  hsa-miR-183-5p  hsa-miR-182-5p  hsa-miR-204-5p     hsa-miR-375 hsa-miR-196a-5p     hsa-miR-429  hsa-miR-141-5p  hsa-miR-141-3p    hsa-miR-4485 
#      0.2198262       0.1839024       0.1560916       0.1549113       0.1536521       0.1535141       0.1458437       0.1313109       0.1090166       0.1067694 
PCs_new@loadings[order(PCs_new@loadings[,2])[1:10],2,drop=F]
#hsa-miR-155-5p   hsa-miR-3196   hsa-miR-4516   hsa-miR-6087 hsa-miR-708-5p   hsa-miR-4532 hsa-miR-183-5p hsa-miR-155-3p   hsa-miR-4508   hsa-miR-4492 
#    -0.1768166     -0.1700055     -0.1603780     -0.1577698     -0.1365831     -0.1338456     -0.1258234     -0.1238152     -0.1223622     -0.1214840 
PCs_new@loadings[order(-PCs_new@loadings[,2])[1:10],2,drop=F]
#hsa-miR-181b-3p   hsa-miR-27a-5p hsa-miR-26a-2-3p   hsa-miR-33b-3p   hsa-miR-23a-5p  hsa-miR-181a-3p  hsa-miR-1304-3p  hsa-miR-4482-3p  hsa-miR-6502-5p     hsa-miR-3161 
#      0.08913468       0.07556809       0.07325892       0.06605275       0.06142554       0.05625587       0.05510889       0.05309950       0.05192095       0.04913008 


#############################################################################################
###### compare PCs of non Read Length adjusted data with read length and miRNA GC        ####
#############################################################################################

old_miRNA=fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/old/miRNA_counts.log2RPM.GC_Batch_corrected.tsv",sep=''))

PCs=pca(t(as.matrix(old_miRNA[,mget(samples)])),nPcs=6)

pdf(sprintf('%s/03_Analysis/miRNA_project/figures/correlation_oldPCs_toGC_and_readLength.pdf',HOME))
plot(as.numeric(PCs@scores[,3]),all_test_var_cont$GC_pct_miRNAaligned,pch=16,xlab='PC3',ylab='GC % miRNA',col=colERC5[1])
cor(as.numeric(PCs@scores[,3]),all_test_var_cont$GC_pct_miRNAaligned) # 0.5761681
plot(as.numeric(PCs@scores[,3]),all_test_var_cont$GC,pch=16,xlab='PC3',ylab='GC % all reads',col=colERC5[1])
cor(as.numeric(PCs@scores[,3]),all_test_var_cont$GC) # 0.0708655
plot(as.numeric(PCs@scores[,2]),all_test_var_cont$MeanReadLength,pch=16,xlab='PC2',ylab='mean read length',col=colERC5[1])
cor(as.numeric(PCs@scores[,2]),all_test_var_cont$MeanReadLength) # -0.6940441
dev.off()



###################################
###     Create table S1A        ###
################################### 

mir_TSS_annot = fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_TSS_annot_DeRie2017.txt',EVO_IMMUNO_POP))
mir_TSS_annot[,has_mirQTL:=hsa_ID%in%miRqtl_cis_annot$miRNA]
mir_TSS_annot[,conserved_TSS:=Conservation>.2]
mir_TSS_annot = mir_TSS_annot[order(hsa_ID, MIMAT),mget(c('hsa_ID','Nb_loci','MIMAT','mat_chr','mat_start','mat_end','mat_strand','MI_ID','hairpin_start','hairpin_end','arm','Promoter.name','TSS','Conservation'))]
fwrite(mir_TSS_annot,file=sprintf('%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable1A_mirRNA_annotation.tsv',EVO_IMMUNO_POP),sep='\t')

