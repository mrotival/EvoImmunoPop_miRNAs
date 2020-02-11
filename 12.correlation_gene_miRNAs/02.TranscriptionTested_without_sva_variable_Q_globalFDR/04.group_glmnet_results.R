# sbatch -p geh --qos=geh -J TRtest_agg -o TRtest_agg.log --mem=99G /pasteur/homes/mrotival/01_scripts/R_runner_tars.sh /pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/12.correlation_gene_miRNAs/TranscriptionTested_without_sva_variable_Q_globalFDR/04.group_glmnet_results.R

###################################################################################################
##Processing of the results of correlation_allmiRNAs_geneFPKM_transcriptionTested_svaadjusted.R##
###################################################################################################
library(data.table)

all_results = list()
all_FP = list()
for (condition in 1:5){
  for (part in 0:99){
    perm=0
    print(paste(condition, part))
    all_results[[paste(condition, part,perm)]] = try(fread(paste(EVO_IMMUNO_POP,
                                                        "/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_transcriptionTested_withoutSVA/perm",perm,"/gene_expression_glmnet_alpha0.5_cond",condition,"_",part,".tsv", sep="")))
        if (class(all_results[[paste(condition, part,perm)]])=='try-error') { 
            all_results[[paste(condition, part,perm)]]=NULL
            }else{
            all_results[[paste(condition, part,perm)]][,Perm:=perm]
            }

    perm=1
    all_FP[[paste(condition, part,perm)]] = try(fread(paste(EVO_IMMUNO_POP,
                                                        "/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_transcriptionTested_withoutSVA/perm",perm,"/gene_expression_glmnet_alpha0.5_cond",condition,"_",part,".tsv", sep="")))
        if (class(all_FP[[paste(condition, part,perm)]])=='try-error') { 
            all_FP[[paste(condition, part,perm)]]=NULL
            }else{
            all_FP[[paste(condition, part,perm)]][,Perm:=perm]
            }
      }
}
all_results = rbindlist(all_results)
all_FP = rbindlist(all_FP)

########################
##  write all results ##
########################
fwrite(all_results[probs_selection>.5,],file=sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_transcriptionTested_withoutSVA/gene_expression_glmnet_alpha0.5_allConds_allQ_ProbOver0.5.tsv",EVO_IMMUNO_POP))

############################################################
##     compute FDR among miRNA/gene pairs, for all Q      ##
############################################################

tim=Sys.time()
FP=all_FP[grepl('hsa-',selected_variables),]
print(Sys.time()-tim);tim=Sys.time()
FP=unique(FP[order(-probs_selection),],by=c('Q','gene','selected_variables'))
print(Sys.time()-tim);tim=Sys.time()
FP=FP[,sapply(seq(0,1,by=.01),function(x){list(length(paste(gene,selected_variables)[probs_selection>x]))}),by=Q]
print(Sys.time()-tim);tim=Sys.time()

print(Sys.time()-tim);tim=Sys.time()
TP=all_results[grepl('hsa-',selected_variables),]
print(Sys.time()-tim);tim=Sys.time()
TP=unique(TP[order(-probs_selection),],by=c('Q','gene','selected_variables'))
print(Sys.time()-tim);tim=Sys.time()
TP=TP[,sapply(seq(0,1,by=.01),function(x){list(length(paste(gene,selected_variables)[probs_selection>x]))}),by=Q]
print(Sys.time()-tim);tim=Sys.time()

#FP=all_FP[grepl('hsa-',selected_variables), sapply(seq(0,1,by=.01),function(x){list(length(unique(paste(gene,selected_variables)[probs_selection>x])))}),by=Q]
#TP=all_results[grepl('hsa-',selected_variables), sapply(seq(0,1,by=.01),function(x){list(length(unique(paste(gene,selected_variables)[probs_selection>x])))}),by=Q]

FP_cutoff=melt(FP,id.var='Q',value.name='NbFP',variable.name='Cutoff')
FP_cutoff[,Cutoff:= (as.numeric(gsub('V','',Cutoff))-1)/100 ]

TP_cutoff=melt(TP,id.var='Q',value.name='NbTP',variable.name='Cutoff')
TP_cutoff[,Cutoff:= (as.numeric(gsub('V','',Cutoff))-1)/100 ]

FDR_cutoff=merge(TP_cutoff,FP_cutoff)
FDR_cutoff[,FDR:= cummin(NbFP/NbTP),by=Q]

# FDR_cutoff[Q==18 & Cutoff==0.47,]
#    Q Cutoff NbTP NbFP        FDR
#1: 18   0.47 2942  137 0.04656696
# CutOff_max
# Qmax

# all_results_max=all_results[probs_selection>CutOff_max & Q==Qmax,]
# all_results_max[substr(selected_variables,1,3)=='hsa',length(unique(gene))]
# [1] 10615


#  FDR for all Q, cutoffs and P-values
fwrite(FDR_cutoff,file=sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_transcriptionTested_withoutSVA/FDR_cutoff.txt",EVO_IMMUNO_POP))

FDR1=unique(FDR_cutoff[FDR<.01,],by='Q')
FDR5=unique(FDR_cutoff[FDR<.05,],by='Q')
FDR10=unique(FDR_cutoff[FDR<.1,],by='Q')
FDR_pair=rbind(FDR1,FDR5,FDR10)
FDR_pair[,FDR_class:= cut(FDR,c(-1,0.01,0.05,0.1))]
levels(FDR_pair$FDR_class)=c('<1%','<5%','<10%')

#  cutoffs and Nb of TP and FP ,for all Q :  at 1%, 5% and 10% FDR
fwrite(FDR_pair,file=sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_transcriptionTested_withoutSVA/FDR_pair.txt",EVO_IMMUNO_POP))

############################################################
##     compute FDR among genes with a miRNA, for all Q    ##
############################################################

FP=all_FP[grepl('hsa-',selected_variables),]
FP=unique(FP[order(-probs_selection),],by=c('Q','gene'))
FP=FP[,sapply(seq(0,1,by=.01),function(x){list(length(gene[probs_selection>x]))}),by=Q]

TP=all_results[grepl('hsa-',selected_variables),]
TP=unique(TP[order(-probs_selection),],by=c('Q','gene'))
TP=TP[,sapply(seq(0,1,by=.01),function(x){list(length(gene[probs_selection>x]))}),by=Q]

FP_cutoff_gene=melt(FP,id.var='Q',value.name='NbFP',variable.name='Cutoff')
FP_cutoff_gene[,Cutoff:= (as.numeric(gsub('V','',Cutoff))-1)/100 ]

TP_cutoff_gene=melt(TP,id.var='Q',value.name='NbTP',variable.name='Cutoff')
TP_cutoff_gene[,Cutoff:= (as.numeric(gsub('V','',Cutoff))-1)/100 ]

FDR_cutoff_gene=merge(TP_cutoff_gene,FP_cutoff_gene)
FDR_cutoff_gene[,FDR:= cummin(NbFP/NbTP),by=Q]
#  FDR for all Q, cutoffs and P-values
fwrite(FDR_cutoff_gene,file=sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_transcriptionTested_withoutSVA/FDR_cutoff_gene.txt",EVO_IMMUNO_POP))

FDR1_gene=unique(FDR_cutoff_gene[FDR<.01,],by='Q')
FDR5_gene=unique(FDR_cutoff_gene[FDR<.05,],by='Q')
FDR10_gene=unique(FDR_cutoff_gene[FDR<.1,],by='Q')
FDR_gene=rbind(FDR1_gene,FDR5_gene,FDR10_gene)

FDR_gene=rbind(FDR1_gene,FDR5_gene,FDR10_gene)
FDR_gene[,FDR_class:= cut(FDR,c(-1,0.01,0.05,0.1))]
levels(FDR_gene$FDR_class)=c('<1%','<5%','<10%')

#  cutoffs and Nb of TP and FP ,for all Q :  at 1%, 5% and 10% FDR
fwrite(FDR_gene,file=sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_transcriptionTested_withoutSVA/FDR_gene.txt",EVO_IMMUNO_POP))

########################################################################
##     write all results that pass a 5% FDR in pairs (miRNA-gene)     ##
########################################################################

imax=FDR_pair[FDR_class=='<5%' ,which.max(NbTP)]
Qmax=FDR_pair[which(FDR_class=='<5%')[imax],Q]
CutOff_max=FDR_pair[which(FDR_class=='<5%')[imax],Cutoff]
fwrite(all_results[probs_selection>CutOff_max & Q==Qmax,],file=sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_transcriptionTested_withoutSVA/gene_expression_glmnet_alpha0.5_allConds.tsv",EVO_IMMUNO_POP))


############################################################
##     plot FDR among miRNA/gene pairs, for all Q         ##
############################################################
library(ggplot2)

#FDR_cutoff_pair = fread(sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_transcriptionAdjusted_withoutSVA_old/FDR_cutoff.txt",EVO_IMMUNO_POP))
#FDR_pair = fread(sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_transcriptionAdjusted_withoutSVA_old/FDR_pair.txt",EVO_IMMUNO_POP))


# plot Q vs TP for various FDR thresholds
FDR_pair[,FDR:=factor(FDR_class,levels=c('<1%','<5%','<10%')),]
library(ggplot2)
pdf(sprintf("%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/variableQ_globalFDR/test_transcription_withoutSVA/FDR_pair.pdf",EVO_IMMUNO_POP),height=2.5,width=5)
p <- ggplot(FDR_pair,aes(x=Q,y= NbTP,col=FDR)) + geom_point() + geom_line()
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
p <- p + scale_colour_manual(values=brewer.pal(4,"YlOrRd")[2:4])
p <- p + xlab('Q') + ylab('Number of miRNA-gene associations') + xlim(c(0,60))+ ylim(c(0,max(FDR_pair$NbTP)))
print(p)
dev.off()

# plot Q vs Cutoff, for various FDR thresholds
pdf(sprintf("%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/variableQ_globalFDR/test_transcription_withoutSVA/cutoff_pair.pdf",EVO_IMMUNO_POP),height=2.5,width=5)
FDR_pair[,FDR:=factor(FDR_class,levels=c('<1%','<5%','<10%')),]
p <- ggplot(FDR_pair,aes(x=Q,y= Cutoff,col=FDR)) + geom_point() + geom_line()
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
p <- p + scale_colour_manual(values=brewer.pal(4,"YlGnBu")[2:4])
p <- p + xlab('Q') + ylab('Probability Cutoff') + xlim(c(0,60))+ ylim(c(0,1))
print(p)
dev.off()

# plot Cutoff VS FDR, for various values of Q
pdf(sprintf("%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/variableQ_globalFDR/test_transcription_withoutSVA/cutoff_pair_2.pdf",EVO_IMMUNO_POP),height=2.5,width=5)
FDR_cutoff = FDR_cutoff[Q%in%c(3,9,15,30),]
FDR_cutoff[,Q:=as.factor(Q)]
p <- ggplot(FDR_cutoff,aes(x=Cutoff,y= 100*FDR,col=Q)) + geom_line() + geom_point(size=.7)
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
p <- p + scale_colour_manual(values=brewer.pal(5,"YlGnBu")[2:5])
p <- p + ylab('FDR') + xlab('Probability Cutoff') + xlim(c(0,1))+ ylim(c(0,100))
print(p)
dev.off()


#################################################
##     plot FDR among genes, for all Q         ##
#################################################

#FDR_gene = fread(sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_transcriptionAdjusted_withoutSVA_old/FDR_gene.txt",EVO_IMMUNO_POP))
#FDR_cutoff_gene = fread(sprintf("%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_transcriptionAdjusted_withoutSVA_old/FDR_cutoff_gene.txt",EVO_IMMUNO_POP))


# plot Q vs TP, for various FDR thresholds
pdf(sprintf("%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/variableQ_globalFDR/test_transcription_withoutSVA/FDR_gene.pdf",EVO_IMMUNO_POP),height=2.5,width=5)
FDR_gene[,FDR:=factor(FDR_class,levels=c('<1%','<5%','<10%')),]
p <- ggplot(FDR_gene,aes(x=Q,y= NbTP,col=FDR)) + geom_point() + geom_line()
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
p <- p + scale_colour_manual(values=brewer.pal(4,"YlOrRd")[2:4])
p <- p + xlab('Q') + ylab('Number of gene with associated miRNAs') + xlim(c(0,60))+ ylim(c(0,max(FDR_gene$NbTP)))
print(p)
dev.off()


# plot Q vs Cutoff, for various FDR thresholds
pdf(sprintf("%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/variableQ_globalFDR/test_transcription_withoutSVA/cutoff_gene.pdf",EVO_IMMUNO_POP),height=2.5,width=5)
FDR_gene[,FDR:=factor(FDR_class,levels=c('<1%','<5%','<10%')),]
p <- ggplot(FDR_gene,aes(x=Q,y= Cutoff,col=FDR)) + geom_point() + geom_line()
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
p <- p + scale_colour_manual(values=brewer.pal(4,"YlGnBu")[2:4])
p <- p + xlab('Q') + ylab('Probability Cutoff') + xlim(c(0,60))+ ylim(c(0,1))
print(p)
dev.off()

# plot Cutoff VS FDR, for various values of Q
pdf(sprintf("%s/Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/variableQ_globalFDR/test_transcription_withoutSVA/cutoff_gene_2.pdf",EVO_IMMUNO_POP),height=2.5,width=5)
FDR_cutoff_gene = FDR_cutoff_gene[Q%in%c(3,9,15,30),]
FDR_cutoff_gene[,Q:=as.factor(Q)]
p <- ggplot(FDR_cutoff_gene,aes(x=Cutoff,y= 100*FDR,col=Q)) + geom_line() + geom_point(size=.7)
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=45, hjust=1))
p <- p + scale_colour_manual(values=brewer.pal(5,"YlGnBu")[2:5])
p <- p + ylab('FDR') + xlab('Probability Cutoff') + xlim(c(0,1))+ ylim(c(0,100))
print(p)
dev.off()


