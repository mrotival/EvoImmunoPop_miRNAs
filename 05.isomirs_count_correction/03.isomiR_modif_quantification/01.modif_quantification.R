##################################################################################################################
###### This script quantifies isomiR moddifications it should be updated and rerun to get the final numbers ######
##################################################################################################################

library(data.table)
library(stringr)
mirFile=sprintf('%s/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv',EVO_IMMUNO_POP)

mirFileRaw=sprintf('%s/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_raw_counts.tsv',EVO_IMMUNO_POP)
mirAnnotFile=sprintf('%s/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_raw_counts.tsv',EVO_IMMUNO_POP)

isoMir=fread(sprintf('%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_counts.log2RPM.GCRL_Batch_lane_corrected.tsv',EVO_IMMUNO_POP))
samples_names=unlist(fread(sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/sample_names_977_highQuality.tsv",EVO_IMMUNO_POP)))
isoMir_mat=as.matrix(as.data.frame(isoMir[,mget(samples_names)]))
rownames(isoMir_mat)=isoMir$V1


######################
##Loading covariates##
######################
all_test_var_cat = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/covariates/all_categorical_variables.tsv", sep=""))
names(all_test_var_cat)[1] = "library_ID"
all_test_var_cat = all_test_var_cat[match(samples_names, library_ID)]

all_test_var_cont<-fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/covariates/all_continuous_variables_V2.0_MR.tsv", sep=""))
names(all_test_var_cont)[1] = "library_ID"
all_test_var_cont = all_test_var_cont[match(samples_names, library_ID)]

###### load miRNA data
miRNA_data=as.data.frame(fread(mirFile))
miRNA_mat=as.matrix(miRNA_data[,-1])
rownames(miRNA_mat)=miRNA_data[,1]

## load miRNA annot
miRNA_Annot=fread(mirAnnotFile)
miRNA_Annot = miRNA_Annot[,.(mir_1=miRNA_name, mir_2=miRNA_name2,hairpin= miRNA_hairpin)]
miRNA_Annot = miRNA_Annot[!duplicated(miRNA_Annot)]

isomiR_annot = as.data.frame(fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_commonOnly_V2.0.tsv',EVO_IMMUNO_POP)))
isomiR_annot=isomiR_annot[match(rownames(isoMir_mat),isomiR_annot$ID),]
mirID=isomiR_annot$mirID

miRNA_arm_miR=ifelse(grepl('3p_MIMAT',sort(unique(mirID))),'3p',ifelse(grepl('5p_MIMAT',sort(unique(mirID))),'5p',NA))

######### Assess the co-occurence of miRNA modifications 
# co-occurence of miRNA modifs
library(corrplot)
library(Hmisc)
COR_mir_modifs=rcorr(as.matrix(isomiR_annot[,which(sapply(isomiR_annot,is.numeric)|sapply(isomiR_annot,is.logical))]))

pdf(sprintf('%s/Maxime/miRNA_V2/figures/04.annotate_miRNAs&isomiRs/cooccurence_of_miRNA_modifs.pdf',EVO_IMMUNO_POP),width=12,height=12)
corrplot(COR_mir_modifs$r,p.mat=COR_mir_modifs$P,sig.level=c(1e-15,1e-10,1e-5),insig='label_sig',pch.cex=.8)
dev.off()

miRNA_modifs=c("miRNA_arm_num", "shift_3p",  "shift_5p","shift_3p_template", "shift_5p_template", 
"is_3p_extended", "is_3p_reduced", "is_3p_cannonical",
 "is_5p_extended", "is_5p_reduced", "is_5p_canonical", 
"is_3p_template_extended", "is_3p_template_reduced", "is_3p_template_canonical", 
 "is_5p_template_extended", "is_5p_template_reduced", "is_5p_template_canonical", 
 "nta_5p_miR", "nta_3p_miR", "nta_3p_a_miR", "nta_3p_u_miR", "nta_3p_o_miR", 
 "nsubs", 
"subs_Pos_5p_1", "subs_Pos_5p_2", "subs_Pos_5p_3", "subs_Pos_5p_4", "subs_Pos_5p_5", "subs_Pos_5p_6", "subs_Pos_5p_7", "subs_Pos_5p_8", "subs_Pos_5p_9", "subs_Pos_5p_10", 
"subs_Pos_3p_10", "subs_Pos_3p_9", "subs_Pos_3p_8", "subs_Pos_3p_7", "subs_Pos_3p_6", "subs_Pos_3p_5", "subs_Pos_3p_4", "subs_Pos_3p_3", "subs_Pos_3p_2", "subs_Pos_3p_1")

######### quantify the frequency of each type of miRNA modification, at any given miRNA and test for differences:
    # between conditions
    # between populations
    # between arms

Pct_modif_permiRNA=list()
RES_cond=list()
RES_pop=list()
RES_arm=list()
condIndex=c("NS", "LPS", "PAM3CSK4", "R848", "IAV")
condPopIndex=c("NS_AFB", "NS_EUB","LPS_AFB","LPS_EUB", "PAM3CSK4_AFB","PAM3CSK4_EUB","R848_AFB", "R848_EUB", "IAV_AFB", "IAV_EUB")

By=function(...){x=by(...);y=names(x);x=as.vector(x);names(x)=y;x}
 for (arm in c('any','3p','5p')){
    dir.create(sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/isomiR_modifs/%s',EVO_IMMUNO_POP,arm))
    dir.create(sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/isomiR_modifs/%s/population',EVO_IMMUNO_POP,arm))
    dir.create(sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/isomiR_modifs/%s/condition',EVO_IMMUNO_POP,arm))
    dir.create(sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/isomiR_modifs/%s/arm_differences_averaging_miRNAs',EVO_IMMUNO_POP,arm))
    dir.create(sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/isomiR_modifs/%s/arm_differences_averaging_samples',EVO_IMMUNO_POP,arm))
}
for (i in miRNA_modifs){
    cat(i,'')
    Pct_modif_permiRNA[[i]]=apply(2^isoMir_mat-1,2,function(x){by(x*isomiR_annot[,i],isomiR_annot$mirID,sum)/by(x,isomiR_annot$mirID,sum)})
    cat('computation done\n')
    fwrite(as.data.frame(Pct_modif_permiRNA[[i]]),file=sprintf('%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/data/%s_allmiRNAs.txt',EVO_IMMUNO_POP,i),sep='\t')
#    mean frequency of T->G substitutions at 2nd nucleotide 
    for (arm in c('any','3p','5p')){
        if(arm=='any'){
            meanFreq=apply(Pct_modif_permiRNA[[i]],2,mean,na.rm=T)
        }else{
            meanFreq=apply(Pct_modif_permiRNA[[i]][miRNA_arm_miR==arm,],2,mean,na.rm=T)
        }
     Mean_cond=By(meanFreq,condIndex[all_test_var_cat$cond],mean)
     Mean_cond=Mean_cond[condIndex]
     names(Mean_cond)=paste('Mean', names(Mean_cond),sep='_')
#    1. test every cond VS NS
     testCond=function(cond){
        P_cond=wilcox.test(meanFreq[all_test_var_cat$cond==1],meanFreq[all_test_var_cat$cond==cond])$p.value
        }
     P_cond=sapply(2:5,testCond)
     names(P_cond)=paste('P', condIndex[2:5],sep='_')
#    2. plot( cond differences )
     pdf(sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/isomiR_modifs/%s/condition/%s_%s_condition_differences.pdf',EVO_IMMUNO_POP,arm,i,arm),height=7,width=4)
     par(mar=c(10,4,1,1))
     boxplot(meanFreq~factor(condIndex[all_test_var_cat$cond],levels=condIndex),col=colERC5,notch=T,las=2)
     dev.off()
     RES_cond[[paste(i,arm,sep=':')]]=c(type=i,Mean_cond,P_cond,arm=arm)
     condPop=paste(condIndex[all_test_var_cat$cond],all_test_var_cat$pop,sep='_')
#    3. test popDiff in every Cond
     pdf(sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/isomiR_modifs/%s/population/%s_%s_population_differences.pdf',EVO_IMMUNO_POP,arm,i,arm),height=4,width=7)
     par(mar=c(10,4,1,1))
     boxplot(meanFreq~factor(condPop,levels=condPopIndex),col=colERC_AE,notch=T,las=2)
     dev.off()
#    4. plot( condPop differences )
     Mean_condpop=By(meanFreq,condPop,mean)[condPopIndex]
     names(Mean_condpop)=paste('Mean', condIndex,sep='_')
     testCond=function(cond){    
        P_cond=wilcox.test(meanFreq[all_test_var_cat$cond==cond & all_test_var_cat$pop=='AFB'],meanFreq[all_test_var_cat$cond==cond & all_test_var_cat$pop=='EUB'])$p.value
        }
     P_condpop=sapply(1:5,testCond)
     names(P_condpop)=paste('P', condIndex,sep='_')
     RES_pop[[paste(i,arm,sep=':')]]=c(type=i,Mean_condpop,P_condpop,arm=arm)
     if(arm=='any'){
#    5. test arm differences
        pdf(sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/isomiR_modifs/%s/arm_differences_averaging_miRNAs/arm_differences_%s_Averaging_miRNAs.pdf',EVO_IMMUNO_POP,arm,i),height=7,width=4)
        par(mar=c(10,4,1,1))
        MeanFreq_permiRNA=apply(Pct_modif_permiRNA[[i]],1,mean)
        boxplot(MeanFreq_permiRNA~miRNA_arm_miR,notch=T,col=brewer.pal(3,'Set3')[-2])
        dev.off()
        P_miRNA=wilcox.test(MeanFreq_permiRNA[miRNA_arm_miR=="3p"],MeanFreq_permiRNA[miRNA_arm_miR=="5p"])$p.value
        
        meanFreq_3p=apply(Pct_modif_permiRNA[[i]][miRNA_arm_miR=="3p",],2,mean,na.rm=T)
        meanFreq_5p=apply(Pct_modif_permiRNA[[i]][miRNA_arm_miR=="5p",],2,mean,na.rm=T)
        pdf(sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/isomiR_modifs/%s/arm_differences_averaging_samples/arm_differences_%s_Averaging_samples.pdf',EVO_IMMUNO_POP,arm,i),height=7,width=4)
        boxplot(list('3p'=meanFreq_3p,'5p'=meanFreq_5p),notch=T,col=brewer.pal(3,'Set3')[-2])
        dev.off()
        P_sample=wilcox.test(meanFreq_3p,meanFreq_5p)$p.value
        RES_arm[[i]]=c(type=i,meanFreq_3p=mean(meanFreq_3p),meanFreq_5p=mean(meanFreq_5p),P_sample,P_miRNA)
        }
    }
}
RES_pop=do.call(rbind,RES_pop)
RES_cond=do.call(rbind,RES_cond)
RES_arm=do.call(rbind,RES_arm)

RES_pop=as.data.frame(RES_pop)
RES_arm=as.data.frame(RES_arm)
RES_cond=as.data.frame(RES_cond)


RES_cond_melt=melt(as.data.table(RES_cond),id.vars=c('type','arm'),measure.vars=c('P_LPS','P_PAM3CSK4','P_R848','P_IAV'))
RES_cond_melt=RES_cond_melt[-grep('subs_Pos',type),]
RES_cond_melt$FDR=p.adjust(RES_cond_melt$value,'fdr')

for (i in grep('Mean|P',colnames(RES_cond))){ RES_cond[,i]=as.numeric(RES_cond[,i]) }

fwrite(RES_pop,file=sprintf('%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/tests/population_differences.txt',EVO_IMMUNO_POP),sep='\t')
fwrite(RES_cond,file=sprintf('%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/tests/stimulation_differences.txt',EVO_IMMUNO_POP),sep='\t')
fwrite(RES_arm,file=sprintf('%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/tests/arm_differences.txt',EVO_IMMUNO_POP),sep='\t')

Pct_modif_perSample=sapply(Pct_modif_permiRNA,apply,2,mean,na.rm=T)
fwrite(cbind(ID=rownames(Pct_modif_perSample),as.data.frame(Pct_modif_perSample)),file=sprintf('%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/Pct_modifs_perSample.txt',EVO_IMMUNO_POP),sep='\t')


    