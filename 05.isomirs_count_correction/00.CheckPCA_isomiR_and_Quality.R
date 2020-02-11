

samples_names=unlist(fread(sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/sample_names_977_highQuality.tsv",EVO_IMMUNO_POP)))
all_test_var_cat = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/covariates/all_categorical_variables.tsv", sep=""))
names(all_test_var_cat)[1] = "library_ID"
all_test_var_cat = all_test_var_cat[match(samples_names, library_ID)]
all_test_var_cat[,condPop:=2*(condition-1)+ifelse(pop=='EUB',1,2)]
all_test_var_cat[,Machine_lane:=as.factor(paste(all_test_var_cat$machine, all_test_var_cat$lane,sep='_'))]

all_test_var_cont<-fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/covariates/all_continuous_variables_V2.0_MR.tsv", sep=""))
names(all_test_var_cont)[1] = "library_ID"
all_test_var_cont = all_test_var_cont[match(samples_names, library_ID)]


colERC=c("#969696", "#525252", "#FB9A99", "#E31A1C", "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4", "#CAB2D6", "#6A3D9A")

Pct_modif=fread(sprintf('%s/Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_modifs/Pct_modifs_perSample.txt',EVO_IMMUNO_POP))


library(RcolorBrewer)
nMachine=luq(substr(levels(as.factor(all_test_var_cat$machine)),1,22))
MachineCol=brewer.pal(nMachine,'Set3')
all_test_var_cat$machine_lane=paste(LETTERS[as.numeric(as.factor(all_test_var_cat$machine))],all_test_var_cat$lane,sep='_')
all_test_var_cat$machine_lane[all_test_var_cat$machine_lane=='L_99']='X_99'
#DF=data.frame(MeanQuality=all_test_var_cont$MeanQuality,
#        Machine_lane=all_test_var_cat$Machine_lane,
#        machine=all_test_var_cat$machine,
#        PC1_isomiR=PCs@scores[,1],
#        G_to_T_4thPos=Pct_modif$G_to_T_4thPos,
#        G_to_T_1stPos=Pct_modif$G_to_T_1stPos,
#        T_to_G_1stPos=Pct_modif$T_to_G_1stPos,
#        any_1stPos=Pct_modif$any_to_any_1stPos,
#        any_2ndPos=Pct_modif$any_to_any_2ndPos,
#        any_4thPos=Pct_modif$any_to_any_4thPos,
#        nta_3p_u=Pct_modif$nta_3p_u_miR,
#        nta_3p_a=Pct_modif$nta_3p_a_miR)
#
#
#DF=data.frame(MeanQuality=all_test_var_cont$MeanQuality,
#        Machine_lane=all_test_var_cat$Machine_lane,
#        machine=all_test_var_cat$machine,
#        PC1_isomiR=PCs@scores[,1],
#        G_to_T_4thPos=Pct_modif$G_to_T_4thPos,
#        G_to_T_1stPos=Pct_modif$G_to_T_1stPos,
#        T_to_G_1stPos=Pct_modif$T_to_G_1stPos,
#        any_1stPos=Pct_modif$any_to_any_1stPos,
#        any_2ndPos=Pct_modif$any_to_any_2ndPos,
#        any_4thPos=Pct_modif$any_to_any_4thPos,
#        nta_3p_u=Pct_modif$nta_3p_u_miR,
#        nta_3p_a=Pct_modif$nta_3p_a_miR)
#

DF=data.frame(MeanQuality=all_test_var_cont$MeanQuality,
        Machine_lane=all_test_var_cat$machine_lane,
        Sequencer=LETTERS[as.numeric(as.factor(all_test_var_cat$machine))],Pct_modif)
DF$Sequencer[DF$Sequencer=='L']='X'
library(lmer4)
library(MuMIn)
library(lmerTest)

R2_lane=c()
Var2test=c("MeanQuality", "nta_3p_o_miR","nta_3p_u_miR", "nta_3p_a_miR","shift_3p","shift_5p","shift_3p_template","shift_5p_template",colnames(Pct_modif)[grep('subs_Pos', colnames(Pct_modif))])
for( variable in Var2test){
    res=try(r.squaredGLMM(lmer(DF[,variable]~(1|Machine_lane),data=DF))[2])
    if(class(res)=='try-error'){
        R2_lane[variable]=0
    }else{
        R2_lane[variable]=res
        }
    }
fwrite(data.frame(var=Var2test,R2_lane),file=sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/QC_laneEffect/Variance_explained_by_laneEffect.txt',EVO_IMMUNO_POP),sep='\t')

dir.create(sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/QC_laneEffect',EVO_IMMUNO_POP))
pdf(sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/QC_laneEffect/MeanQuality_perlane.pdf',EVO_IMMUNO_POP),width=9,height=5)
p <- ggplot(DF,aes(y=MeanQuality,x=Machine_lane,fill=Sequencer))+geom_boxplot()
  p <- p + xlab('sequencer_lane')
  p <- p + ylab('Mean quality') 
  p <- p + theme_bw() + scale_fill_manual(values = MachineCol) 
  p <- p + theme(axis.text.x = element_text(angle = 90))
  p
dev.off()


    VarLabels=Var2test=c("Mean quality", "mean % of non template addition (C of G) across miRNAs","mean % of terminal Uridylation across miRNAs", "mean % of terminal Adenylation across miRNAs",
                                        "mean shift of 3' end site", "mean shift of 5' start site",
                                        "mean shift of template 3' end site","mean shift of template 5' start site",
                                        paste("mean % of substitution at position",1:10,"from start site"),paste("mean % of substitution at position",1:10,"from end site"))
                                        
for (i in 1:length(Var2test)){
    pdf(sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/QC_laneEffect/%s_perlane.pdf',EVO_IMMUNO_POP,Var2test[i]),width=9,height=5)
    p <- ggplot(DF,aes(y=DF[[Var2test[i]]],x=Machine_lane,fill=Sequencer))+geom_boxplot()
      p <- p + xlab('Sequencer lane')
      p <- p + ylab(VarLabels[i]) 
      p <- p + theme_bw() + scale_fill_manual(values = MachineCol) 
      p <- p + theme(axis.text.x = element_text(angle = 90))
      print(p)
    dev.off()
}

fwrite(DF,sprintf('%s/Maxime/miRNA_V2/figures/05.isomirs_count_correction/QC_laneEffect/data_figures.txt',EVO_IMMUNO_POP),sep='\t')
