
##############
##Libariries##
##############
require(ggplot2)
luq = function(x){length(unique(x))}

###########################
##### get GO info       ###
###########################

source(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/scripts/00a.GOannotation/01.Load_and_format_databases.R", sep=""))
GO_informations = load_BHF_UCL_miRNA()

########
##Data##
########
mir_DE = fread(sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/miRNA_differentially_expressed_in_conditions.tsv",EVO_IMMUNO_POP))
mir_DE=mir_DE[order(pvalue,-log2FC)]
mir_DE_diff= unique(mir_DE[fdr<.01 ,miRNA],by='miRNA')
mir_DE_diff_up= unique(mir_DE[fdr<.01 & log2FC>0,miRNA],by='miRNA')
mir_DE[,minP:=min(fdr),by=miRNA]
mir_DE[,maxFC:=max(abs(log2FC)),by=miRNA]

luq(mir_DE_diff) # 340

luq(unique(mir_DE[fdr<.01 ,miRNA],by='miRNA')) # 233

mir_DE[,.(signif=sum(fdr<.01),signifUp=sum(fdr<.01 & log2FC>0),signifDown=sum(fdr<.01 & log2FC<0),PctUp=sum(fdr<.01 & log2FC>0)/sum(fdr<.01)) ,by=cond]
#    condition_test signif signifUp signifDown     PctUp
# 2:            LPS    126       91         35 0.7222222
# 4:       PAM3CSK4    159      118         41 0.7421384
# 1:           R848    242      157         85 0.6487603
# 3:            IAV    214      124         90 0.5794393

######################################################
######       get best model informations         #####
######################################################

mir_DE_likelihood = fread(sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/miRNA_differentially_expressed_in_conditions_Likelihoods.tsv",EVO_IMMUNO_POP),colClass=c('character','character','numeric','numeric','numeric'))
mir_DE_likelihood[,bestModel:=condition_test]
mir_DE_likelihood[,condition_test:=NULL]
table(mir_DE_likelihood[ProbModel>.99, bestModel])
# 00001 00010 00011 00110 00111 01001 01011 01100 01101 01110 01111 
#    49    22    29     5     3     1     1     7     1    32    16 

nullModel = mir_DE_likelihood[bestModel=='00000',mget(c('miRNA','logLik'))][order(miRNA)]
bestModel = unique(mir_DE_likelihood[order(miRNA,-ProbModel)],by='miRNA')
bestModel=merge(bestModel,nullModel,suffixes=c('','.null'),by='miRNA',all.x=TRUE)
bestModel[,logDiff:=logLik-logLik.null]
bestModel[,LR_pvalue:=pchisq(2*logDiff,1,low=F)]
bestModel[,LR_fdr:=p.adjust(LR_pvalue,'fdr')]
bestModel=bestModel[,diff_miR:=miRNA%in%mir_DE_diff]

# get Best model excluding null 

# bestModel=bestModel[LR_fdr<.01]

#######################################################################
######       table with strong differences (annotated)            #####
#######################################################################

# mir_DE_lFC1=merge(mir_DE[log2FoldChange>1 & padj<0.01,],GO_informations[,mget(c('miRNA','GO_meaning'))],by='miRNA',all.x=TRUE)
# fwrite(mir_DE_lFC1,sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/paired/miRNA_differential_expression_detail_lFC1.tsv",EVO_IMMUNO_POP),sep='\t')


#######################################################################
######       table with all FC + best model 			           ####
#######################################################################

mir_DE_cast=dcast(mir_DE,miRNA+baseMean+minP+maxFC~condition_test,value.var=c('log2FC','pvalue','fdr'))

# add likelihood informations for DE miRs (other miRs are NA)

mir_DE_table = merge(mir_DE_cast,bestModel[which(diff_miR),mget(c("miRNA","bestModel","ProbModel","LR_fdr","diff_miR"))],by='miRNA',all.x=TRUE)

cols=c("miRNA","baseMean","minP", "maxFC","log2FC_LPS", "pvalue_LPS","fdr_LPS", 
			"log2FC_PAM3CSK4", "pvalue_PAM3CSK4",   "fdr_PAM3CSK4",
			"log2FC_R848","pvalue_R848",  "fdr_R848",
			 "log2FC_IAV","pvalue_IAV", "fdr_IAV","bestModel","ProbModel")
mir_DE_table=mir_DE_table[order(minP),mget(cols)]



# order by minP
mir_DE_table=mir_DE_table[order(minP)]
mir_DE_table[,minP:=NULL]

# write
fwrite(mir_DE_table,file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable2A_mirDE.tsv",EVO_IMMUNO_POP),sep='\t')


###################################################################
####    create a plot with sharing of miRNA across conditions   ####
####################################################################

Nb_miR=luq(mir_DE_table$miRNA)

allModels_STIM=cbind(NS=rep(0,8),
				LPS=rep(rep(0:1,e=8),1),
				PAM3CSK4=rep(rep(0:1,e=4),2),
				R848=rep(rep(0:1,e=2),4),
				IAV=rep(rep(0:1,e=1),8))
rownames(allModels_STIM)=apply(allModels_STIM,1,paste,collapse='')


code_values_ordered = c("01000","00100","00010","00001","01100","01010","01001",
                        "00110","00101","00011","01110","01101","01011","00111","01111")

pdf(sprintf("%s/Maxime/miRNA_V2/figures/06.response_to_stimulation/miRNA_differentially_expressed_in_conditions_Likelihoods_part1.pdf",EVO_IMMUNO_POP),height=3,width=5)
par(mar=c(4,7,1,1))
Image(t(allModels_STIM[code_values_ordered,-1]))
dev.off()

tab=table(mir_DE_table[, bestModel])
tab=tab[code_values_ordered]
tab[is.na(tab)]=0
names(tab)=code_values_ordered

pdf(sprintf("%s/Maxime/miRNA_V2/figures/06.response_to_stimulation/miRNA_differentially_expressed_in_conditions_Likelihoods_part2.pdf",EVO_IMMUNO_POP),height=3,width=5)
par(mar=c(4,8,3,3))
	allModels=allModels_STIM[code_values_ordered,-1]
    Prob_diff=sapply(2:5,function(i){sum(tab[substr(names(tab),i,i)==1])/Nb_miR})
    ExpectedPct=exp(allModels%*%log(matrix(Prob_diff,4))+(1-allModels)%*%log(1-matrix(Prob_diff,4)))
	x=barplot(tab,col='white',ylab='Model frequency',las=2)
	for (i in 2:5){
    	barplot(tab*allModels[,i-1],col=gsub('AA$','66',colERC5[i]),las=2,add=T,xaxt='n',yaxt='n')
	}
    x=barplot(tab,col='white',ylab='Model frequency',las=2)
	points(x,ExpectedPct*sum(tab),pch='-',col='grey',cex=3)
	for (i in 2:5){
    	barplot(tab*allModels[,i-1],col=gsub('AA$','66',colERC5[i]),las=2,add=T,xaxt='n',yaxt='n')
	}
    pbinom(as.numeric(tab),Nb_miR,as.numeric(ExpectedPct),low=F)
# 1.000000e+00 1.000000e+00 1.000000e+00 4.493919e-01 7.024789e-01 1.000000e+00 9.999934e-01 9.928960e-01 1.000000e+00 5.826827e-02 6.249472e-41 9.677411e-01 3.108674e-01 7.081362e-01 1.768179e-26
 dev.off()


#################################
##    load expression Data     ##
#################################

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

plot_miRNA_cond=function(miRNA,cond=1:5,add.violin=T){
    temp=mir_count_melted[ID==miRNA & condition %in% cond,]
    temp[,condition := factor(condIndex[as.numeric(condition)],levels=condIndex),]
    p <- ggplot(temp,aes(x=condition,y=count,fill=condition))+theme_bw()
    p <- p + scale_fill_manual(values=colERC5[cond]) + geom_jitter(width=0.2,size=0.5,colour='darkgrey')
  if(add.violin){
        p <- p + geom_violin() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot(fill="#FFFFFF88", outlier.size=0, notch=TRUE,width=0.4)
    }else{
        p <- p +  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot( outlier.size=0, notch=TRUE)
    }
    print(p)
}



#################################
##    plot miR response 1 by 1 ##
#################################
for (i in 1:length(mir_DE_diff)){
	mymiR=mir_DE_diff[i]
	pdf(sprintf("%s/Maxime/miRNA_V2/figures/06.response_to_stimulation/exemples/%s_%s-ALL.pdf",EVO_IMMUNO_POP,i,mymiR),height=2.3,width=5)
	plot_miRNA_cond(mymiR,1:5)
	dev.off()
}

#########################################
##   plot Number of Up/down regulated  ##
#########################################
to_plot=mir_DE[,.(signif=sum(fdr<.01),signifUp=sum(fdr<.01 & log2FC>0),signifDown=sum(fdr<.01 & log2FC<0),PctUp=sum(fdr<.01 & log2FC>0)/sum(fdr<.01)) ,by=condition_test]

to_plot_mat=t(as.matrix(to_plot[,mget(c('signifUp','signifDown'))]))
colnames(to_plot_mat)=to_plot[,condition_test]
to_plot_mat=to_plot_mat[,condIndex[-1]]

pdf(paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/figures/06.response_to_stimulation//number_miRNA_updownregulated_per_condition_global.pdf",sep=""),height=3.5,width=3)
par(mar=c(7,4,4,4))
barplot(to_plot_mat,beside=T,col=colERC[-(1:2)],las=2,space=c(0,.5),ylab='Number of miRNAs')
legend('top',fill=grey(c(.8,.3)),legend=c('up-regulated','down-regulated'),bty='n',inset=-.4,xpd=T)
dev.off()

