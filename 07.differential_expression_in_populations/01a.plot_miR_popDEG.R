
##############
##Libariries##
##############
require(ggplot2)
luq = function(x){length(unique(x))}
colERC5 = c("#525252AA","#E31A1CAA","#33A02CAA","#1F78B4AA","#6A3D9AAA")
colERC = c("#969696", "#525252", "#FB9A99", "#E31A1C", "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4", "#CAB2D6", "#6A3D9A")
condIndex = c("NS","LPS","PAM3CSK4","R848","IAV" )

names(colERC5) = condIndex
names(colERC)=paste(rep(condIndex,e=2),'-',rep(c('EUB','AFB'),5))

###########################
##### get GO info       ###
###########################

source(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/scripts/00a.GOannotation/01.Load_and_format_databases.R", sep=""))
GO_informations = load_BHF_UCL_miRNA()

########
##Data##
########
mir_popDE = fread(sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/miRNA_differentially_expressed_between_populations.tsv",EVO_IMMUNO_POP))
mir_popDE=mir_popDE[order(pvalue,-log2FC_EminusA)]
mir_popDE_diff= unique(mir_popDE[fdr<.01 ,miRNA],by='miRNA')
#mir_popDE_diff_up= unique(mir_popDE[fdr<.01 & log2FC_EminusA>0,miRNA],by='miRNA')
mir_popDE[,minP:=min(fdr),by=miRNA]
mir_popDE[,maxFC:=max(abs(log2FC_EminusA)),by=miRNA]
mir_popDE[,condition_test:=factor(condition_test,levels=condIndex,ordered=TRUE)]
luq(mir_popDE_diff) # 244

######################################################
######       get best model informations         #####
######################################################

mir_popDE_likelihood = fread(sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/miRNA_differentially_expressed_between_populations_Likelihoods.tsv",EVO_IMMUNO_POP),colClass=c('character','character','numeric','numeric','numeric','numeric','numeric'))
mir_popDE_likelihood[,bestModel:=condition_test]
mir_popDE_likelihood[,condition_test:=NULL]
table(mir_popDE_likelihood[ProbModel>.99, bestModel])

nullModel = mir_popDE_likelihood[bestModel=='00000',mget(c('miRNA','logLik'))][order(miRNA)]
bestModel = unique(mir_popDE_likelihood[order(miRNA,-ProbModel)],by='miRNA')
bestModel=merge(bestModel,nullModel,suffixes=c('','.null'),by='miRNA',all.x=TRUE)
bestModel[,logDiff:=logLik-logLik.null]
bestModel[,LR_pvalue:=pchisq(2*logDiff,1,low=F)]
bestModel[,LR_fdr:=p.adjust(LR_pvalue,'fdr')]
bestModel=bestModel[,diff_miR:=miRNA%in%mir_popDE_diff]

mir_popDE_likelihood[ProbModel>.99 & bestModel!='11111', ]
#            miRNA    logLik df log2FC_EminusA         pval ProbModel bestModel
#1: hsa-miR-222-5p -699.6780  7     -0.4921022 1.525728e-30 0.9999852     01110
#2: hsa-miR-33b-3p -478.2233  7     -0.2648282 1.948098e-15 0.9956105     11100
#3:   hsa-miR-4773 -786.5139  7     -0.3906876 3.474624e-12 0.9992335     01100

# get Best model excluding null 
# bestModel=bestModel[LR_fdr<.01]

#######################################################################
######       table with strong differences (annotated)            #####
#######################################################################

# mir_popDE_lFC.2=merge(mir_popDE[log2FC_EminusA>.2 & fdr<0.01,],GO_informations[,mget(c('miRNA','GO_meaning'))],by='miRNA',all.x=TRUE)
# fwrite(mir_popDE_lFC1,sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/miRNA_differentially_expressed_between_populations_detail_lFC.2.tsv",EVO_IMMUNO_POP),sep='\t')
#mir_popDE_lFC.2[grep('inflam',mir_popDE_lFC.2$GO_meaning),]
#           miRNA condition_test  MeanAFB  MeanEUB log2FC_EminusA       pvalue          fdr        pbonf         minP     maxFC                                   GO_meaning
#1: hsa-miR-20a-5p             NS 9.710502 9.934890      0.2247284 1.545092e-08 3.505761e-07 4.860860e-05 3.505761e-07 0.2247284 negative regulation of inflammatory response
#2: hsa-miR-20a-5p            LPS 9.711057 9.915615      0.2025061 4.397037e-07 7.232287e-06 1.359124e-03 3.505761e-07 0.2247284 negative regulation of inflammatory response
#3: hsa-miR-590-3p       PAM3CSK4 3.946583 4.261816      0.3225051 1.861186e-11 8.874351e-10 5.996742e-08 8.874351e-10 0.3262782 negative regulation of inflammatory response
#4: hsa-miR-590-3p             NS 3.968366 4.292158      0.3262782 4.675313e-09 1.260802e-07 1.481607e-05 8.874351e-10 0.3262782 negative regulation of inflammatory response
#5: hsa-miR-590-3p           R848 3.894982 4.177980      0.2813603 1.950257e-07 3.487145e-06 6.059450e-04 8.874351e-10 0.3262782 negative regulation of inflammatory response
#6: hsa-miR-590-3p            LPS 3.955812 4.240473      0.2709563 1.224659e-06 1.692910e-05 3.738883e-03 8.874351e-10 0.3262782 negative regulation of inflammatory response
#7: hsa-miR-590-3p            IAV 4.011369 4.284814      0.2683177 2.258548e-06 2.902587e-05 6.854692e-03 8.874351e-10 0.3262782 negative regulation of inflammatory response


#######################################################################
######       table with all FC + best model 			           ####
#######################################################################

mir_popDE_cast=dcast(mir_popDE,miRNA + minP + maxFC ~ condition_test,value.var=c('MeanAFB','MeanEUB','log2FC_EminusA','pvalue','fdr'))

# add likelihood informations for DE miRs (other miRs are NA)

mir_popDE_table = merge(mir_popDE_cast,bestModel[which(diff_miR),mget(c("miRNA","bestModel","ProbModel","LR_fdr","diff_miR"))],by='miRNA',all.x=TRUE)

cols=c("miRNA","minP", "maxFC",
			"MeanAFB_NS", "MeanEUB_NS", "pvalue_NS", "fdr_NS",
			"MeanAFB_LPS", "MeanEUB_LPS", "pvalue_LPS", "fdr_LPS",
			"MeanAFB_PAM3CSK4", "MeanEUB_PAM3CSK4", "pvalue_PAM3CSK4", "fdr_PAM3CSK4",
			"MeanAFB_R848", "MeanEUB_R848", "pvalue_R848", "fdr_R848",
			"MeanAFB_IAV", "MeanEUB_IAV", "pvalue_IAV", "fdr_IAV",
			"bestModel", "ProbModel")

mir_popDE_table=mir_popDE_table[order(minP),mget(cols)]

# order by minP
mir_popDE_table=mir_popDE_table[order(minP)]
mir_popDE_table[,minP:=NULL]
mir_popDE_table[,maxFC:=NULL]

# write
fwrite(mir_popDE_table,file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4A_miR_popDE_V2.tsv",EVO_IMMUNO_POP),sep='\t')
mir_popDE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4A_miR_popDE_V2.tsv",EVO_IMMUNO_POP),colClass=c('character',rep('numeric',20),'character','numeric'))

######################################################################
######       table with miR response + best model 			       ####
#######################################################################

mir_popDE_resp=fread(sprintf("%s/Maxime/miRNA_V2/data/07.differential_expression_in_populations/miRNA_differential_response_between_populations.tsv",EVO_IMMUNO_POP))
mir_popDE_resp[,minP:=min(fdr),by=miRNA]
mir_popDE_resp[,maxFC:=max(abs(log2FC_EminusA)),by=miRNA]
mir_popDE_resp[,condition_test:=factor(condition_test,levels=condIndex[-1],ordered=TRUE)]

mir_popDE_resp_cast=dcast(mir_popDE_resp,miRNA + minP + maxFC ~ condition_test,value.var=c('MeanAFB','MeanEUB','log2FC_EminusA','pvalue','fdr'))

mir_popDE_resp_table = merge(mir_popDE_resp_cast,bestModel[which(diff_miR),mget(c("miRNA","bestModel","ProbModel","LR_fdr","diff_miR"))],by='miRNA',all.x=TRUE)
mir_popDE_resp_table[order(minP),]
fwrite(mir_popDE_resp_table[order(minP),],file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4C_miR_popDE_resp.tsv",EVO_IMMUNO_POP),sep='\t')


cols=c("miRNA","MeanAFB_LPS", "MeanEUB_LPS", "pvalue_LPS", "fdr_LPS",
			"MeanAFB_PAM3CSK4", "MeanEUB_PAM3CSK4", "pvalue_PAM3CSK4", "fdr_PAM3CSK4",
			"MeanAFB_R848", "MeanEUB_R848", "pvalue_R848", "fdr_R848",
			"MeanAFB_IAV", "MeanEUB_IAV", "pvalue_IAV", "fdr_IAV",
			"bestModel", "ProbModel")
			
mir_popDE_resp_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4C_miR_popDE_resp.tsv",EVO_IMMUNO_POP),sep='\t')

fwrite(mir_popDE_resp_table[order(minP),mget(cols)],file=sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4C_miR_popDE_resp_V2.tsv",EVO_IMMUNO_POP),sep='\t')

###################################################################
####    create a plot with sharing of miRNA across conditions   ####
####################################################################

Nb_miR=luq(mir_popDE_table$miRNA)

allModels_STIM=cbind(NS=rep(0:1,16),
				LPS=rep(rep(0:1,e=16),1),
				PAM3CSK4=rep(rep(0:1,e=8),2),
				R848=rep(rep(0:1,e=4),4),
				IAV=rep(rep(0:1,e=2),8))
rownames(allModels_STIM)=apply(allModels_STIM,1,paste,collapse='')

code_values_ordered = c("00000","10000","01000","00100","00010","00001",
                                "11000","10100","10010","10001",
                                "01100","01010","01001","00110","00101","00011",
                                "11100","11010","11001","10110","10101","10011",
                                "01110","01101","01011","00111",
                                "01111","10111","11011","11101","11110","11111")

pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/miRNA_differentially_expressed_between_populations_Likelihoods_part1.pdf",EVO_IMMUNO_POP),height=3,width=5)
par(mar=c(4,7,1,1))
Image(t(allModels_STIM[code_values_ordered,]))
dev.off()

tab=table(mir_popDE_table[, bestModel])
tab=tab[code_values_ordered]
tab[is.na(tab)]=0
names(tab)=code_values_ordered


#### Check this: how do we plot sharing ? 
# 1. number of conditions + which stimulations among stimulation specific ?  

pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/miRNA_differentially_expressed_between_populations_Likelihoods_part2.pdf",EVO_IMMUNO_POP),height=3,width=5)
par(mar=c(4,8,3,3))
	allModels=allModels_STIM[code_values_ordered,]
    Prob_diff=sapply(1:5,function(i){sum(tab[substr(names(tab),i,i)==1])/Nb_miR})
    ExpectedPct=exp(allModels%*%log(matrix(Prob_diff,5))+(1-allModels)%*%log(1-matrix(Prob_diff,5)))
	x=barplot(tab,col='white',ylab='Model frequency',las=2)
	for (i in 1:5){
    	barplot(tab*allModels[,i],col=gsub('AA$','66',colERC5[i]),las=2,add=T,xaxt='n',yaxt='n')
	}
    x=barplot(tab,col='white',ylab='Model frequency',las=2)
	points(x,ExpectedPct*sum(tab),pch='-',col='grey',cex=3)
	for (i in 1:5){
    	barplot(tab*allModels[,i],col=gsub('AA$','66',colERC5[i]),las=2,add=T,xaxt='n',yaxt='n')
	}
    pbinom(as.numeric(tab),Nb_miR,as.numeric(ExpectedPct),low=F)

 dev.off()

# create pie charts of sharing between condition 

tab_NbCond_miR_popDE=table(apply(allModels_STIM,1,sum)[mir_popDE_table$bestModel])
#  1   2   3   4   5 
#  5   8  24  44 163 

pieCol=c("#41B6C4","#A1DAB4","#FFFFB2","#FECC5C","#E31A1C")
par(mar=c(3,3,3,3))
tabpct=paste(round(100*tab_NbCond_miR_popDE/sum(tab_NbCond_miR_popDE), 1),'%')
pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/pie_NbCond_miR_popDE.pdf",EVO_IMMUNO_POP),width=4,height=4)
pie(tab_NbCond_miR_popDE,col=pieCol,init.angle=90,labels=tabpct) # 5.5 x 5.5 inches
pie(tab_NbCond_miR_popDE,col=pieCol,init.angle=90,labels=rep(' ',5)) # 4 x 4 inches
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
mir_count_melted[, condition := condIndex[as.numeric(substr(sample, 8,8))]]
mir_count_melted[, population := substr(sample, 1,3)]
inverseNormalRankTransform=function(x){n=length(x);
									qnorm(frank(x,na.last=FALSE,ties.method='random')/(n+1),mean(x,na.rm=T),sd(x,na.rm=T))}
mir_count_melted[, count_transformed := inverseNormalRankTransform(count),by=.(ID)]
mir_count_melted[, condition := factor(condition,levels=condIndex,ordered=TRUE)]

########################################
##    plot population differences     ##
########################################

plot_miRNA_pop=function(miRNA,cond=1:5,add.violin=T){
    temp=mir_count_melted[ID==miRNA & condition %in% condIndex[cond],]
    temp[,condition := factor(condition,levels=condIndex[cond],ordered=TRUE),]

    temp[,cond_pop := factor(paste(condition,'-',population),levels=paste(rep(condIndex,e=2),'-',rep(c('AFB','EUB'),5)))]
    # plot
    p <- ggplot(temp,aes(x=population,y=count,fill=cond_pop))+theme_bw()+facet_grid(~condition)
    p <- p + scale_fill_manual(values=colERC) 
  if(add.violin){
        p <- p + geom_violin() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot(fill="#FFFFFF88", outlier.size=0, notch=TRUE,width=0.4)
    }else{
        p <- p +  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
        p <- p + geom_boxplot( outlier.size=0, notch=TRUE)
    }
    p <- p + geom_jitter(width=0.2,size=0.5,colour="#52525266")
    print(p)
}
#################################
##    plot miR popDE 1 by 1    ##
#################################
for (i in 1:length(mir_popDE_diff)){
	mymiR=mir_popDE_diff[i]
	cat(i,mymiR,'\n')
	pdf(sprintf("%s/Maxime/miRNA_V2/figures/07.differential_expression_in_populations/exemples/%s_%s-ALL.pdf",EVO_IMMUNO_POP,i,mymiR),height=2.8,width=7)
	plot_miRNA_pop(mymiR,1:5)
	dev.off()
}

###############################################################
##   plot Number of higher in europeans/ africans regulated  ##
###############################################################
Nb_popDE = mir_popDE[,.(signif=sum(fdr<.01),signif_higherEUB=sum(fdr<.01 & log2FC_EminusA>0),signif_higherAFB=sum(fdr<.01 & log2FC_EminusA<0),Pct_higherEUB=sum(fdr<.01 & log2FC_EminusA>0)/sum(fdr<.01)) ,by=condition_test][order(condition_test)]
Nb_popDE[,pval:=pbinom(signif_higherEUB ,signif,0.5,low=F)]
#   condition_test signif signif_higherEUB signif_higherAFB Pct_higherEUB         pval
#1:             NS    133               75               58     0.5639098 0.0591200552
#2:            LPS    133               65               68     0.4887218 0.5687982544
#3:       PAM3CSK4    134               63               71     0.4701493 0.7272328168
#4:           R848    152               80               72     0.5263158 0.2327544991
#5:            IAV    153               96               57     0.6274510 0.0005752276

to_plot_mat=t(as.matrix(Nb_popDE[,mget(c('signif_higherAFB','signif_higherEUB'))]))
colnames(to_plot_mat)=Nb_popDE[,condition_test]

pdf(paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/figures/07.differential_expression_in_populations/number_miRNA_HigherInEUBorAFB_per_condition_global.pdf",sep=""),height=3.5,width=3)
par(mar=c(7,4,4,4))
barplot(to_plot_mat,beside=T,col=colERC[paste(rep(condIndex,e=2),rep(c('AFB','EUB'),5),sep=' - ')],las=2,space=c(0,.5),ylab='Number of miRNAs')
legend('top',fill=grey(c(.3,.8)),legend=c('Higher in Africans','Higher in Europeans'),bty='n',inset=-.4,xpd=T)
dev.off()

