
require(ggplot2)
luq = function(x){length(unique(x))}
colERC5 = c("#525252AA","#E31A1CAA","#33A02CAA","#1F78B4AA","#6A3D9AAA")
colERC = c("#969696", "#525252", "#FB9A99", "#E31A1C", "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4", "#CAB2D6", "#6A3D9A")
condIndex = c("NS","LPS","PAM3CSK4","R848","IAV" )


FPKM=fread(sprintf('%s/Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/FPKM_matrix.txt',EVO_IMMUNO_POP),sep='\t')
FPKM_melted=melt(FPKM)
colnames(FPKM_melted)=c('ID','sample','FPKM')
FPKM_melted[, individual := substr(sample, 1,6)]
FPKM_melted[, condition := condIndex[as.numeric(substr(sample, 8,8))]]
FPKM_melted[, population := substr(sample, 1,3)]
inverseNormalRankTransform=function(x){n=length(x);
									qnorm(frank(x,na.last=FALSE,ties.method='random')/(n+1),mean(x,na.rm=T),sd(x,na.rm=T))}
FPKM_melted[, FPKM_transformed := inverseNormalRankTransform(FPKM),by=.(ID)]
FPKM_melted[, condition := factor(condition,levels=condIndex,ordered=FALSE)]

eQTL_bestSNP=fread(sprintf("%s/Maxime/evo_immuno_pop_QTLs/eQTL/eQTL_assoc_bestSNP.txt",EVO_IMMUNO_POP))

####################################################
##	       find best QTL sharing model         	  ##
####################################################

make_modelisation_joint_eQTL <-function(data.table.line){
	temp = make_data_frame_to_study(data.table.line)
	temp[, condition := factor(condition,levels=condIndex)]

	allModels_STIM=cbind(NS=rep(0:1,16),
				LPS=rep(rep(0:1,e=16),1),
				PAM3CSK4=rep(rep(0:1,e=8),2),
				R848=rep(rep(0:1,e=4),4),
				IAV=rep(rep(0:1,e=2),8))
	res=list()
	for (num_model in 1:32){
		model_to_test=allModels_STIM[num_model,]
		temp[, condition_snp:= (genotype)*model_to_test[as.numeric(condition)]]
		snp_name = data.table.line$snps
		gene = data.table.line$gene
		if ( !all(temp$FPKM == 0)){#We test if the gene is truly expressed in that condition
    		mod = lm(temp, formula = FPKM_transformed ~ condition + population + condition_snp)
		    Likelihood = logLik(mod)
		    if( sum(model_to_test)>0 ){
        		Coeffs=summary(mod)$coeff
		        estimate=Coeffs["condition_snp","Estimate"]
        		pvalue=Coeffs["condition_snp","Pr(>|t|)"]
		    }else{
        		estimate=NA
		        pvalue=NA
        	}
		    res[[num_model]] = data.table(gene = gene, condition_test = paste(model_to_test,collapse=''), logLik = as.numeric(Likelihood), df = attr(Likelihood,'df'), Beta = estimate, pval = pvalue)
    	}else{
		    res[[num_model]] = data.table(gene = character(0), condition_test = character(0), logLik = numeric(0), df = numeric(0), Beta=numeric(0), pval=numeric(0))
   		}
   	}
   	res=rbindlist(res)
   	res[,ProbModel:=exp(logLik-min(logLik))/sum(exp(logLik-min(logLik)))]
   	res
}

make_data_frame_to_study <- function(data.table.line,alleles=F){
 snp_name = data.table.line$snps[1]
  gene = data.table.line$gene[1]
  data_temp=FPKM_melted[ID==gene,.(ID,individual,condition,population,FPKM_transformed,FPKM)]
  map=getMapInfo(snp_name)
  genotypes=getSNP(snp_name)
  if(alleles){
    alleles=c(map[,'allele.1'],map[,'allele.2'])
    if(nchar(alleles)[1]>1){
        alleles=c('ins','-')
    }
    if(nchar(alleles)[2]>1){
        alleles=c('del','-')
    }
    geno_allele=paste(c(alleles[1],alleles[2],alleles[2]),c(alleles[1],alleles[1],alleles[2]),sep='/')
    inds=names(genotypes)
    genotypes = geno_allele[1+ genotypes]
    data_temp[, genotype := factor(genotypes[match(data_temp$individual, inds)],rev(geno_allele))]
    }else{
    genotypes=2-genotypes
    data_temp[, genotype := genotypes[match(data_temp$individual, names(genotypes))]]
    }
#    data_temp=data_temp[which(!is.na(genotypes)),]
  return(data_temp)
}

######################
## Main computation ##
######################

results = list()
eQTL_bestSNP=eQTL_bestSNP[FDR<0.05,]
tim=Sys.time()
for (line_n in 1:nrow(eQTL_bestSNP)){
    print(paste(line_n))
    results[[paste(line_n)]] = make_modelisation_joint_eQTL(eQTL_bestSNP[line_n])
  }
print(Sys.time()-tim)

results = rbindlist(results)
fwrite(results, sprintf("%s/Maxime/miRNA_V2/data/15.eQTL_comparisons/sharing_eQTL_conditions_likelihoods.tsv",EVO_IMMUNO_POP), sep="\t")
# results=fread(sprintf("%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/sharing_mirQTL_conditions_likelihoods.tsv",EVO_IMMUNO_POP),colClass=c('character','character','numeric','numeric','numeric','numeric','numeric'))

bestModel = unique(results[order(gene, -logLik)],by= 'gene')
bestModel[,bestModel:=condition_test]

tab_NbCond_eQTL=table(sapply(strsplit(bestModel$bestModel,''),function(x){sum(x=='1')}))
tab_NbCond_eQTL=tab_NbCond_eQTL[as.character(1:5)]
names(tab_NbCond_eQTL)=1:5
tab_NbCond_eQTL[is.na(tab_NbCond_eQTL)]=0


###### piechart of mirQTL sharing
pieCol=c("#41B6C4","#A1DAB4","#FFFFB2","#FECC5C","#E31A1C")
par(mar=c(3,3,3,3))
tabpct=paste(round(100*tab_NbCond_eQTL/sum(tab_NbCond_eQTL), 1),'%')
pdf(sprintf("%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/pie_NbCond_eqtl.pdf",EVO_IMMUNO_POP),width=4,height=4)
pie(tab_NbCond_eQTL,col=pieCol,init.angle=90,labels=tabpct) # 5.5 x 5.5 inches
pie(tab_NbCond_eQTL,col=pieCol,init.angle=90,labels=rep(' ',5)) # 4 x 4 inches
dev.off()


meanFPKM=FPKM_melted[,.(meanFPKM=mean(FPKM)),by=.(ID,condition)]
ubigenes=meanFPKM[,.(ubiquitous=all(meanFPKM>log2(1+10))),by=ID]
ubigenes=ubigenes[ubiquitous==TRUE,ID]


tab_NbCond_eQTL_ubi=table(sapply(strsplit(bestModel$bestModel[bestModel$gene%in%ubigenes],''),function(x){sum(x=='1')}))
pdf(sprintf("%s/Maxime/miRNA_V2/figures/15.eQTL_comparisons/pie_NbCond_eqtl_ubiquitouslyexpressed.pdf",EVO_IMMUNO_POP),width=4,height=4)
tabpct=paste(round(100*tab_NbCond_eQTL_ubi/sum(tab_NbCond_eQTL_ubi), 1),'%')
pie(tab_NbCond_eQTL_ubi,col=pieCol,init.angle=90,labels=tabpct) # 5.5 x 5.5 inches
pie(tab_NbCond_eQTL_ubi,col=pieCol,init.angle=90,labels=rep(' ',5)) # 4 x 4 inches
dev.off()
