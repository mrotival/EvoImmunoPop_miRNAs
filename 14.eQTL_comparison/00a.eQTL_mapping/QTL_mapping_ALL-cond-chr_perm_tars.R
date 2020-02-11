
TEST=FALSE

condIndex=c("NS","LPS","PAM3CSK4","R848","IAV")
options(stringsAsFactors=FALSE,max.print=9999)

pop='ALL'
CisDist=1e6
pvCis=0.05
minFreq=0.05
library(snpStats)
library(MatrixEQTL)

if(!is.null(CHR)){
	Geno=read.plink(paste(EVO_IMMUNO_POP,'/DATA_FREEZE/ERC_Main_Genotyping_24022015/Imputation/EvoImmunoPop_imputation_200x19619457_chr',CHR,'.bed',sep=''))
	GenoNum=as(Geno$genotype,'numeric') 
	load(file=paste(EVO_IMMUNO_POP,'/Maxime/SNP_annotation/imputed_Genotypes/MapFiles/Map_imputation_200x19619457_chr',CHR,'.Rdata',sep=''))
	rm(Geno);gc()
	chr.dist=5e7
	toRemove=read.table(paste(EVO_IMMUNO_POP,'/Maxime/SNP_annotation/imputed_Genotypes/multiallelic_Sites.txt',sep=''),header=T)
	GenoNum=2-GenoNum[,!Map$snp.name%in%toRemove$snp.name]
	Map=Map[which(!Map$snp.name%in%toRemove$snp.name),]
	Map$MAF_AFB=pmin(Map$SNPfreq_AF,1-Map$SNPfreq_AF)
	Map$MAF_EUB=pmin(Map$SNPfreq_EU,1-Map$SNPfreq_EU)
}

SampleAnnotFile=paste(EVO_IMMUNO_POP,'/Maxime/Evo_Immuno_pop_data/SampleAnnot.txt',sep='')
SampleAnnot=fread( SampleAnnotFile )
SampleAnnot=as.data.frame(SampleAnnot)
SampleAnnot=SampleAnnot[ SampleAnnot$sample_ID %in% colnames(Feature_mat) , ]

Feature_mat=Feature_mat[,SampleAnnot$sample_ID]
Feature_annot=as.data.frame(Feature_annot)
colnames(Feature_annot)[1]='feature_id'

TSS=ifelse(Feature_annot$strand=='+',Feature_annot$start,Feature_annot$end)
Feature_annot$start=TSS
Feature_annot$end=TSS


########################################################################
##				Genes by Pop - Standard - Cis 						  ##
########################################################################

dir.create(paste(EVO_IMMUNO_POP,'/Maxime/miRNA_V2/data/15.eQTL_comparisons/eQTL_mapping/',QTL_type,sep=''))
dir.create(paste(EVO_IMMUNO_POP,'/Maxime/miRNA_V2/data/15.eQTL_comparisons/eQTL_mapping/',QTL_type,'/permutations',sep=''))
dir.create(paste(EVO_IMMUNO_POP,'/Maxime/miRNA_V2/data/15.eQTL_comparisons/eQTL_mapping/',QTL_type,'/permutations/perm', perm, sep=''))

RESCIS=list()
QTL_pvalues=list()
tim=Sys.time()
indiv=sort(unique(SampleAnnot$individu))
if(perm>0){
	set.seed(perm)
	indiv[grep('AFB',indiv)]=sample(indiv[grep('AFB',indiv)])
	indiv[grep('EUB',indiv)]=sample(indiv[grep('EUB',indiv)])
	}

for(cond in 1:5){
	if(pop=='ALL'){
		wSNP=which(Map$MAF_AFB > minFreq | Map$MAF_EUB > minFreq)#[1:10000]
		if(TEST){wSNP=wSNP[1:10000]}
		PCs=read.table(paste(HOME,'/02_ERC_Data/SampleCaseReports/EvoImmunoPop_',pop,'_ancestryPCA.txt',sep=''),header=T)
        PCs$population=ifelse(substr(PCs$IID,1,3) == "EUB", 0, 1)
		wCondPop=which(SampleAnnot$condition==condIndex[cond])
		indiv=indiv[indiv%in%SampleAnnot[wCondPop,'individu']]
		wCondPop=match(paste(sort(indiv),cond,sep='-'),SampleAnnot$sample_ID)
	}
	GenoNum_Mat = SlicedData$new()
	GenoNum_Mat$CreateFromMatrix(t(GenoNum)[wSNP,match(indiv,rownames(GenoNum))])
	GenoNum_Mat$ResliceCombined(sliceSize = 5000)
    rankTransform = function(x){
        percentile=rank(x,ties.method='average')/(length(x)+1)
        mean_level=mean(x,na.rm=T)
        sd_level=sd(x,na.rm=T)
        qnorm(percentile,mean_level,sd_level)
        }

	RPKM_Mat = SlicedData$new()
	RPKM_Mat$CreateFromMatrix(t(apply(Feature_mat[,wCondPop],1,rankTransform)))
	RPKM_Mat$ResliceCombined(sliceSize = 5000)
	show(RPKM_Mat)
	
	Cvrt=SlicedData$new()
	Cvrt$CreateFromMatrix(t(PCs[match(indiv,PCs$IID),5,drop=F]))

	res=Matrix_eQTL_main(GenoNum_Mat,
	                     RPKM_Mat,
	                     Cvrt,
	                     output_file_name=NULL,
	                     pvOutputThreshold=0,
	                     output_file_name.cis=NULL,
	                     pvOutputThreshold.cis=pvCis,
	                     cisDist=CisDist,
	                     snpspos=Map[wSNP,c("snp.name","chromosome","position")],
	                     genepos=Feature_annot[,c('feature_id','seqnames','start','end')],
	                     min.pv.by.genesnp=TRUE)

	if(nrow(res$cis$eqtls)>0){
		res$cis$eqtls$condition=condIndex[cond]
		res$cis$eqtls$population=pop
		res$cis$eqtls$isCis=TRUE
		RESCIS[[cond]]=res$cis$eqtls
	}
    Pvalues=merge(data.table(feature_id=names(res$cis$min.pv.gene),
                                    pvalue=res$cis$min.pv.gene,
                                    condition=condIndex[cond],
                                    QTL_type=QTL_type),Feature_annot)
    QTL_pvalues[[cond]]=Pvalues[seqnames==CHR,mget(c('feature_id','pvalue','condition','QTL_type'))]
	print(Sys.time()-tim)
}

QTL_pvalues=rbindlist(QTL_pvalues)
#QTL_pvalues=dcast(QTL_pvalues,feature_id+QTL_type~condition,value.var='pvalue')
#QTL_pvalues[,Min_Pvalue:=min(pvalue_NS,pvalue_LPS,pvalue_PAM3CSK4,pvalue_R848,pvalue_IAV))
fwrite(QTL_pvalues,file=paste(EVO_IMMUNO_POP,'/Maxime/miRNA_V2/data/15.eQTL_comparisons/eQTL_mapping/',QTL_type,'/permutations/perm',perm,'/Cis-',QTL_type,'_ALL_allCond_chr',CHR,'_perm',perm,'_bestP_perFeature.txt',sep=''),sep='\t')

if(perm==0){
    RESCIS=do.call(rbind,RESCIS)
    RESCIS$snps=as.character(RESCIS$snps)
    RESCIS$gene=as.character(RESCIS$gene)
    RESCIS$QTL_type=QTL_type
    wSNP=match(RESCIS$snp,rownames(Map))
    wFeature=match(RESCIS$gene,rownames(Feature_annot))
    RESCIS$CisDist=ifelse(Map[wSNP,'chromosome']!=Feature_annot[wFeature,'seqnames'],Inf,
    		ifelse(Feature_annot[wFeature,'strand']=='+', 1, -1)* ifelse(Map[wSNP,'position'] < Feature_annot[wFeature,'start'], 
    															Feature_annot[wFeature,'start']-Map[wSNP,'position'],
											    				ifelse(Feature_annot[wFeature,'end']<Map[wSNP,'position'],
											    					Map[wSNP,'position']-Feature_annot[wFeature,'end'],0))
    		)
	fwrite(RESCIS,file=paste(EVO_IMMUNO_POP,'/Maxime/miRNA_V2/data/15.eQTL_comparisons/eQTL_mapping/',QTL_type,'/permutations/perm',perm,'/Cis-',QTL_type,'_ALL_allCond_chr',CHR,'_perm',perm,'_asssoc.txt',sep=''),sep='\t')
	}

q('no')