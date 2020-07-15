miR_QTL=fread(sprintf('%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/SupTable4A_mirQTLs_raw.tsv',EVO_IMMUNO_POP))

isomiR_QTL=fread(sprintf('%s/Maxime/miRNA_V2/data/11.response_miRQTL_and_other_analyses/SupTable4C_isomirQTLs_raw.tsv',EVO_IMMUNO_POP))

snps_data = fread(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/08.snp_data/general_informations/all_snps.tsv", sep=""))
snps_data[, snp_name := snp.name]
snps_data[, snp.name := NULL]

snp_data=fread(paste(EVO_IMMUNO_POP, "/Maxime/Evo_Immuno_pop_data/SNP_annotations/Map_imputed_essential_informations.txt", sep=""))
snp_data=merge(snps_data,snp_data)

miR_snps=snps_data[(MAF_EUB>.05 | MAF_AFB>.05) & (in_5p_miRNA | in_3p_miRNA)]


rankTransform = function(x){
        percentile=rank(x,ties.method='average')/(length(x)+1)
        mean_level=mean(x,na.rm=T)
        sd_level=sd(x,na.rm=T)
        qnorm(percentile,mean_level,sd_level)
        }




condIndex=c("NS","LPS","PAM3CSK4","R848","IAV")
options(stringsAsFactors=FALSE,max.print=9999)

pop='ALL'
CisDist=1e6
pvCis=0.05
minFreq=0.05
library(snpStats)
library(MatrixEQTL)

SampleAnnotFile=paste(EVO_IMMUNO_POP,'/Maxime/Evo_Immuno_pop_data/SampleAnnot.txt',sep='')
SampleAnnot=fread( SampleAnnotFile )
SampleAnnot=as.data.frame(SampleAnnot)

library(data.table)

Gene_annot=as.data.frame(fread(sprintf('%s/Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/GeneAnnot_expressed.txt',EVO_IMMUNO_POP)))
colnames(Gene_annot)[1:5]=c('gene_id','seqnames','strand','start','end')
Gene_annot$strand=ifelse(Gene_annot$strand>0,'+','-')

# read Gene data
Gene_data=fread(sprintf('%s/Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/FPKM_matrix.txt',EVO_IMMUNO_POP))
Gene_mat=as.matrix(as.data.frame(Gene_data)[grep('ENSG',Gene_data$ID),grep('AFB|EUB',colnames(Gene_data))])
rownames(Gene_mat)=Gene_data[grep('ENSG',Gene_data$ID),get('ID')]
rm(Gene_data);gc()

Feature_annot=Gene_annot[,c('gene_id','seqnames','start','end','strand')]
rownames(Feature_annot)=Feature_annot$gene_id
Feature_mat=Gene_mat

Feature_mat=Feature_mat[,SampleAnnot$sample_ID]
colnames(Feature_annot)[1]='feature_id'

RESperm=list()
for (perm in 0:100){
    cat(perm,':\n')
    if(perm==0){
        snplist=unique(c(miR_QTL$snps,isomiR_QTL$snps))
        Geno=getGenos(snplist)
        }
    else{
        resamp=snp_data[(MAF_EUB>.05 | MAF_AFB>.05) ,  sample(snp.name, sum(snp.name%in%snplist),replace=FALSE), by=cut(SNPfreq,seq(0,.5,.05))]
        Geno=getGenos(resamp$V1)
        }

    Map=Geno[,1:5]
    GenoNum=2-Geno[,-(1:5)]
    rownames(GenoNum)=Map$snp.name
    Map$SNPfreq=apply(GenoNum,1,mean,na.rm=T)/2

    TEST=FALSE
    pop='ALL'
    indiv=colnames(GenoNum)
    RES=NULL
    for(cond in 1:5){
            wSNP=which(Map$SNPfreq > 0)#[1:10000]
            if(TEST){wSNP=wSNP[1:10000]}
            PCs=read.table(sprintf('%s/02_ERC_Data/SampleCaseReports/EvoImmunoPop_%s_ancestryPCA.txt',HOME,pop),header=T)
            wCondPop=which(SampleAnnot$condition==condIndex[cond])
            indiv=indiv[indiv%in%SampleAnnot[wCondPop,'individu']]
            wCondPop=match(paste(sort(indiv),cond,sep='-'),SampleAnnot$sample_ID)

        GenoNum_Mat = SlicedData$new()
        GenoNum_Mat$CreateFromMatrix(as.matrix(GenoNum[wSNP,match(indiv,colnames(GenoNum))]))
        GenoNum_Mat$ResliceCombined(sliceSize = 5000)

        RPKM_Mat = SlicedData$new()
        RPKM_Mat$CreateFromMatrix(t(apply(Feature_mat[,wCondPop],1,rankTransform)))
        RPKM_Mat$ResliceCombined(sliceSize = 5000)
        show(RPKM_Mat)
    
        Cvrt=SlicedData$new()
        Cvrt$CreateFromMatrix(t(PCs[match(indiv,PCs$IID),3:4]))

        res=Matrix_eQTL_main(GenoNum_Mat,
                             RPKM_Mat,
                             Cvrt,
                             output_file_name=NULL,
                             pvOutputThreshold=1,
                             output_file_name.cis=NULL,
                             pvOutputThreshold.cis=1,
                             cisDist=CisDist,
                             snpspos=Map[wSNP,c("snp.name","chromosome","position")],
                             genepos=Feature_annot[,c('feature_id','seqnames','start','end')],
                             min.pv.by.genesnp=TRUE)

    #save(res,file=sprintf('%s/Maxime/miRNA_V2/data/16_trans_effects/transEffects_miR_cond%s.Rdata',EVO_IMMUNO_POP,cond))
        res$trans$eqtls$cond=cond
        RES=rbind(RES,data.table(res$trans$eqtls[,c('snps','gene','beta','pvalue','cond')]))
    }
    RES[,Perm:=perm]
    fwrite(RES, file=sprintf('%s/Maxime/miRNA_V2/data/16_trans_effects/perm/miRQTL_transEffects_%s.txt',EVO_IMMUNO_POP,perm))
#    RESperm[[perm+1]]=RES
}

#RESperm = rbindlist(RESperm)
#fwrite(RESperm, file=sprintf('%s/Maxime/miRNA_V2/data/16_trans_effects/perm/allPerm_miRQTL_transEffects.txt',EVO_IMMUNO_POP))

pi1=c()
for (perm in 0:100){
    cat(perm,'')
    trans_0=fread(sprintf('%s/Maxime/miRNA_V2/data/16_trans_effects/perm/miRQTL_transEffects_%s.txt',EVO_IMMUNO_POP,perm))
    pi1[perm+1]=1-pi0est(p = trans_0$p)$pi0
}

pi1
#  [1] 5.149773e-03 6.538302e-03 0.000000e+00 1.467121e-03 0.000000e+00
#  [6] 0.000000e+00 9.853807e-04 4.392215e-03 2.087316e-03 1.197675e-02
# [11] 0.000000e+00 0.000000e+00 4.614022e-03 0.000000e+00 0.000000e+00
# [16] 0.000000e+00 0.000000e+00 1.618161e-02 8.102693e-03 1.372052e-03
# [21] 0.000000e+00 2.331855e-02 0.000000e+00 4.563552e-03 9.423678e-04
# [26] 1.584876e-03 2.667240e-03 1.409373e-03 0.000000e+00 5.097479e-04
# [31] 0.000000e+00 9.719979e-03 3.318294e-03 0.000000e+00 1.445233e-02
# [36] 1.354946e-02 0.000000e+00 1.723956e-02 1.605995e-02 2.387260e-03
# [41] 1.264008e-02 0.000000e+00 6.449587e-03 1.426385e-03 0.000000e+00
# [46] 7.472382e-03 1.804666e-03 2.669939e-02 1.362624e-02 1.505778e-02
# [51] 0.000000e+00 1.357643e-02 0.000000e+00 1.132600e-02 0.000000e+00
# [56] 6.564793e-03 0.000000e+00 5.528671e-03 2.045293e-03 1.438578e-03
# [61] 0.000000e+00 3.755734e-03 7.279908e-03 5.112682e-03 0.000000e+00
# [66] 0.000000e+00 0.000000e+00 1.732266e-02 0.000000e+00 2.133541e-02
# [71] 1.957470e-03 0.000000e+00 2.026725e-02 6.214139e-03 1.043848e-02
# [76] 0.000000e+00 4.349861e-03 1.334765e-02 5.680720e-03 2.779909e-05
# [81] 0.000000e+00 5.833369e-03 1.337938e-02 2.288531e-02 5.012387e-03
# [86] 1.792375e-02 7.357817e-03 0.000000e+00 0.000000e+00 0.000000e+00
# [91] 6.249598e-03 0.000000e+00 3.878577e-03 0.000000e+00 0.000000e+00
# [96] 0.000000e+00 5.340314e-03 0.000000e+00 0.000000e+00 1.520660e-02
#[101] 0.000000e+00
> mean(pi1[-1]>pi1[1])
#[1] 0.36



miRNABS = fread(paste(EVO_IMMUNO_POP, "Martin/miRNA/11.miRNA_binding_sites/results/miRNA_binding_sites_on_protein_coding_3primeUTR_simplified.tsv", sep=""))
miRNABS[, code := paste(EnsemblGeneID, gsub("-", "_", miRNA))]