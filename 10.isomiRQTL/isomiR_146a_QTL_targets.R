
############################################
### impact of isomiR-QTL on miRNA switch ###
############################################

targets=fread(sprintf('%s/Martin/miRNA/11.miRNA_binding_sites/results/common_isomiR_binding_sites_on_protein_coding_3primeUTR_simplified.tsv',EVO_IMMUNO_POP))
allGenes=unique(targets$EnsemblGeneID)

miR_146a_3p_targets=unique(targets[grep('146a-3p',targets$isomirID),],by=c('isomirID','GeneNames'))
miR_146a_3p_targets[,PctmiRNA_match:=as.numeric(gsub('%','',per1))]
miR_146a_3p_targets[,Pcttarget_match:=as.numeric(gsub('%','',per2))]
miR_146a_3p_targets_cast=dcast( miR_146a_3p_targets,EnsemblGeneID +GeneNames~gsub('(.*);_hsa.*','\\1',isomirID),value.var=c("PctmiRNA_match","Pcttarget_match"))
sum(!is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_-2_-2_+_GACCTCTGAAATTCAGTTCTTC_0') & !is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_0_-1_+_CCTGTGAAATTCAGTTCTTCA_1;4C->G')) # 873 targets in common
sum(!is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_-2_-2_+_GACCTCTGAAATTCAGTTCTTC_0') & is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_0_-1_+_CCTGTGAAATTCAGTTCTTCA_1;4C->G')) # 2273 targets lost 
sum(is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_-2_-2_+_GACCTCTGAAATTCAGTTCTTC_0') & !is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_0_-1_+_CCTGTGAAATTCAGTTCTTCA_1;4C->G')) # 2352 targets gained
# 873/(873+2273)

REF_targets=miR_146a_3p_targets_cast$EnsemblGeneID[!is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_-2_-2_+_GACCTCTGAAATTCAGTTCTTC_0')]
ALT_targets=miR_146a_3p_targets_cast$EnsemblGeneID[!is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_0_-1_+_CCTGTGAAATTCAGTTCTTCA_1;4C->G')]
length(unique(REF_targets)) #3146
length(unique(ALT_targets)) #3225
resGO=list()


getGO=function(set,background,setname,backgoundname,...){
    resGO=GOSeq(set,background,...)
    if(nrow(resGO)>0){
    resGO$set=rep(setname,nrow(resGO))
    resGO$background=rep(backgoundname,nrow(resGO))
    }
    resGO
}

resGO_all=rbindlist(resGO)
fwrite(resGO_all,file=sprintf('%s/Maxime/miRNA_V2/data/10.isomiRQTL/GO_enrich_miR-146a-3p_targets.txt',EVO_IMMUNO_POP),sep='\t')

resGO[['REF_all']]=getGO(REF_targets,allGenes,'REF_targets','allGenes',FDR=.2,overOnly=F)
resGO[['ALT_all']]=getGO(ALT_targets,allGenes,'ALT_targets','allGenes',FDR=.2,overOnly=F)
resGO[['lost_all']]=getGO(setdiff(REF_targets,ALT_targets),allGenes,'lost_targets','allGenes',FDR=.2,overOnly=F)
resGO[['gained_all']]=getGO(setdiff(ALT_targets,REF_targets),allGenes,'gained_targets','allGenes',FDR=.2,overOnly=F)
resGO[['common_all']]=getGO(intersect(ALT_targets,REF_targets),allGenes,'common_targets','allGenes',FDR=.2,overOnly=F)

resGO[['REF_expressed']]=getGO(REF_targets,rownames(FPKM_gene),'REF_targets','expressed',FDR=.2,overOnly=F)
resGO[['ALT_expressed']]=getGO(ALT_targets,rownames(FPKM_gene),'ALT_targets','expressed',FDR=.2,overOnly=F)
resGO[['lost_expressed']]=getGO(setdiff(REF_targets,ALT_targets),rownames(FPKM_gene),'lost_targets','expressed',FDR=.2,overOnly=F)
resGO[['gained_expressed']]=getGO(setdiff(ALT_targets,REF_targets),rownames(FPKM_gene),'gained_targets','expressed',FDR=.2,overOnly=F)
resGO[['common_expressed']]=getGO(intersect(ALT_targets,REF_targets),rownames(FPKM_gene),'common_targets','expressed',FDR=.2,overOnly=F)

miR_146a_3p_targets_80=miR_146a_3p_targets[PctmiRNA_match>80, ]

miR_146a_3p_targets_cast=dcast(miR_146a_3p_targets_80,EnsemblGeneID +GeneNames~gsub('(.*);_hsa.*','\\1',isomirID),value.var=c("PctmiRNA_match","Pcttarget_match"))
sum(!is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_-2_-2_+_GACCTCTGAAATTCAGTTCTTC_0') & !is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_0_-1_+_CCTGTGAAATTCAGTTCTTCA_1;4C->G')) # 54 targets in common
sum(!is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_-2_-2_+_GACCTCTGAAATTCAGTTCTTC_0') & is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_0_-1_+_CCTGTGAAATTCAGTTCTTCA_1;4C->G')) # 651 targets lost 
sum(is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_-2_-2_+_GACCTCTGAAATTCAGTTCTTC_0') & !is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_0_-1_+_CCTGTGAAATTCAGTTCTTCA_1;4C->G')) # 836 targets gained


REF_targets=miR_146a_3p_targets_cast$EnsemblGeneID[!is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_-2_-2_+_GACCTCTGAAATTCAGTTCTTC_0')]
ALT_targets=miR_146a_3p_targets_cast$EnsemblGeneID[!is.na(miR_146a_3p_targets_cast$'PctmiRNA_match_chr5_0_-1_+_CCTGTGAAATTCAGTTCTTCA_1;4C->G')]
resGO[['REF_80_all']]=getGO(REF_targets,allGenes,'REF_targets','allGenes')
resGO[['ALT_80_all']]=getGO(ALT_targets,allGenes,'ALT_targets','allGenes')
resGO[['lost_80_all']]=getGO(setdiff(REF_targets,ALT_targets),allGenes,'lost_targets','allGenes')
resGO[['gained_80_all']]=getGO(setdiff(ALT_targets,REF_targets),allGenes,'gained_targets','allGenes')
resGO[['common_80_all']]=getGO(intersect(ALT_targets,REF_targets),allGenes,'common_targets','allGenes')

resGO[['REF_80_expressed']]=getGO(REF_targets,rownames(FPKM_gene),'REF_targets','expressed')
resGO[['ALT_80_expressed']]=getGO(ALT_targets,rownames(FPKM_gene),'ALT_targets','expressed')
resGO[['lost_80_expressed']]=getGO(setdiff(REF_targets,ALT_targets),rownames(FPKM_gene),'lost_targets','expressed',...)
resGO[['gained_80_expressed']]=getGO(setdiff(ALT_targets,REF_targets),rownames(FPKM_gene),'gained_targets','expressed')
resGO[['common_80_expressed']]=getGO(intersect(ALT_targets,REF_targets),rownames(FPKM_gene),'common_targets','expressed')

Require("GeneAnnot")
Require("allGOterms")
GeneAnnot=as.data.table(GeneAnnot)
allGOterms=as.data.table(allGOterms)
GeneAnnot[,isImmune:=Ensembl.Gene.ID%in%allGOterms[go=='GO:0006955',gene]]

gained=GeneAnnot[match(setdiff(ALT_targets,REF_targets),GeneAnnot$Ensembl.Gene.ID),c('Associated.Gene.Name','Ensembl.Gene.ID','NS_mean','isImmune')][order(-NS_mean)]
lost=GeneAnnot[match(setdiff(REF_targets,ALT_targets),GeneAnnot$Ensembl.Gene.ID),c('Associated.Gene.Name','Ensembl.Gene.ID','NS_mean','isImmune')][order(-NS_mean)]

