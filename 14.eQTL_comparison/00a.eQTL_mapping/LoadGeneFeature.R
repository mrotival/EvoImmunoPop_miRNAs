QTL_type='gene_eQTL'

# read Gene Annot
library(data.table)

Gene_annot=as.data.frame(fread(sprintf('%s/Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/GeneAnnot_expressed.txt',EVO_IMMUNO_POP)))
colnames(Gene_annot)[1:5]=c('gene_id','seqnames','strand','start','end')
Gene_annot$strand=ifelse(Gene_annot$strand>0,'+','-')

# read Gene data
Gene_data=fread(sprintf('%s/Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/FPKM_matrix.txt',EVO_IMMUNO_POP))
Gene_mat=as.matrix(as.data.frame(Gene_data)[grep('ENSG',Gene_data$ID),grep('AFB|EUB',colnames(Gene_data))])
rownames(Gene_mat)=Gene_data[grep('ENSG',Gene_data$ID),get('ID')]
rm(Gene_data);gc()

Feature_mat=Gene_mat
Feature_annot=Gene_annot[,c('gene_id','seqnames','strand','start','end')]
