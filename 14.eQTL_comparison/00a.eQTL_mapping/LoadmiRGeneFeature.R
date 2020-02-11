QTL_type='joint_eQTL'

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

# read miRNA data
corrected_count_transformed = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""), sep="\t")

#miRNA positions
miRNA_coordinate = fread(paste(EVO_IMMUNO_POP, "ERCPilot_SharedDBs/mirbase20/miRNA_mature_coordinates_strandinfo.bed", sep=""))
names(miRNA_coordinate) = c("chromosome", "start", "end", "miRNA_name", "V5", "V6", "V7", "V8")

## filter miRNAs with multiple sites on the genome
 miRNA_coordinate = miRNA_coordinate[miRNA_name %in% corrected_count_transformed$ID] #not enough
 duplicated_miRNAs = miRNA_coordinate[duplicated(miRNA_name), unique(miRNA_name)]
 miRNA_coordinate = miRNA_coordinate[!(miRNA_name %in% duplicated_miRNAs)]


miR_annot=miRNA_coordinate[,c('miRNA_name','chromosome','start','end','V6')]
rownames(miR_annot)=miR_annot$miRNA_name
colnames(miR_annot)=c('gene_id','seqnames','start','end','strand')
miR_annot$seqnames=gsub('chr','',miR_annot$seqnames)
miR_annot$seqnames[miR_annot$seqnames=='X']=23

miR_mat = corrected_count_transformed[match(miRNA_coordinate$miRNA_name,ID),]
colnames(miR_mat)=substr(colnames(miR_mat),1,8)

common=intersect(colnames(miR_mat),colnames(Gene_mat))
miR_ID=miR_mat$ID
miR_mat=as.matrix(miR_mat[,mget(common)])
rownames(miR_mat)=miR_ID

Gene_mat=Gene_mat[,common]

Feature_mat=rbind(miR_mat,Gene_mat)
Feature_annot=rbind(miR_annot,Gene_annot[,c('gene_id','seqnames','start','end','strand')])
