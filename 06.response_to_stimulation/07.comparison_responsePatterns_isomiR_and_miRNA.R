miR_QTL=fread(sprintf('%s/Maxime/miRNA_V2/data/09.mirQTL/cis_mirQTLs_with_FDR_filtered_best_miRNA_snp_association.tsv',EVO_IMMUNO_POP))
isomiR_QTL=fread(sprintf('%s/Maxime/miRNA_V2/data/10.isomiRQTL/cis_isomiRQTLs_with_FDR_filtered_best_isomiRNA_snp_association.tsv',EVO_IMMUNO_POP))

#isomiR_popDE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable3B_isomir_popDE.tsv",EVO_IMMUNO_POP),colClass=c('character','character','logical',rep('numeric',16),'character','numeric'))
#miR_popDE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable3A_mir_popDE.tsv",EVO_IMMUNO_POP),colClass=c('character',rep('numeric',16),'character','numeric'))
miR_popDE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4A_miR_popDE_V2.tsv",EVO_IMMUNO_POP),colClass=c('character',rep('numeric',20),'character','numeric'))
isomiR_popDE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4B_isomir_popDE_V2.tsv",EVO_IMMUNO_POP),colClass=c(rep('character',7),'logical',rep('numeric',20),'character','numeric','character'))
isomiR_popDE_table[,miRNA:=hsa_ID]

miR_DE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable2A_mirDE.tsv",EVO_IMMUNO_POP),colClass=c('character',rep('numeric',13),'character','numeric'))
isomiR_DE_table=fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable2B_V2_isomirDE.tsv",EVO_IMMUNO_POP),colClass=c('character','character','character','numeric','numeric','character','character','logical',rep('numeric',13),'character','numeric','character'))


merged=merge(isomiR_DE_table,miR_DE_table,by.x='hsa_ID',by.y='miRNA',suffix=c('.isomiR','.miR'))

table(merged$bestModel.isomiR,merged$bestModel.miR)

possible_models=unique(c(merged$bestModel.isomiR,merged$bestModel.miR))
OR=list()
for (i in possible_models){
    myOR=odds.ratio(table(merged$bestModel.isomiR==i,merged$bestModel.miR==i))
    myOR$model=i
    myOR$N_isomiR=sum(merged$bestModel.isomiR==i)
    myOR$N_miR=sum(merged$bestModel.miR==i)
    myOR$Pct_of_isomiRs_with_same_pattern_at_miRNA_level=mean(merged[bestModel.isomiR==i,bestModel.miR]==i)
    OR[[i]]=myOR
    }
    
    OR = rbindlist(OR)
    fwrite(OR,file=sprintf("%s/Maxime/miRNA_V2/data/06.response_to_stimulation/DE_miR_vs_isomiR_overlap.txt",EVO_IMMUNO_POP),sep='\t')
    
    mean(merged$bestModel.isomiR==merged$bestModel.miR)
#    0.3450464
    mean(merged[bestModel.isomiR!="" & bestModel.miR!="",bestModel.isomiR]== merged[bestModel.isomiR!="" & bestModel.miR!="",bestModel.miR])
#    0.3090024
table(ifelse(merged$bestModel.isomiR!="",'DE_isomiR', ''),ifelse(merged$bestModel.miR!="",'DE_miR',''))
#                DE_miR
#            453    458
#  DE_isomiR 316    822


#72% of DE isomiRs are observed for miRNA that change their 


tab=table(ifelse(unique(merged$hsa_ID)%in%merged[bestModel.isomiR!="",hsa_ID],'DE_isomiR',''),ifelse(unique(merged$hsa_ID)%in%merged[bestModel.miR!="",hsa_ID],'DE_miR',''))
tab
#                DE_miR
#             83     62
#  DE_isomiR  87    203

# 70% of miRNAs that change their isomiR ratio also change their expression (OR=3.1, P< 8.1 x 10-8), but only 30% of isomiRs display the same pattern of response across stimuli as their associated miRNA. 
# This percentage varied greatly according to the pattern of isomiR response, from up to 49% among isomiRs that respond in a TLR-specific manner to only 21% among isomiRs that display a viral-specific response.

odds.ratio(tab)
#            LowerCI       OR  UpperCI alpha            P
#odds ratio 2.021487 3.114816 4.827144  0.05 8.051599e-08
