####################################################################################################
##Aim : compute correlation between miRNA expression and gene expression adjusted on transcription##
####################################################################################################
require(Hmisc)
require(sva)
require(MatrixEQTL)

#############
##Arguments##
#############
# args = commandArgs(trailingOnly=TRUE)
# cond_to_investigate = as.numeric(args[1])
#miRNA_number_to_investigate = as.numeric(args[2])

inverseNormalRankTransform=function(x){n=length(x);
									qnorm(frank(x,na.last=FALSE,ties.method='random')/(n+1),mean(x,na.rm=T),sd(x,na.rm=T))}

#############
##Load data##
#############

miRNA_expression= fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""))
names(miRNA_expression)[1] = "miRNA"
miRNA_expression = melt(miRNA_expression, id = "miRNA", variable.name = "sample", value.name = "expression")
miRNA_expression[, individual := substr(sample, 1,6)]
miRNA_expression[, condition := substr(sample, 8,8)]
miRNA_expression[, population := substr(sample, 1,3)]
miRNA_expression[, expression_inrt := inverseNormalRankTransform(expression),by=.(miRNA,condition)]

# miRNA_studied = sort(unique(miRNA_expression$miRNA))
# miRNA_studied = miRNA_studied[miRNA_number_to_investigate]
# miRNA_expression = miRNA_expression[miRNA == miRNA_studied]

 genes_expression = fread(paste(EVO_IMMUNO_POP, "Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/FPKM_matrix.txt", sep=""))
 genes_expression = melt(genes_expression, id = "ID", variable.name = "sample", value.name = "expression")
 genes_expression[, individual := substr(sample, 1,6)]
 genes_expression[, condition := substr(sample, 8,8)]
 genes_expression[, population := substr(sample, 1,3)]
 genes_expression[, expression_inrt := inverseNormalRankTransform(expression),by=.(ID,condition)]
 names(genes_expression)[1] = "gene"

intronic_expression = fread(paste(EVO_IMMUNO_POP, "Maxime/Evo_Immuno_pop_data/05_IntronicReadCounts_HTSeq/AllSamples_Intronic_Gene_count.txt", sep=""))
intronic_expression[, gene := substr(intron_id, 1,nchar(intron_id)-2)]
intronic_expression = intronic_expression[gene %in% genes_expression$gene]
genes_expression = genes_expression[gene %in% intronic_expression$gene]
intronic_expression = melt(intronic_expression, id = "gene", variable.name = "sample", value.name = "transcription_rate")
intronic_expression = intronic_expression[grepl("RPKM", sample)]
intronic_expression[, sample := as.character(sample)]
intronic_expression[, individual := substr(sample, 1,6)]
intronic_expression[, condition := substr(sample, 8,8)]
intronic_expression[, population := substr(sample, 1,3)]
intronic_expression[, transcription_rate:=as.numeric(transcription_rate)]
intronic_expression[, transcription_rate_inrt := inverseNormalRankTransform(transcription_rate),by=.(gene,condition)]
intronic_expression[, gene := paste(gene, "intron", sep="_")]

genes_names = genes_expression[,unique(gene)]
genes_introns_names =intronic_expression[,unique(gene)]

geneAnnot=fread(sprintf("%s/Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/GeneAnnot_expressed.txt",EVO_IMMUNO_POP))
miRAnnot=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_TSS_annot_DeRie2017.txt',EVO_IMMUNO_POP),sep='\t')

######### SANITY TEST ############
my_miR='hsa_miR_146b_3p'
my_gene='ENSG00000184557' #(SOCS3)
my_cond=4

my_data=merge(genes_expression[gene== my_gene & condition== my_cond,],miRNA_expression[miRNA== gsub('_','-', my_miR) & condition== my_cond,],suffix=c('.gene','.miRNA'),by=c('individual','population','condition'))
cor(my_data$expression_inrt.gene,my_data$expression_inrt.miRNA)^2

my_data=merge(intronic_expression[gene== paste(my_gene,'intron',sep='_') & condition== my_cond,],miRNA_expression[miRNA== gsub('_','-', my_miR) & condition== my_cond,],suffix=c('.gene','.miRNA'),by=c('individual','population','condition'))
cor(my_data$transcription_rate_inrt,my_data$expression_inrt)^2

##################################


##################################################################################################
##Let's create a great data.table with one line per individual miRNA expression, gene_expression##
##Population and intronic expression                                                            ##
##################################################################################################

for(cond_to_investigate in 1:5){

cat(cond_to_investigate,':\n')
##Compute sva
intronic_expression_cond = intronic_expression[condition == cond_to_investigate]

intronic_expression_matrix = dcast(intronic_expression_cond, gene~individual, value.var = "transcription_rate_inrt", identity, fill = NA)
intron_list = intronic_expression_matrix$gene
intronic_expression_matrix[, gene := NULL]
intronic_expression_matrix = as.matrix(intronic_expression_matrix)
#intronic_expression_matrix=matrix(as.numeric(intronic_expression_matrix),nrow(intronic_expression_matrix))
rownames(intronic_expression_matrix) = intron_list

genes_expression_cond = genes_expression[condition == cond_to_investigate]

gene_expression_matrix = dcast(genes_expression_cond, gene~individual, value.var = "expression_inrt", identity, fill = NA)
gene_list = gene_expression_matrix$gene
gene_expression_matrix[, gene := NULL]
gene_expression_matrix = as.matrix(gene_expression_matrix)
rownames(gene_expression_matrix) = gene_list

miRNA_expression_cond = miRNA_expression[condition == cond_to_investigate]

miRNA_expression_matrix = dcast(miRNA_expression_cond, miRNA~individual, value.var = "expression_inrt", identity, fill = NA)
miR_list= miRNA_expression_matrix$miRNA
miRNA_expression_matrix[, miRNA := NULL]
miRNA_expression_matrix = as.matrix(miRNA_expression_matrix)
rownames(miRNA_expression_matrix) = miR_list

ii=intersect(colnames(gene_expression_matrix),colnames(intronic_expression_matrix))
iii=intersect(ii,colnames(miRNA_expression_matrix))


#cor(miRNA_expression_matrix[gsub('_','-',my_miR),iii], gene_expression_matrix[my_gene,iii])
#cor(miRNA_expression_matrix[gsub('_','-',my_miR),iii], intronic_expression_matrix[paste(my_gene,'intron',sep='_'),iii])

indiv_in_genes = data.table(individual = colnames(gene_expression_matrix))
indiv_in_genes[, population := substr(individual, 1,3)]


mod = model.matrix(~as.factor(population), data=indiv_in_genes)
sva_matrix=sva(gene_expression_matrix, mod = mod)$sv
rownames(sva_matrix) = colnames(gene_expression_matrix)
colnames(sva_matrix) = paste('SV',1:ncol(sva_matrix),sep='')

sva_results = as.data.table(sva_matrix)
sva_col_names = colnames(sva_results[1])
sva_results[,individual := indiv_in_genes$individual]
 

RPKM_Mat=SlicedData$new()
RPKM_Mat$CreateFromMatrix(gene_expression_matrix[,match(iii,colnames(gene_expression_matrix))])
RPKM_Mat$ResliceCombined(6000)

miR_Mat=SlicedData$new()
miR_Mat$CreateFromMatrix(miRNA_expression_matrix[,match(iii,colnames(miRNA_expression_matrix))])
miR_Mat$ResliceCombined(6000)

miR_pos = miRAnnot[,mget(c('hsa_ID','mat_chr','mat_start'))]
miR_pos[,mat_chr:=gsub('chr','', mat_chr)]
miR_pos=as.data.frame(miR_pos)

gene_pos=as.data.frame(geneAnnot)[c(1,2,4,5)]
colnames(gene_pos)= c('gene','chr','start','end') 

Covar_Mat=SlicedData$new()
Covar_Mat$CreateFromMatrix(t(cbind(sva_matrix,pop=ifelse(indiv_in_genes$pop=='AFB',0,1)))[,match(iii,rownames(sva_matrix))])
Covar_Mat$ResliceCombined(6000)

res_miRNA_gene_SV = Matrix_eQTL_main(miR_Mat, RPKM_Mat, cvrt = Covar_Mat, 
      output_file_name = NULL, pvOutputThreshold = 0.05,
       output_file_name.cis = NULL, 
       pvOutputThreshold.cis = 0.05,
       snpspos = miR_pos, 
       genepos = gene_pos,
       cisDist = 1e3)


res_miRNA_gene_SV$trans$eqtls$r = res_miRNA_gene_SV$trans$eqtls$statistic/sqrt(res_miRNA_gene_SV$trans$eqtls$statistic^2+res_miRNA_gene_SV$param$dfFull)
cor_miRNA_gene_SV=data.table(res_miRNA_gene_SV$trans$eqtls)


Covar_Mat=SlicedData$new()
Covar_Mat$CreateFromMatrix(matrix(ifelse(indiv_in_genes$pop=='AFB',0,1),nrow=1)[,match(iii,rownames(sva_matrix)),drop=F])
Covar_Mat$ResliceCombined(6000)

res_miRNA_gene_noSV = Matrix_eQTL_main(miR_Mat, RPKM_Mat, cvrt = Covar_Mat, 
       output_file_name = NULL, pvOutputThreshold = 0.05,
        output_file_name.cis = NULL, 
        pvOutputThreshold.cis = 0.05,
        snpspos = miR_pos, 
        genepos = gene_pos,
        cisDist = 1e3)


res_miRNA_gene_noSV$trans$eqtls$r = res_miRNA_gene_noSV$trans$eqtls$statistic/sqrt(res_miRNA_gene_noSV$trans$eqtls$statistic^2+res_miRNA_gene_noSV$param$dfFull)
cor_miRNA_gene_noSV=data.table(res_miRNA_gene_noSV$trans$eqtls)

annotCols=c('Ensembl.Gene.ID','Chromosome.Name','Gene Start (bp)','Gene End (bp)','Associated Gene Name','NS_mean')
cor_miRNA_gene_SV = merge(cor_miRNA_gene_SV, geneAnnot[,mget(annotCols)],by.x='gene',by.y='Ensembl.Gene.ID')
cor_miRNA_gene_SV=cor_miRNA_gene_SV[order(r)]

cor_miRNA_gene_noSV = merge(cor_miRNA_gene_noSV, geneAnnot[,mget(annotCols)],by.x='gene',by.y='Ensembl.Gene.ID')
cor_miRNA_gene_noSV=cor_miRNA_gene_noSV[order(r)]

fwrite(cor_miRNA_gene_SV,file=sprintf('%s/Maxime/miRNA_V2/data/13.correlation_transcription_miRNA/cor_miRNA_gene_MatrixEQTL_%s_SVadj.txt',EVO_IMMUNO_POP,cond_to_investigate),sep='\t')
fwrite(cor_miRNA_gene_noSV,file=sprintf('%s/Maxime/miRNA_V2/data/13.correlation_transcription_miRNA/cor_miRNA_gene_MatrixEQTL_%s_noSV.txt',EVO_IMMUNO_POP,cond_to_investigate),sep='\t')



RPKM_Mat=SlicedData$new()
RPKM_Mat$CreateFromMatrix(intronic_expression_matrix[,match(iii,colnames(intronic_expression_matrix))])
RPKM_Mat$ResliceCombined(6000)

gene_pos=as.data.frame(geneAnnot)[c(1,2,4,5)]
colnames(gene_pos)= c('gene','chr','start','end') 
gene_pos$gene=paste(gene_pos$gene,'intron',sep='_')

Covar_Mat=SlicedData$new()
Covar_Mat$CreateFromMatrix(t(cbind(sva_matrix,pop=ifelse(indiv_in_genes$pop=='AFB',0,1)))[,match(iii,rownames(sva_matrix))])
Covar_Mat$ResliceCombined(6000)

res_miRNA_intron_SV = Matrix_eQTL_main(miR_Mat, RPKM_Mat, cvrt = Covar_Mat, 
       output_file_name = NULL, pvOutputThreshold = 0.05,
        output_file_name.cis = NULL, 
        pvOutputThreshold.cis = 0.05,
        snpspos = miR_pos, 
        genepos = gene_pos,
        cisDist = 1e3)

res_miRNA_intron_SV$trans$eqtls$r = res_miRNA_intron_SV$trans$eqtls$statistic/sqrt(res_miRNA_intron_SV$trans$eqtls$statistic^2+res_miRNA_intron_SV$param$dfFull)
cor_miRNA_intron_SV=data.table(res_miRNA_intron_SV$trans$eqtls)


Covar_Mat=SlicedData$new()
Covar_Mat$CreateFromMatrix(matrix(ifelse(indiv_in_genes$pop=='AFB',0,1),nrow=1)[,match(iii,rownames(sva_matrix)),drop=F])
Covar_Mat$ResliceCombined(6000)


res_miRNA_intron_noSV = Matrix_eQTL_main(miR_Mat, RPKM_Mat, cvrt = Covar_Mat, 
        output_file_name = NULL, pvOutputThreshold = 0.05,
        output_file_name.cis = NULL,
        pvOutputThreshold.cis = 0.05,
        snpspos = miR_pos, 
        genepos = gene_pos,
        cisDist = 1e3)

res_miRNA_intron_noSV$trans$eqtls$r = res_miRNA_intron_noSV$trans$eqtls$statistic/sqrt(res_miRNA_intron_noSV$trans$eqtls$statistic^2+res_miRNA_intron_noSV$param$dfFull)
cor_miRNA_intron_noSV=data.table(res_miRNA_intron_noSV$trans$eqtls)

geneAnnot[,intron_name:=paste(Ensembl.Gene.ID,'intron',sep='_')]
annotCols=c('intron_name','Chromosome.Name','Gene Start (bp)','Gene End (bp)','Associated Gene Name','NS_mean')
cor_miRNA_intron_SV = merge(cor_miRNA_intron_SV, geneAnnot[,mget(annotCols)],by.x='gene',by.y='intron_name')
cor_miRNA_intron_SV = cor_miRNA_intron_SV[order(-abs(r))]

cor_miRNA_intron_noSV = merge(cor_miRNA_intron_noSV, geneAnnot[,mget(annotCols)],by.x='gene',by.y='intron_name')
cor_miRNA_intron_noSV=cor_miRNA_intron_noSV[order(-abs(r))]

fwrite(cor_miRNA_intron_SV,file=sprintf('%s/Maxime/miRNA_V2/data/13.correlation_transcription_miRNA/cor_miRNA_intron_MatrixEQTL_%s_SVadj.txt',EVO_IMMUNO_POP,cond_to_investigate),sep='\t')
fwrite(cor_miRNA_intron_noSV,file=sprintf('%s/Maxime/miRNA_V2/data/13.correlation_transcription_miRNA/cor_miRNA_intron_MatrixEQTL_%s_noSV.txt',EVO_IMMUNO_POP,cond_to_investigate),sep='\t')
}
# all_individuals = unique(c(miRNA_expression$individual, genes_expression$individual))
# data_to_analyse = data.table(individual = all_individuals)

# data_to_analyse[, population := ifelse(substr(individual,1,3)== "AFB", 1, 0)]
# data_to_analyse[, miRNA_expression := miRNA_expression[match(data_to_analyse$individual, individual), expression]]

# for (cn in sva_col_names){
#   data_to_analyse[ ,eval(cn) := sva_results[match(data_to_analyse$individual, individual), get(eval(cn))]]
# }
# 
# formula = paste("gene_expression ~ miRNA_expression + gene_transcription + population", paste(sva_col_names, collapse = " + "), sep=" + ")
