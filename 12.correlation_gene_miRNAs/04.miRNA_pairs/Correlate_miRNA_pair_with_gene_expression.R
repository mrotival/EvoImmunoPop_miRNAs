
DATA_DIR=sprintf("%s/Maxime/miRNA_V2/data/",EVO_IMMUNO_POP)
FIG_DIR=sprintf("%s/Maxime/miRNA_V2/figures/Revisions",EVO_IMMUNO_POP)

library(rtracklayer)
######################################################################################################################
##Aim : compute correlation between expression of a pair of miRNAs with colocalized targets and gene expression		##
##								 adjusted on transcription															##
######################################################################################################################
require(Hmisc)
require(sva)

#############
##Arguments##
#############
args = commandArgs(trailingOnly=TRUE)
cond_to_investigate = as.numeric(args[1])
miRNA_pair_to_investigate = as.numeric(args[2])

#############
##Load data##
#############

# mRNA expression
genes_expression = fread(paste(EVO_IMMUNO_POP, "Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/FPKM_matrix.txt", sep=""))
genes_expression = melt(genes_expression, id = "ID", variable.name = "sample", value.name = "expression")
genes_expression[, individual := substr(sample, 1,6)]
genes_expression[, condition := substr(sample, 8,8)]
genes_expression[, population := substr(sample, 1,3)]
genes_expression = genes_expression[condition == cond_to_investigate]
names(genes_expression)[1] = "gene"

# intronic expression
# intronic_expression = fread(paste(EVO_IMMUNO_POP, "Maxime/Evo_Immuno_pop_data/05_IntronicReadCounts_HTSeq/AllSamples_Intronic_Gene_count.txt", sep=""))
# intronic_expression[, gene := substr(intron_id, 1,nchar(intron_id)-2)]
# intronic_expression = intronic_expression[gene %in% genes_expression$gene]
# genes_expression = genes_expression[gene %in% intronic_expression$gene]
# intronic_expression = melt(intronic_expression, id = "gene", variable.name = "sample", value.name = "transcription_rate")
# intronic_expression = intronic_expression[grepl("RPKM", sample)]
# intronic_expression[, sample := as.character(sample)]
# intronic_expression[, individual := substr(sample, 1,6)]
# intronic_expression[, condition := substr(sample, 8,8)]
# intronic_expression[, population := substr(sample, 1,3)]
# intronic_expression = intronic_expression[condition == cond_to_investigate]
# intronic_expression[, gene := paste(gene, "intron", sep="_")]

genes_names = genes_expression[,unique(gene)]
# genes_introns_names =intronic_expression[,unique(gene)]

# pairs of colocalized miRNAs with >50 common targets expressed in monocytes
miR_pair_gene=fread(sprintf('%s/19_miR_targets/mir_pairs_colocalization/miR_pairs_genecount.txt',DATA_DIR))
miR_pair_gene=miR_pair_gene[Padj_Coloc_enrich<.01 & NbColoc_gene>100,]

miR_pair=fread(sprintf('%s/19_miR_targets/mir_pairs_colocalization/miR_pairs.txt',DATA_DIR))
miR_pair=miR_pair[paste(miRNA.x,miRNA.y)%chin%paste(miR_pair_gene$miRNA.x,miR_pair_gene$miRNA.y),]
miR_pair=miR_pair[EnsemblGeneID %chin% genes_names,]

miR_pair_gene=miR_pair[, .(NbColoc_gene=sum(NbColoc>0)),by=.(miRNA.x,miRNA.y)]
miR_pair_gene=miR_pair_gene[NbColoc_gene>100,]
miR_pair=miR_pair[paste(miRNA.x,miRNA.y)%chin%paste(miR_pair_gene$miRNA.x,miR_pair_gene$miRNA.y),]
nrow(miR_pair_gene) # 390

# miRNAs
miRNA_expression= fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""))
names(miRNA_expression)[1] = "miRNA"
miRNA_expression = melt(miRNA_expression, id = "miRNA", variable.name = "sample", value.name = "expression")
miRNA_expression[, individual := substr(sample, 1,6)]
miRNA_expression[, condition := substr(sample, 8,8)]
miRNA_expression[, population := substr(sample, 1,3)]
miRNA_studied_1 = miR_pair_gene[miRNA_pair_to_investigate,miRNA.x]
miRNA_studied_2 = miR_pair_gene[miRNA_pair_to_investigate,miRNA.y]
miRNA_expression = miRNA_expression[condition == cond_to_investigate,]
miR1_expression = miRNA_expression[miRNA == miRNA_studied_1,]
miR2_expression = miRNA_expression[miRNA == miRNA_studied_2,]

##################################################################################################
##Let's create a great data.table with one line per individual, expression of both miRNAs,		##
## gene_expression, Population and intronic expression                                          ##	
##################################################################################################
##Compute sva
indiv_in_genes = data.table(individual = unique(genes_expression$individual))
indiv_in_genes[, population := substr(individual, 1,3)]
gene_matrix = dcast(genes_expression, gene ~ individual, value.var = "expression", identity, fill = NA)
gene_list = gene_matrix$gene
gene_matrix[, gene := NULL]
gene_matrix = as.matrix(gene_matrix)
rownames(gene_matrix) = gene_list
# mod = model.matrix(~as.factor(population), data=indiv_in_genes)
# sva_results = as.data.table(sva(gene_matrix, mod = mod)$sv)
# sva_col_names = colnames(sva_results[1])
# sva_results[,individual := colnames(gene_matrix)]

all_individuals = unique(c(miRNA_expression$individual, genes_expression$individual))
data_to_analyse = data.table(individual = all_individuals)
gene_matrix = dcast(genes_expression, individual ~ gene, value.var = "expression", identity, fill = NA)
# intronic_expression_matrix = dcast(intronic_expression, individual ~ gene, value.var = "transcription_rate", identity, fill = NA)

data_to_analyse[, population := ifelse(substr(individual,1,3)== "AFB", 1, 0)]
x=miR1_expression[match(data_to_analyse$individual, individual), expression]
data_to_analyse[, miR1_expression := x]
x=miR2_expression[match(data_to_analyse$individual, individual), expression]
data_to_analyse[, miR2_expression := x]
# 
# for (cn in sva_col_names){
#   data_to_analyse[ ,eval(cn) := sva_results[match(data_to_analyse$individual, individual), get(eval(cn))]]
# }

# formula = paste("gene_expression ~ miR1_expression + miR2_expression + miR1_expression:miR2_expression + gene_transcription + population", paste(sva_col_names, collapse = " + "), sep=" + ")
formula = "gene_expression ~ miR1_expression + miR2_expression + miR1_expression:miR2_expression + population"

compute_jointEffect <- function(gene_name){
  require(car)
  require(care)
  print(gene_name)
  data_to_analyse[, gene_expression := gene_matrix[match(data_to_analyse$individual, individual), get(eval(gene_name))]]
#   data_to_analyse[, gene_transcription := as.numeric(intronic_expression_matrix[match(data_to_analyse$individual, individual), get(paste(eval(gene_name), "intron",sep="_"))])]
  data_to_analyse[, interactionTerm:=miR1_expression*miR2_expression]
  data_to_analyse=data_to_analyse[!is.na(gene_expression) & !is.na(interactionTerm),]
  mod = lm(data_to_analyse,formula =  formula)
  stats=linearHypothesis(mod,c('miR1_expression=0','miR2_expression=0','miR1_expression:miR2_expression=0'))
  Pval=stats[['Pr(>F)']][2]
  Pval_int=summary(mod)$coeff['miR1_expression:miR2_expression',4]
#  car_w=carscore(data_to_analyse[,gene_expression],data_to_analyse[,mget(c('miR1_expression','miR2_expression','interactionTerm','gene_transcription','population',sva_col_names))])
  car_w=carscore(data_to_analyse[,gene_expression],data_to_analyse[,mget(c('miR1_expression','miR2_expression','interactionTerm','population'))])
  data_to_analyse[, gene_expression := NULL]
#   data_to_analyse[, gene_transcription := NULL]
  return(data.table(gene=gene_name, miR1 = miRNA_studied_1, miR2 = miRNA_studied_2, condition = cond_to_investigate, miR1_car=car_w[1],miR2_car=car_w[2], miR12_car=car_w[3], pvalue = Pval,Pval_int=Pval_int))
}

tim=Sys.time()
results = rbindlist(lapply(genes_names, compute_jointEffect))
print(Sys.time()-tim)

fwrite(results, file=paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/19_miR_targets/mir_pairs_colocalization/CorrelationGeneExpression/JointEffect_geneRPKM_miR1_", miRNA_studied_1, "_miR2_",miRNA_studied_2,"_cond", cond_to_investigate, "_withInt.tsv", sep=""), sep="\t")




