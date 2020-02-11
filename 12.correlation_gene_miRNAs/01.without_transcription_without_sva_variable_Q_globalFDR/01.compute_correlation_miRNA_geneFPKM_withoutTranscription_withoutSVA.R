####################################################################################################
##Aim : compute correlation between miRNA expression and gene expression adjusted on transcription##
####################################################################################################
require(Hmisc)

###############################
##WARNING : SCRIPT NOT TESTED##
###############################


#############
##Arguments##
#############
args = commandArgs(trailingOnly=TRUE)
cond_to_investigate = as.numeric(args[1])
miRNA_number_to_investigate = as.numeric(args[2])
#############
##Load data##
#############

miRNA_expression= fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""))
names(miRNA_expression)[1] = "miRNA"
miRNA_expression = melt(miRNA_expression, id = "miRNA", variable.name = "sample", value.name = "expression")
miRNA_expression[, individual := substr(sample, 1,6)]
miRNA_expression[, condition := substr(sample, 8,8)]
miRNA_expression[, population := substr(sample, 1,3)]
miRNA_studied = sort(unique(miRNA_expression$miRNA))
miRNA_studied = miRNA_studied[miRNA_number_to_investigate]
miRNA_expression = miRNA_expression[condition == cond_to_investigate]
miRNA_expression = miRNA_expression[miRNA == miRNA_studied]

genes_expression = fread(paste(EVO_IMMUNO_POP, "Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/FPKM_matrix.txt", sep=""))
genes_expression = melt(genes_expression, id = "ID", variable.name = "sample", value.name = "expression")
genes_expression[, individual := substr(sample, 1,6)]
genes_expression[, condition := substr(sample, 8,8)]
genes_expression[, population := substr(sample, 1,3)]
genes_expression = genes_expression[condition == cond_to_investigate]
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
intronic_expression = intronic_expression[condition == cond_to_investigate]
intronic_expression[, gene := paste(gene, "intron", sep="_")]

genes_names = genes_expression[,unique(gene)]
genes_introns_names =intronic_expression[,unique(gene)]
##################################################################################################
##Let's create a great data.table with one line per individual miRNA expression, gene_expression##
##Population and intronic expression                                                            ##
##################################################################################################
all_individuals = unique(c(miRNA_expression$individual, genes_expression$individual))
data_to_analyse = data.table(individual = all_individuals)
gene_matrix = dcast(genes_expression, individual ~ gene, value.var = "expression", identity, fill = NA)
intronic_expression_matrix = dcast(intronic_expression, individual ~ gene, value.var = "transcription_rate", identity, fill = NA)

data_to_analyse[, population := ifelse(substr(individual,1,3)== "AFB", 1, 0)]
data_to_analyse[, miRNA_expression := miRNA_expression[match(data_to_analyse$individual, individual), expression]]


formula = "gene_expression ~ miRNA_expression + population"

compute_correlation <- function(gene_name){
  print(gene_name)
  data_to_analyse[, gene_expression := gene_matrix[match(data_to_analyse$individual, individual), get(eval(gene_name))]]
  data_to_analyse[, gene_transcription := as.numeric(intronic_expression_matrix[match(data_to_analyse$individual, individual), get(paste(eval(gene_name), "intron",sep="_"))])]
  mod = glm(data_to_analyse,formula =  formula, family = "gaussian")
  tvalue = summary(mod)$coefficients["miRNA_expression", "t value"]
  df = mod$df.residual
  r = tvalue/sqrt(tvalue**2 + df)
  data_to_analyse[, gene_expression := NULL]
  data_to_analyse[, gene_transcription := NULL]
  return(data.table(gene_name, miRNA = miRNA_studied, condition = cond_to_investigate, r , pvalue = summary(mod)$coefficients["miRNA_expression", "Pr(>|t|)"]))
}

results = rbindlist(lapply(genes_names, compute_correlation))

write.table(results, paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/correlation_miRNA_geneRPKM_transcriptionAdjusted_withoutSVA/correlation_miRNA_geneRPKM_transcriptionAdjusted_withoutSVA_", miRNA_studied, "_", cond_to_investigate, ".tsv", sep=""), sep="\t", quote = F, row.names = F)
