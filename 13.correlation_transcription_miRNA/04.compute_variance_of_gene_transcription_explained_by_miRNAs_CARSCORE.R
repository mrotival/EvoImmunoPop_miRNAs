require(care)
require(sva)
args = commandArgs(trailingOnly=TRUE)
cond_to_investigate = as.numeric(args[1])
##############
##Load datas##
##############
begin = proc.time()

inrt=function(x){ qnorm(rank(x,ties='random',na.last = NA)/(length(x)+1),mean(x,na.rm=T),sd(x,na.rm=T))}
#as.data.table(t(apply(corrected_count_transformed, 2, inrt)))
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


miRNA_expression = miRNA_expression[condition == cond_to_investigate]
miRNA_expression = miRNA_expression[, expression := inrt(expression), by = miRNA]
miRNA_expression_matrix = dcast(miRNA_expression, formula = individual ~ miRNA, identity, value.var = "expression", fill = NA)


genes_expression = fread(paste(EVO_IMMUNO_POP, "Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/FPKM_matrix.txt", sep=""))
genes_expression = melt(genes_expression, id = "ID", variable.name = "sample", value.name = "expression")
genes_expression[, individual := substr(sample, 1,6)]
genes_expression[, condition := substr(sample, 8,8)]
genes_expression[, population := substr(sample, 1,3)]
genes_expression = genes_expression[condition == cond_to_investigate]
names(genes_expression)[1] = "gene"
genes_expression[, expression := inrt(expression), by = gene]


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
intronic_expression[, transcription_rate := as.numeric(transcription_rate)]
intronic_expression[, transcription_rate := inrt(transcription_rate), by = gene]

#We need to have variation within intronic expression, otherwise we can do nothing
variance_intronic_expression = intronic_expression[, .(vari = var(transcription_rate, na.rm = T)), by = c("gene")]
variance_intronic_expression[, genes_bis := substr(gene, 1,15)]
genes_to_remove = variance_intronic_expression[ vari == 0, gene]
genes_names = variance_intronic_expression[ vari != 0, genes_bis]
intronic_expression = intronic_expression[!(gene %in% genes_to_remove)]

# genes_names = genes_expression[,unique(gene)]
genes_introns_names =intronic_expression[,unique(gene)]

glmnet_data = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/glmnet_correlation_miRNA_geneRPKM_transcriptionAdjusted_svaAdjusted/gene_expression_glmnet_alpha0.5_allConds.tsv", sep=""))
glmnet_data = glmnet_data[selected_variables %in% miRNA_expression$miRNA]
glmnet_data = glmnet_data[condition ==  cond_to_investigate]


##################################################################################################
##Let's create a great data.table with one line per individual miRNA expression, gene_expression##
##Population and intronic expression                                                            ##
##################################################################################################
##Compute sva
all_individuals = unique(c(miRNA_expression$individual, intronic_expression$individual))
data_to_analyse = data.table(individual = all_individuals)
gene_matrix = dcast(genes_expression, individual ~ gene, value.var = "expression", identity, fill = NA)
intronic_expression_matrix = dcast(intronic_expression, individual ~ gene, value.var = "transcription_rate", identity, fill = NA)

data_to_analyse[, population := ifelse(substr(individual,1,3)== "AFB", 1, 0)]



individual_to_keep = intersect(gene_matrix$individual, intersect(data_to_analyse$individual, miRNA_expression_matrix$individual))
data_to_analyse = data_to_analyse[individual %in% individual_to_keep]
data_to_analyse = data_to_analyse[order(individual)]
miRNA_matrix = miRNA_expression_matrix[match(data_to_analyse$individual, individual)]
miRNA_matrix[, individual := NULL]

# data_to_analyse = cbind(data_to_analyse, miRNA_expression_matrix)


compute_CAR_score <- function(gene_name){
  print(paste(gene_name, cond_to_investigate))
  data_to_analyse_t = data_to_analyse[1:.N]
  # data_to_analyse_t[, gene_expression := gene_matrix[match(data_to_analyse_t$individual, individual), get(eval(gene_name))]]
  data_to_analyse_t[, gene_transcription := as.numeric(intronic_expression_matrix[match(data_to_analyse_t$individual, individual), get(paste(eval(gene_name), "intron",sep="_"))])]
  miRNA_to_include = glmnet_data[gene == gene_name, unique(selected_variables)]
  miRNA_to_include_new = c()
  for (mi in miRNA_to_include){
    new_name = gsub("-", "_", mi)
    miRNA_to_include_new = c(miRNA_to_include_new, new_name)
    data_to_analyse_t[ ,eval(new_name) := miRNA_matrix[match(data_to_analyse$individual, individual), get(eval(mi))]]
  }
  data_to_analyse_t = data_to_analyse_t[complete.cases(data_to_analyse_t)]
  gene_expression_vector = data_to_analyse_t[, gene_transcription]
  data_to_analyse_t[, gene_transcription := NULL]
  data_to_analyse_t[, individual := NULL]

  car = carscore( data_to_analyse_t, gene_expression_vector, lambda = 0)
  variance_explained = car**2
  variance_explained_miRNAs = sort(variance_explained[miRNA_to_include_new], decreasing = T)
  carscores_miRNAs = sort(car[miRNA_to_include_new], decreasing = T)
  print(carscores_miRNAs)
  # variance_explained_SVA = variance_explained[sva_col_names]
  variance_explained_population = variance_explained["population"]
  # variance_explained_transcription = variance_explained["gene_transcription"]

  miRNA_max_variance_explained = names(variance_explained_miRNAs)[1]
  variance_explained_by_max_miRNA = variance_explained_miRNAs[1]
  res = data.table(gene = gene_name, condition = cond_to_investigate,
                    variance_explained_by_pop = variance_explained_population,
                   # variance_explained_by_transcription = variance_explained_transcription,
                    variance_explained_by_miRNAs = sum(variance_explained_miRNAs),
                    # variance_explained_by_SVAs = sum(variance_explained_SVA),
                    number_miRNA_in_model = length(miRNA_to_include_new))

  miRNA_variance = as.data.table(t(data.table(variance_explained_miRNAs[1:10])))
  names(miRNA_variance) = paste("miRNA_variance", 1:10, sep = "_")
  miRNAs = as.data.table(t(data.table(names(variance_explained_miRNAs)[1:10])))
  names(miRNAs) = paste("miRNA", 1:10, sep = "_")
  carscore_sign = ifelse(carscores_miRNAs<0, -1, 1)
  carscore_sign = as.data.table(t(data.table(carscore_sign[1:10])))
  names(carscore_sign) = paste("miRNA_carscore_sign", 1:10, sep="_")
  res = cbind(res, miRNA_variance, miRNAs, carscore_sign)
  return(res)
}
results = rbindlist(lapply(genes_names, compute_CAR_score))
write.table(results, file = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/13.correlation_transcription_miRNA/variance_explained_by_miRNA_CARScore/",
                                    "variance_explained_by_miRNA_CARScore_", cond_to_investigate, ".tsv", sep=""),
                                    sep="\t", quote = F, row.names = F)
