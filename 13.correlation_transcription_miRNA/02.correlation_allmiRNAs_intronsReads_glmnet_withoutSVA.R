#########################################################################################################
##Aim : compute correlation between all miRNAs expression and gene expression adjusted on transcription##
#########################################################################################################
require(Hmisc)
require(glmnet)
require(stabs)
#############
##Arguments##
#############
args = commandArgs(trailingOnly=TRUE)
cond_to_investigate = as.numeric(args[1])
group_of_genes = as.numeric(args[2])
#Functions to normalize everything#
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
print("test")
variance_intronic_expression = intronic_expression[, .(vari = var(transcription_rate, na.rm = T)), by = c("gene")]
variance_intronic_expression[, genes_bis := substr(gene, 1,15)]
genes_to_remove = variance_intronic_expression[ vari == 0, gene]
genes_names = variance_intronic_expression[ vari != 0, genes_bis]
intronic_expression = intronic_expression[!(gene %in% genes_to_remove)]
print("pouloup")

##################################################################################################
##Let's create a great data.table with one line per individual miRNA expression, gene_expression##
##Population and intronic expression                                                            ##
##################################################################################################


all_individuals = unique(c(miRNA_expression$individual, intronic_expression$individual))
data_to_analyse = data.table(individual = all_individuals)
gene_matrix = dcast(genes_expression, individual ~ gene, value.var = "expression", identity, fill = NA)
intronic_expression_matrix = dcast(intronic_expression, individual ~ gene, value.var = "transcription_rate", identity, fill = NA)

data_to_analyse[, population := ifelse(substr(individual,1,3)== "AFB", 1, 0)]
data_to_analyse[, miRNA_expression := miRNA_expression[match(data_to_analyse$individual, individual), expression]]



individual_to_keep = intersect(gene_matrix$individual, intersect(data_to_analyse$individual, miRNA_expression_matrix$individual))
data_to_analyse = data_to_analyse[individual %in% individual_to_keep]
data_to_analyse = data_to_analyse[order(individual)]
miRNA_expression_matrix = miRNA_expression_matrix[match(data_to_analyse$individual, individual)]
miRNA_expression_matrix[, individual := NULL]

data_to_analyse = cbind(data_to_analyse, miRNA_expression_matrix)

glmnet.WithPenalty=function (x, y, q, type = c("conservative", "anticonservative"), penalty.factor,...){
    glmnet.lasso(x, y, q=q+sum(1-penalty.factor), type, penalty.factor=penalty.factor,...)
    }

run_glmnet <- function(gene_name, PFER = 1, Q = 30, alpha = 0.5){
  print(gene_name)
  temp = data_to_analyse[1:.N]
  gene_transcription = as.numeric(intronic_expression_matrix[match(temp$individual, individual), get(paste(eval(gene_name), "intron",sep="_"))])
  temp[, individual := NULL]
  pen_factor = rep(1, ncol(temp))
  names(pen_factor) = names(temp)
  pen_factor["population"] = 0
  # pen_factor["gene_transcription"] = 0

  stab.glmnet <- stabsel(x = as.matrix(temp), y = gene_transcription,
              fitfun = glmnet.WithPenalty, PFER=PFER, q=Q,
              args.fitfun=list(penalty.factor = pen_factor, alpha=alpha)
              )

  res = data.table(selected_variables = names(temp)[stab.glmnet$selected], probs_selection = stab.glmnet$max[stab.glmnet$selected])
  res[, gene := gene_name]
  res[, condition := cond_to_investigate]
  setcolorder(res, c("gene", "condition", "selected_variables"))
  return(res)
}

results = rbindlist(lapply(genes_names[1:length(genes_names) %% 20 == group_of_genes], run_glmnet))
write.table(results, paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/13.correlation_transcription_miRNA/glmnet_miRNA_intronReads_withoutSVA/gene_transcription_glmnet_alpha0.5_cond", cond_to_investigate, "_", group_of_genes, ".tsv", sep=""), sep="\t", quote=F , row.names = F)
print("done")
