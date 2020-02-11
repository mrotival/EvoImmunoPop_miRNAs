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
perm = as.numeric(args[2])
group_of_genes = as.numeric(args[3])
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
# intronic_expression[, transcription_rate := as.numeric(transcription_rate)]
# intronic_expression[, transcription_rate := inrt(transcription_rate), by = gene]


genes_names = genes_expression[,unique(gene)]
##################################################################################################
##Let's create a great data.table with one line per individual miRNA expression, gene_expression##
##Population and intronic expression                                                            ##
##################################################################################################


all_individuals = unique(c(miRNA_expression$individual, genes_expression$individual))
data_to_analyse = data.table(individual = all_individuals)
gene_matrix = dcast(genes_expression, individual ~ gene, value.var = "expression", identity, fill = NA)

data_to_analyse[, population := ifelse(substr(individual,1,3)== "AFB", 1, 0)]
# data_to_analyse[, miRNA_expression := miRNA_expression[match(data_to_analyse$individual, individual), expression]]



individual_to_keep = intersect(gene_matrix$individual, intersect(data_to_analyse$individual, miRNA_expression_matrix$individual))
data_to_analyse = data_to_analyse[individual %in% individual_to_keep]
data_to_analyse = data_to_analyse[order(individual)]


if(perm==0){
	miRNA_expression_matrix = miRNA_expression_matrix[match(data_to_analyse$individual, individual)]
}else{
	set.seed(perm)
	miRNA_expression_matrix = miRNA_expression_matrix[match(sample(data_to_analyse$individual,replace=F), individual)]
}

miRNA_expression_matrix[, individual := NULL]

data_to_analyse = cbind(data_to_analyse, miRNA_expression_matrix)

Qlist=c(3*1:20)

glmnet.WithPenalty=function (x, y, q, type = c("conservative", "anticonservative"), penalty.factor,...){
    glmnet.lasso(x, y, q=q+sum(1-penalty.factor), type, penalty.factor=penalty.factor,...)
    }

run_glmnet <- function(gene_name, PFER = 1, Q = 30, alpha = 0.5){
  print(gene_name)
  temp = data_to_analyse[1:.N]
  gene_expression = gene_matrix[match(data_to_analyse$individual, individual), get(eval(gene_name))]
  # temp[, gene_transcription := as.numeric(intronic_expression_matrix[match(temp$individual, individual), get(paste(eval(gene_name), "intron",sep="_"))])]
  temp[, individual := NULL]
  pen_factor = rep(1, ncol(temp))
  names(pen_factor) = names(temp)
  pen_factor["population"] = 0
  # pen_factor["gene_transcription"] = 0

  res=list()
  tim0=Sys.time()

for(Q in Qlist){
	cat(Q,'')
  tim=Sys.time()

  stab.glmnet <- stabsel(x = as.matrix(temp), y = gene_expression,
              fitfun = glmnet.WithPenalty, PFER=PFER, q=Q,
              args.fitfun=list(penalty.factor = pen_factor, alpha=alpha)
              )
  Qchar = as.character(Q)

  res[[Qchar]] = data.table(selected_variables = names(temp), probs_selection = stab.glmnet$max,Q=Q)
  res[[Qchar]][, gene := gene_name]
  res[[Qchar]][, condition := cond_to_investigate]
  setcolorder(res[[Qchar]], c("gene", "condition", "Q","selected_variables"))
  print(Sys.time()-tim)
                }
  res=rbindlist(res)
  print(Sys.time()-tim0)
  return(res)
}

results = rbindlist(lapply(genes_names[1:length(genes_names) %% 100 == group_of_genes], run_glmnet))

dir.create(sprintf('%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/',EVO_IMMUNO_POP))
dir.create(sprintf('%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_withoutTranscription_withoutSVA/',EVO_IMMUNO_POP))
dir.create(sprintf('%s/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_withoutTranscription_withoutSVA/perm%s/',EVO_IMMUNO_POP,perm))

write.table(results, paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variableQ_globalFDR/glmnet_correlation_miRNA_geneRPKM_withoutTranscription_withoutSVA/perm",perm,"/gene_expression_glmnet_alpha0.5_cond", cond_to_investigate, "_", group_of_genes, ".tsv", sep=""), sep="\t", quote=F , row.names = F)
print("done")
