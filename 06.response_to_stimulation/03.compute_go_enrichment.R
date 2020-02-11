#####################################################################################################
##Aim of the script : look for functionnal differences between genes upregulated in each conditions##
#####################################################################################################

######################################
##Librairies and necessary functions##
######################################
source(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/scripts/00a.GOannotation/02.run_GO_enrichment.R", sep=""))

#############
##Functions##
#############
GO_enrichment_in_genes_upregulated_in_one_condition <- function(padj_cutoff, condition, population,  pairing = "unpaired", type){
  expression_change_data = fread(paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/data/06.response_to_stimulation/", pairing, "/response_to_stimulation_basalVScond",condition, ifelse(population == "", "", "_"), population, ".tsv", sep=""))

  background_genes = expression_change_data[, unique(miRNA)]
  if(type == "overexpressed"){
    genes_impacted = expression_change_data[(padj < padj_cutoff) & (log2FoldChange>0), unique(miRNA)]
  }else if(type == "underexpressed"){
    genes_impacted = expression_change_data[(padj < padj_cutoff) & (log2FoldChange<0), unique(miRNA)]
  }else if(type == "over or under expressed"){
    genes_impacted = expression_change_data[(padj < padj_cutoff), unique(miRNA)]
  }else{
    stop("unknown type")
  }

  to_return = GO_analysis(genes_impacted, background_genes, GO_data_base_to_use ="BHF_UCL")
  to_return[, type := type]
  to_return[,  population := ifelse(population == "", "both", population)]
  to_return[,  condition := condition]

  return(to_return)
}
results = list()
for (c in 2:5){
  for (p in c("AFB", "EUB", "")){
    for (t in c("overexpressed", "underexpressed", "over or under expressed")){
      print(paste("classic", c, p, t))
      results[[paste(c, p, t)]] =GO_enrichment_in_genes_upregulated_in_one_condition(0.05, c, p, pairing='paired',type = t)
    }
  }
}

results = rbindlist(results)
write.table(results, paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/data/06.response_to_stimulation/go_enrichment_on_differentially_expressed_genes.paired.tsv",sep=""), sep="\t", quote = F, row.names = F)


GO_enrichment_in_genes_upregulated_in_one_condition_BS <- function(padj_cutoff, condition, population,  pairing = "unpaired", type){
  expression_change_data = fread(paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/data/06.response_to_stimulation/", pairing, "/response_to_stimulation_basalVScond",condition, ifelse(population == "", "", "_"), population, ".tsv", sep=""))

  background_genes = expression_change_data[, unique(miRNA)]
  if(type == "overexpressed"){
    genes_impacted = expression_change_data[(padj < padj_cutoff) & (log2FoldChange>0), unique(miRNA)]
  }else if(type == "underexpressed"){
    genes_impacted = expression_change_data[(padj < padj_cutoff) & (log2FoldChange<0), unique(miRNA)]
  }else if(type == "over or under expressed"){
    genes_impacted = expression_change_data[(padj < padj_cutoff), unique(miRNA)]
  }else{
    stop("unknown type")
  }

  to_return = GO_analysis(genes_impacted, background_genes, GO_data_base_to_use ="miRNABS")
  to_return[, type := type]
  to_return[,  population := ifelse(population == "", "both", population)]
  to_return[,  condition := condition]

  return(to_return)
}

results_BS = list()
for (c in 2:5){
  for (p in c("AFB", "EUB", "")){
    for (t in c("overexpressed", "underexpressed", "over or under expressed")){
      print(paste("BS", c, p, t))
      results_BS[[paste(c, p, t)]] =GO_enrichment_in_genes_upregulated_in_one_condition_BS(0.05, c, p, pairing='paired', type = t)
    }
  }
}
results_BS = rbindlist(results_BS)
results_BS[order(pvalue)][1]
write.table(results, paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/data/06.response_to_stimulation/go_enrichment_on_differentially_expressed_genes.paired.miRNABSBased.tsv", sep="\t", quote = F, row.names = F)
