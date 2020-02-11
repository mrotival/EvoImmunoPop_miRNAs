#################################################################
##Aim of the script : group all the expression data in one file##
#################################################################

#############
##Load data##
#############
sample_informations = fread(paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/data/02.pre_processing/filtered_and_combined_samples.tsv", sep=""))

#################
##Main Function##
#################

get_after_transpose_mature_mirna_count <- function(library_ID){
  print(library_ID)
  results = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03b.isomirs_alignment/alignement/after_transpose_and_substitutions/",library_ID, ".maturemiRNA.intersect.transposed.substituted.bed", sep="" ))
  names(results) = c("chromosome", "start", "end", "strand", "sequence", "count", "subsitutions", "mir1", "mir2", "hairpin")
  results[, library_ID := library_ID]
  return(results)
}

results = rbindlist(lapply(sample_informations$library_ID, get_after_transpose_mature_mirna_count))

write.table(results, paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03b.isomirs_alignment/miRNA_transcripts_raw_counts_melded.tsv",sep=""), quote = F, row.names = F , sep="\t")
