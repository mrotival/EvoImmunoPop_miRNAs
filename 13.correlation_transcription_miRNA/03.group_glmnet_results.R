###################################################################################################
##Processing of the results of correlation_allmiRNAs_geneFPKM_transcriptionadjusted_svaadjusted.R##
###################################################################################################
all_results = list()
for (condition in 1:5){
  for (part in 0:19){
    print(paste(condition, part))
    all_results[[paste(condition, part)]] = fread(paste(EVO_IMMUNO_POP,
                                                        "/Maxime/miRNA_V2/data/13.correlation_transcription_miRNA/glmnet_miRNA_intronReads_withoutSVA/",
                                                        "gene_transcription_glmnet_alpha0.5_cond",condition,"_",part,".tsv", sep=""))
  }
}
all_results = rbindlist(all_results)
write.table(all_results, paste(EVO_IMMUNO_POP,
                  "/Maxime/miRNA_V2/data/13.correlation_transcription_miRNA/glmnet_miRNA_intronReads_withoutSVA/",
                  "gene_transcription_glmnet_alpha0.5_allConds.tsv", sep=""), sep="\t", quote = F)
