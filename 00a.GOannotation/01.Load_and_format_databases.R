######################################################################################
##This script countain all the functions to load, and pretreat databases of GO terms##
######################################################################################
##All of the functions should return something like this
##mirna  GO_term  ontology_type GO_meaning
##hsa-miR-3 GO:XXXXX BP "affects the response to almond in bloodstreem"


load_BHF_UCL_miRNA <- function(){
  raw_data_set = suppressWarnings(fread(paste(EVO_IMMUNO_POP, "Martin/miRNA/Ressources/GOA/goa_human_rna.gaf", sep=""), showProgress = F))
  colnames(raw_data_set) = c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_term",
                            "DB_Reference", "Evidence Code", "with_or_from",
                            "ontology_type", "DB_Object_Name", "DB_Object_Synonym",
                              "DB_Object_Type", "Taxon", "Interacting_Taxon_ID",
                              "Date", "Assigned_By", "Annotation_Extension")

  ##We keep only miRNA
  miRNA_data_set = raw_data_set[DB_Object_Type == "miRNA"]



  #We create a "miRNA" column from the "DB_Object_Name" column
  miRNA_data_set = miRNA_data_set[, miRNA := substr(DB_Object_Name, 22, nchar(DB_Object_Name))]
  #We have a few that sill contain "microRNA"
  miRNA_data_set = miRNA_data_set[grepl("microRNA", miRNA), miRNA := substr(miRNA, 10, nchar(miRNA))]


  #We get a column GO_meaning
  GO_terms_and_names = fread(paste(EVO_IMMUNO_POP, "Martin/miRNA/Ressources/GO_annotations/GO_term_and_meaning.tsv", sep=""))
  miRNA_data_set[, GO_meaning := GO_terms_and_names[match(miRNA_data_set$GO_term, GO_term), GO_meaning]]
  miRNA_data_set = miRNA_data_set[, .(miRNA, GO_term, ontology_type, GO_meaning)]

  miRNA_data_set = miRNA_data_set[ !(GO_meaning %in% c("molecular_function", "cellular_component", "biological_process"))]
  return(miRNA_data_set)
}

load_Annotation_based_on_correlation_with_gene_corrected_by_transcription <- function(){
  raw_data_set = fread(paste(EVO_IMMUNO_POP, "Martin/miRNA/18.Gene_regulation_by_miRNA/results/GOEnrichments/significantGOEnrichments_correlation_with_expression_corrected_pop_transcription.tsv", sep=""))
  miRNA_data_set = raw_data_set[, .(miRNA, GO_term = category, ontology_type = ontology, GO_meaning = term)]
}

load_Annotation_based_on_correlation_with_gene_corrected_by_transcription_and_miRNABS <- function(){
  raw_data_set = fread(paste(EVO_IMMUNO_POP, "Martin/miRNA/18.Gene_regulation_by_miRNA/results/GOEnrichments/significantGOEnrichments_correlation_with_expression_and_BS_corrected_pop_transcription.tsv", sep=""))
  miRNA_data_set = raw_data_set[, .(miRNA, GO_term = category, ontology_type = ontology, GO_meaning = term)]
}

load_Annotation_based_on_miRNABS <- function(){
  raw_data_set = fread(paste(EVO_IMMUNO_POP, "Martin/miRNA/18.Gene_regulation_by_miRNA/results/GOEnrichments/significantGOEnrichments_BS.tsv", sep=""))
  miRNA_data_set = raw_data_set[, .(miRNA, GO_term = category, ontology_type = ontology, GO_meaning = term)]
  return(miRNA_data_set)
}
