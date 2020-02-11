################################################################
##This script will contain a super_function to run go analysis##
################################################################

GO_analysis <- function(genes_to_inspect,
                        background_genes,
                        GO_data_base_to_use = "BHF_UCL",
                        covariables = NULL,
                        keep_genes = F){
    #genes_to_inspect = vector of the genes to inspect
    #background = vector of all the genes to be taking into the background (including genes_to_inspect)
    #GO_data_base_to_use = string matching to one or several GO_data_base_to_use
    #covariables : data.table containing in the a column "gene" with each genes in background and the values
    #of each covariables to be taken into account
    source(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/scripts/00a.GOannotation/01.Load_and_format_databases.R", sep=""))

    if(GO_data_base_to_use == "BHF_UCL"){
      GO_informations = load_BHF_UCL_miRNA()
    }else if(GO_data_base_to_use == "miRNABS"){
      GO_informations = load_Annotation_based_on_miRNABS()
    }else{
      stop("GO_data_base_to_use variable does not correspond to any known variable")
    }

    if(is.null(covariables)){
      data_to_measure = data.table(genes = background_genes)
      data_to_measure[, genes_in_group := ifelse(genes %in% genes_to_inspect, 1, 0)]
      #We only keep the lines were we have at least one annotation
      data_to_measure = data_to_measure[genes %in% GO_informations$miRNA]
      print(paste("we kept only", nrow(data_to_measure), "who have a GO annotation"))
      GO_informations = GO_informations[miRNA %in% background_genes]
      print(paste("we investigate", GO_informations[, length(unique(GO_term))], "GO annotations"))
      GO_to_tests = GO_informations[, unique(GO_term)]
      final = list()
      for (GO_term_to_test in GO_to_tests){
        data_to_measure[, inGo := ifelse(genes%in% GO_informations[GO_term %in% GO_term_to_test, miRNA], 1, 0)]
        results = summary(glm(data_to_measure, family = "binomial", formula = genes_in_group ~ inGo))
        coeffictient = results$coefficients["inGo", "Estimate"]
        standard_deviation =results$coefficients["inGo", "Std. Error"]
        pvalue = results$coefficients["inGo", "Pr(>|z|)"]
        oddsratio = exp(coeffictient)
        oddsratio_inf = exp(coeffictient-(2*standard_deviation))
        oddsratio_sup = exp(coeffictient+(2*standard_deviation))
        if( keep_genes){
          final[[GO_term_to_test]] = data.table(GO_term = GO_term_to_test, GO_name = GO_informations[GO_term %in% GO_term_to_test, unique(GO_meaning)], oddsratio, oddsratio_inf, oddsratio_sup, pvalue, genes_influenced = paste(data_to_measure[(inGo==1) & (genes_in_group == 1), unique(genes)], collapse = "//"))
        }else{
          final[[GO_term_to_test]] = data.table(GO_term = GO_term_to_test, GO_name = GO_informations[GO_term %in% GO_term_to_test, unique(GO_meaning)], oddsratio, oddsratio_inf, oddsratio_sup, pvalue)
        }
      }
      final = rbindlist(final)
      final = final[order(pvalue)]
      final[, fdr := p.adjust(pvalue, method = "fdr")]
      return(final)
    }
}
