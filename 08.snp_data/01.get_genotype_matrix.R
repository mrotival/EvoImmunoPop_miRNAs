###############################################################################
##Aim of this script : get a readable genotype matrix for each chromosome    ##
##We also want to get some filters (bi_allelic, MAF > 5% in at least one pop)##
###############################################################################

#############
##Libraries##
#############
require(snpStats)
total = list()
for (chromosome in 1:22){
  print(paste("treating chromosome :", chromosome))
  #################################
  ##Get genotypes from plink data##
  #################################
  plink_information=read.plink(paste(EVO_IMMUNO_POP,"/DATA_FREEZE/ERC_Main_Genotyping_24022015/Imputation/EvoImmunoPop_imputation_200x19619457_chr",chromosome,".bed",sep=""))
  genotype_matrix = 2-as(plink_information$genotype,'numeric')
  genotype_matrix = t(genotype_matrix)
  genotype_matrix = cbind(data.table(snp_id = row.names(genotype_matrix)), as.data.table(genotype_matrix))
  snps_information = as.data.table(plink_information$map)

  ########################
  ##Remove multi allelic##
  ########################
  duplicated_positions = snps_information[, .N, by = c("position")]
  duplicated_positions = duplicated_positions[N>1]
  snps_information = snps_information[!(position %in% duplicated_positions[, position])]
  snps_information[,cM := NULL]
  genotype_matrix = genotype_matrix[snp_id %in% snps_information[, snp.name]]


  ################################################
  ##Compute frequency of SNPs in each population##
  ################################################
  melted_genotypes = melt(genotype_matrix, id.vars = "snp_id", variable.name = "individual", value.name = "genotype")
  melted_genotypes[, population := substr(individual, 1,3)]
  MAF_information = melted_genotypes[, .(MAF = sum(genotype, na.rm = T)/(2*(.N - sum(is.na(genotype))))), by = c("snp_id", "population")]
  MAF_information[, MAF := pmin(MAF, 1-MAF)]
  MAF_information = dcast(MAF_information, snp_id ~ population, value.var = "MAF")
  snps_information[, MAF_EUB := MAF_information[match(snps_information$snp.name, snp_id), EUB]]
  snps_information[, MAF_AFB := MAF_information[match(snps_information$snp.name, snp_id), AFB]]

  write.table(genotype_matrix, sprintf("%s/Maxime/miRNA_V2/data/08.snp_data/genotypes/genotype_matrix_chr%s.tsv",EVO_IMMUNO_POP,chromosome), quote = F, row.names = F, sep="\t")
  write.table(snps_information, sprintf("%s/Maxime/miRNA_V2/data/08.snp_data/general_informations/snps_chr%s.tsv",EVO_IMMUNO_POP,chromosome), quote = F, row.names = F, sep="\t")
  total[[chromosome]] = snps_information
}
total = rbindlist(total)
write.table(total,sprintf("%s/Maxime/miRNA_V2/data/08.snp_data/general_informations/all_snps.tsv",EVO_IMMUNO_POP), quote = F, row.names = F, sep="\t" )
