inrt=function(x){ qnorm(rank(x,ties='random',na.last = T)/(length(x)+1),mean(x,na.rm=T),sd(x,na.rm=T))}
for (condition in 1:5){
  for (chromosome in 1:22){
    print(paste(condition, chromosome))
    genotypes = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/08.snp_data/genotypes/genotype_matrix_chr", chromosome, ".tsv", sep=""))
    snps_informations = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/08.snp_data/general_informations/snps_chr", chromosome, ".tsv", sep=""))
    #let's filter snps for MAF>5 in at least one population
    snps_informations = snps_informations[(MAF_EUB > 0.05) | (MAF_AFB > 0.05)]
    genotypes = genotypes[snp_id %in% snps_informations[, snp.name]]
    SNP_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/files_for_computation/genotype_matrix", chromosome,"_",condition,".tsv", sep="")

    #expression
    isomir_ratio_corrected = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/05.isomirs_count_correction/isomiR_ratios_aggregated_nosubs.GCRL_Batch_lane_corrected.tsv", sep=""))

    isomir_names = isomir_ratio_corrected[, .(mirID, isomir_ID)]

    #we filter on expression
    isomir_ratio_corrected = isomir_ratio_corrected[, mget(colnames(isomir_ratio_corrected)[substr(names(isomir_ratio_corrected), 8,8) == as.character(condition)])]
    names(isomir_ratio_corrected) = substr(names(isomir_ratio_corrected), 1,6)
    isomir_ratio_corrected = as.data.table(t(apply(isomir_ratio_corrected, 1, inrt)))
    expression = cbind(isomir_names, isomir_ratio_corrected)

    ##miRNA positions
    miRNA_coordinate = fread(paste(EVO_IMMUNO_POP, "ERCPilot_SharedDBs/mirbase20/miRNA_mature_coordinates_strandinfo.bed", sep=""))
    names(miRNA_coordinate) = c("chromosome", "start", "end", "miRNA_name", "V5", "V6", "V7", "V8")

    ##Filter all of these files to only keep the correct genes and individuals
    genotypes = genotypes[, mget(c("snp_id", colnames(expression)[3:ncol(expression)]))]

    ##Les coordonées de certains miRNAs sont peut claires, on les retire ce cette analyse (on les utilisera peut être pour les trans)

    duplicated_miRNAs = miRNA_coordinate[duplicated(miRNA_name), unique(miRNA_name)]
    miRNA_coordinate = miRNA_coordinate[!(miRNA_name %in% duplicated_miRNAs)]

    coordinate_isomirs = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_FULL_V2.0.tsv", sep=""), select = c("ID", "mirID", "hsa_ID"))

    coordinate_isomirs[, chromosome := miRNA_coordinate[match(coordinate_isomirs$hsa_ID, miRNA_name), chromosome]]
    coordinate_isomirs[, start := miRNA_coordinate[match(coordinate_isomirs$hsa_ID, miRNA_name), start]]
    coordinate_isomirs[, end := miRNA_coordinate[match(coordinate_isomirs$hsa_ID, miRNA_name), end]]
    coordinate_isomirs = coordinate_isomirs[complete.cases(coordinate_isomirs)]

    expression = expression[isomir_ID %in% coordinate_isomirs$ID]
    coordinate_isomirs = coordinate_isomirs[ID %in% expression$isomir_ID]



    #Covariates
    covariate_information = data.table(individual = colnames(expression[, 3:ncol(expression)]))
    covariate_information[, population := ifelse(substr(individual,1,3) == "EUB", 0, 1)]
    covariate_information = dcast(covariate_information, . ~individual, value.var = "population")
    colnames(covariate_information)[1] = "ID"
    covariate_information[1, ID := "population"]


    snps_informations = snps_informations[, .(snp.name, chromosome = paste("chr", chromosome, sep=""), position)]
    expression[, mirID := NULL]
    coordinate_isomirs = coordinate_isomirs[, .(ID, chromosome, start, end)]

    expression_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/files_for_computation/expression", chromosome,"_",condition,".tsv", sep="")
    covariate_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/files_for_computation/covariate", chromosome,"_",condition,".tsv", sep="")
    snps_informations_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/files_for_computation/snps_informations", chromosome,"_",condition,".tsv", sep="")
    miRNA_coordinate_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/files_for_computation/miRNA_coordinate", chromosome,"_",condition,".tsv", sep="")
    write.table( genotypes, SNP_file_name, sep="\t", quote = F, row.names = F)
    write.table( coordinate_isomirs, miRNA_coordinate_file_name, sep="\t", quote = F, row.names = F)
    write.table( snps_informations, snps_informations_file_name, sep="\t", quote = F, row.names = F)
    write.table( expression, expression_file_name, sep="\t", quote = F, row.names = F)
    write.table( covariate_information, covariate_file_name, sep="\t", quote = F, row.names = F)
  }
}
