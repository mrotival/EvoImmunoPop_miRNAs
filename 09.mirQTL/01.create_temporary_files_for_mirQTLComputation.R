inrt=function(x){ qnorm(rank(x,ties='random',na.last = NA)/(length(x)+1),mean(x,na.rm=T),sd(x,na.rm=T))}
for (condition in 1:5){
for (chromosome in 1:22){
    print(paste(condition, chromosome))
    genotypes = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/08.snp_data/genotypes/genotype_matrix_chr", chromosome, ".tsv", sep=""))
    snps_informations = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/08.snp_data/general_informations/snps_chr", chromosome, ".tsv", sep=""))
    #let's filter snps for MAF>5 in at least one population
    snps_informations = snps_informations[(MAF_EUB > 0.05) | (MAF_AFB > 0.05)]
    genotypes = genotypes[snp_id %in% snps_informations[, snp.name]]
    SNP_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/09.mirQTL/files_for_computation/genotype_matrix", chromosome,"_",condition,".tsv", sep="")
    ##expression
    corrected_count_transformed = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""), sep="\t")

    #We filter the names of column based on the condition
    miRNA_ID = corrected_count_transformed[, .(miRNA_ID = ID)]
    corrected_count_transformed = corrected_count_transformed[, (substr(colnames(corrected_count_transformed), 8,8) == as.character(condition)), with = F ]

    #Applying rank correction
    corrected_count_transformed_after_transformation = as.data.table(t(apply(corrected_count_transformed, 1, inrt)))
    colnames(corrected_count_transformed_after_transformation) = substr(colnames(corrected_count_transformed_after_transformation), 1,6)

    #recreating a proper expression table
    expression = cbind(miRNA_ID, corrected_count_transformed_after_transformation)

    ##miRNA positions
    miRNA_coordinate = fread(paste(EVO_IMMUNO_POP, "ERCPilot_SharedDBs/mirbase20/miRNA_mature_coordinates_strandinfo.bed", sep=""))
    names(miRNA_coordinate) = c("chromosome", "start", "end", "miRNA_name", "V5", "V6", "V7", "V8")


    ##Filter all of these files to only keep the correct genes and individuals
    genotypes = genotypes[, mget(c("snp_id", colnames(expression)[2:ncol(expression)]))]

    ##Les coordonées de certains miRNAs sont peut claires, on les retire ce cette analyse (on les utilisera peut être pour les trans)
    miRNA_coordinate = miRNA_coordinate[miRNA_name %in% expression$miRNA_ID] #not enough
    duplicated_miRNAs = miRNA_coordinate[duplicated(miRNA_name), unique(miRNA_name)]
    miRNA_coordinate = miRNA_coordinate[!(miRNA_name %in% duplicated_miRNAs)]

    expression = expression[!(miRNA_ID %in% duplicated_miRNAs)]

    #Covariates
    covariate_information = data.table(individual = colnames(expression[, 2:ncol(expression)]))
    covariate_information[, population := ifelse(substr(individual,1,3) == "EUB", 0, 1)]
    covariate_information = dcast(covariate_information, . ~individual, value.var = "population")
    colnames(covariate_information)[1] = "ID"
    covariate_information[1, ID := "population"]

    write.table(genotypes, SNP_file_name, sep="\t", quote = F, row.names = F)

    snps_informations = snps_informations[, .(snp.name, chromosome = paste("chr", chromosome, sep=""), position)]
    miRNA_coordinate = miRNA_coordinate[, .(miRNA_name, chromosome, start, end)]


    expression_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/09.mirQTL/files_for_computation/expression", chromosome,"_",condition,".tsv", sep="")
    covariate_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/09.mirQTL/files_for_computation/covariate", chromosome,"_",condition,".tsv", sep="")
    snps_informations_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/09.mirQTL/files_for_computation/snps_informations", chromosome,"_",condition,".tsv", sep="")
    miRNA_coordinate_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/09.mirQTL/files_for_computation/miRNA_coordinate", chromosome,"_",condition,".tsv", sep="")

    write.table(genotypes, SNP_file_name, sep="\t", quote = F, row.names = F)
    write.table( miRNA_coordinate, miRNA_coordinate_file_name, sep="\t", quote = F, row.names = F)
    write.table( snps_informations, snps_informations_file_name, sep="\t", quote = F, row.names = F)
    write.table( expression, expression_file_name, sep="\t", quote = F, row.names = F)
    write.table( covariate_information, covariate_file_name, sep="\t", quote = F, row.names = F)
  }
}
