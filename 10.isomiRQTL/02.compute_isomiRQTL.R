##############
##Librairies##
##############
require(MatrixEQTL)
require(snpStats)

args = commandArgs(trailingOnly=TRUE)
#######################
##Parameters variable##
#######################
##choice of the model
useModel = modelLINEAR # modelANOVA or modelLINEAR or modelLINEAR_CROSS

##Association with pvalue over this THreshold won't be outputted
pvOutputThreshold = 1e-2;

##Covariance matrix for error term
errorCovariance = numeric()

##distance for cis definition
cisDist = 1e6

##chromosome impacted
chromosome = as.numeric(args[1])

##condition tested
condition = as.numeric(args[2])


##permutation or not
permutation = as.numeric(args[3])

##################################################
##Loading data and format it fo matrixeQTL usage##
##################################################
SNP_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/files_for_computation/genotype_matrix", chromosome,"_",condition,".tsv", sep="")
expression_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/files_for_computation/expression", chromosome,"_",condition,".tsv", sep="")
covariate_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/files_for_computation/covariate", chromosome,"_",condition,".tsv", sep="")
snps_informations_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/files_for_computation/snps_informations", chromosome,"_",condition,".tsv", sep="")
miRNA_coordinate_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/files_for_computation/miRNA_coordinate", chromosome,"_",condition,".tsv", sep="")

snps_informations = fread(snps_informations_file_name)
miRNA_coordinate = fread(miRNA_coordinate_file_name)

######################
##Main eQTL analysis##
######################
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name );

## Load gene expression data
if (permutation == 0){ #Classic run
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
}else{ #permutation_run
  set.seed(permutation)
  non_shuffled_genes = fread(expression_file_name)
  colnames_to_shuffle = names(non_shuffled_genes[1:.N])[2:ncol(non_shuffled_genes)]
  shuffled_colnames = c("isomir_ID", sample(colnames_to_shuffle))
  shuffled_gene_expression = non_shuffled_genes[, mget(shuffled_colnames)]
  names(shuffled_gene_expression) = names(non_shuffled_genes)
  expression_file_name = paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/files_for_computation/miRNA_coordinate", chromosome,"_",condition, "_", permutation,".tsv", sep="")
  write.table(shuffled_gene_expression, expression_file_name, quote = F, row.names = F, sep = "\t")
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
}


## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariate_file_name)>0) {
cvrt$LoadFile(covariate_file_name);
}

me = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name  = tempfile(),
      pvOutputThreshold  = 1e-8,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = FALSE,
      output_file_name.cis = tempfile(),
      pvOutputThreshold.cis = 1e-2,
      snpspos = as.data.frame(snps_informations),
      genepos = as.data.frame(miRNA_coordinate),
      cisDist = cisDist,
      pvalue.hist = "qqplot",
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE);

#2146, 8
#1555, 1
#1303, 6
print(length(unique((me$cis$eqtls[me$cis$eqtls$pvalue<1e-4,]$gene))))
print(length(unique((me$cis$eqtls[me$cis$eqtls$pvalue<1e-5,]$gene))))

if(permutation!=0){
  file.remove(expression_file_name)
  write.table(me$cis$eqtls, paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/cis_isomiRQTL/permutations/isomirqtls", chromosome,"_",condition, "_", permutation,"_matrixeQTL_output.tsv", sep=""), quote = F, row.names = F, sep="\t")
  write.table(me$trans$eqtls, paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/trans_isomiRQTL/permutations/isomirqtls", chromosome,"_",condition, "_", permutation,"_matrixeQTL_output.tsv", sep=""), quote = F, row.names = F, sep="\t")
}else{
  write.table(me$cis$eqtls, paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/cis_isomiRQTL/isomirqtls", chromosome,"_",condition, "_", permutation,"_matrixeQTL_output.tsv", sep=""), quote = F, row.names = F, sep="\t")
  write.table(me$trans$eqtls, paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/trans_isomiRQTL/isomirqtls", chromosome,"_",condition, "_", permutation,"_matrixeQTL_output.tsv", sep=""), quote = F, row.names = F, sep="\t")
}
