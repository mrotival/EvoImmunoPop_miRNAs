# [MR4]Violin plot:Comparison of % explained by transcription across conditions
# [MR5]Violin plot:Comparison of % explained by miRNAs across conditions

#WARNING Must be modified

condition_name = c("1"="NS", "2"="LPS", "3" = "PAM3CSK4", "4" ="R848", "5" = "IAV")
colERC5=c("#525252AA", "#E31A1CAA", "#33A02CAA", "#1F78B4AA", "#6A3D9AAA")
names(colERC5) = c("NS", "LPS", "PAM3CSK4", "R848", "IAV")
require(ggplot2)


all_CARSCORES = list()
for (cond in 1:5){
  all_CARSCORES[[cond]] = fread(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/12.correlation_gene_miRNAs/variance_explained_by_miRNA_CARScore_withoutSVA/variance_explained_by_miRNA_CARScore_",cond,".tsv", sep=""))
}
all_CARSCORES = rbindlist(all_CARSCORES)

all_CARSCORES[, condition := condition_name[condition]]
all_CARSCORES[, condition := factor(condition, levels = c("NS", "LPS", "PAM3CSK4", "R848", "IAV"))]

all_CARSCORES = all_CARSCORES[!is.na(variance_explained_by_pop)]

miRNABS = fread(paste(EVO_IMMUNO_POP, "Martin/miRNA/11.miRNA_binding_sites/results/miRNA_binding_sites_on_protein_coding_3primeUTR_simplified.tsv", sep=""))
miRNABS[, code := paste(EnsemblGeneID, gsub("-", "_", miRNA))]

colnames_to_delete = c(paste("miRNA_variance", 9:10, sep="_"),
                        paste("miRNA", 9:10, sep="_"),
                        paste("miRNA_carscore_sign", 9:10, sep="_"))

all_CARSCORES[, eval(colnames_to_delete) := NULL]
all_CARSCORES[, code_miRNA_1 := paste(gene, miRNA_1)]
all_CARSCORES[, code_miRNA_2 := paste(gene, miRNA_2)]
all_CARSCORES[, code_miRNA_3 := paste(gene, miRNA_3)]
all_CARSCORES[, code_miRNA_4 := paste(gene, miRNA_4)]
all_CARSCORES[, code_miRNA_5 := paste(gene, miRNA_5)]
all_CARSCORES[, code_miRNA_6 := paste(gene, miRNA_6)]
all_CARSCORES[, code_miRNA_7 := paste(gene, miRNA_7)]
all_CARSCORES[, code_miRNA_8 := paste(gene, miRNA_8)]

all_CARSCORES[is.na(miRNA_carscore_sign_1), miRNA_carscore_sign_1 := 0]
all_CARSCORES[is.na(miRNA_carscore_sign_2), miRNA_carscore_sign_2 := 0]
all_CARSCORES[is.na(miRNA_carscore_sign_3), miRNA_carscore_sign_3 := 0]
all_CARSCORES[is.na(miRNA_carscore_sign_4), miRNA_carscore_sign_4 := 0]
all_CARSCORES[is.na(miRNA_carscore_sign_5), miRNA_carscore_sign_5 := 0]
all_CARSCORES[is.na(miRNA_carscore_sign_6), miRNA_carscore_sign_6 := 0]
all_CARSCORES[is.na(miRNA_carscore_sign_7), miRNA_carscore_sign_7 := 0]
all_CARSCORES[is.na(miRNA_carscore_sign_8), miRNA_carscore_sign_8 := 0]



all_CARSCORES[, variance_explained_by_negative_miRNA := ifelse(miRNA_carscore_sign_1 == -1, miRNA_variance_1, 0) +
                                                         ifelse(miRNA_carscore_sign_2 == -1, miRNA_variance_2, 0) +
                                                         ifelse(miRNA_carscore_sign_3 == -1, miRNA_variance_3, 0)+
                                                         ifelse(miRNA_carscore_sign_4 == -1, miRNA_variance_4, 0)+
                                                         ifelse(miRNA_carscore_sign_5 == -1, miRNA_variance_5, 0)+
                                                         ifelse(miRNA_carscore_sign_6 == -1, miRNA_variance_6, 0)+
                                                         ifelse(miRNA_carscore_sign_7 == -1, miRNA_variance_7, 0)+
                                                         ifelse(miRNA_carscore_sign_8 == -1, miRNA_variance_8, 0)]

all_CARSCORES[, variance_explained_by_negative_miRNA_and_miRNABS := ifelse((miRNA_carscore_sign_1 == -1) & (code_miRNA_1 %in% miRNABS$code) , miRNA_variance_1, 0) +
                                                                    ifelse((miRNA_carscore_sign_2 == -1)& (code_miRNA_2 %in% miRNABS$code), miRNA_variance_2, 0) +
                                                                    ifelse((miRNA_carscore_sign_3 == -1)& (code_miRNA_3 %in% miRNABS$code), miRNA_variance_3, 0)+
                                                                    ifelse((miRNA_carscore_sign_4 == -1)& (code_miRNA_4 %in% miRNABS$code), miRNA_variance_4, 0)+
                                                                    ifelse((miRNA_carscore_sign_5 == -1)& (code_miRNA_5 %in% miRNABS$code), miRNA_variance_5, 0)+
                                                                    ifelse((miRNA_carscore_sign_6 == -1)& (code_miRNA_6 %in% miRNABS$code), miRNA_variance_6, 0)+
                                                                    ifelse((miRNA_carscore_sign_7 == -1)& (code_miRNA_7 %in% miRNABS$code), miRNA_variance_7, 0)+
                                                                    ifelse((miRNA_carscore_sign_8 == -1)& (code_miRNA_8 %in% miRNABS$code), miRNA_variance_8, 0)]


# #Transcriptions
p <- ggplot(all_CARSCORES, aes(x = condition, y = variance_explained_by_transcription, fill = condition))
p <- p + geom_violin(scale = "width")
p <- p + geom_boxplot(width = 0.1, outlier.shape = NA, notch = T)
p <- p + scale_fill_manual(values = colERC5)
p <- p + ylab("expression variation explained by transcription")
p <- p + guides(fill=FALSE)
p <- p + theme_bw()

pdf(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/with_transcription_withoutSVA/", "expression_variation_explaind_by_transcription.pdf", sep=""), width = 4, height =4)
print(p)
dev.off()

#miRNAs
p <- ggplot(all_CARSCORES, aes(x = condition, y = variance_explained_by_miRNAs, fill = condition))
p <- p + geom_violin(scale = "width")
p <- p + geom_boxplot(width = 0.1, outlier.shape = NA, notch = T)
p <- p + scale_fill_manual(values = colERC5)
p <- p + ylab("expression variation explained by miRNAs")
p <- p + guides(fill=FALSE)
p <- p + theme_bw()

pdf(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/with_transcription_withoutSVA/", "expression_variation_explaind_by_miRNAs.pdf", sep=""), width = 4, height =4)
print(p)
dev.off()


all_CARSCORES_with_miRNAs = all_CARSCORES[number_miRNA_in_model!=0]
p <- ggplot(all_CARSCORES_with_miRNAs, aes(x = condition, y = variance_explained_by_miRNAs, fill = condition))
p <- p + geom_violin(scale = "width")
p <- p + geom_boxplot(width = 0.1, outlier.shape = NA, notch = T)
p <- p + scale_fill_manual(values = colERC5)
p <- p + ylab("expression variation explained by miRNAs")
p <- p + guides(fill=FALSE)
p <- p + theme_bw()

pdf(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/with_transcription_withoutSVA/", "expression_variation_explaind_by_miRNAs_only_genes_with_miRNA.pdf", sep=""), width = 4, height =4)
print(p)
dev.off()

all_CARSCORES_with_miRNAs_negative = all_CARSCORES_with_miRNAs[variance_explained_by_negative_miRNA != 0]
p <- ggplot(all_CARSCORES_with_miRNAs_negative, aes(x = condition, y = variance_explained_by_negative_miRNA, fill = condition))
p <- p + geom_violin(scale = "width")
p <- p + geom_boxplot(width = 0.1, outlier.shape = NA, notch = T)
p <- p + scale_fill_manual(values = colERC5)
p <- p + ylab("expression variation explained by miRNAs")
p <- p + guides(fill=FALSE)
p <- p + theme_bw()

pdf(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/with_transcription_withoutSVA/", "expression_variation_explaind_by_negative_miRNAs_only_genes_with_miRNA.pdf", sep=""), width = 4, height =4)
print(p)
dev.off()


all_CARSCORES_with_miRNAs_negative_binding = all_CARSCORES_with_miRNAs[variance_explained_by_negative_miRNA_and_miRNABS != 0]
p <- ggplot(all_CARSCORES_with_miRNAs_negative_binding, aes(x = condition, y = variance_explained_by_negative_miRNA_and_miRNABS, fill = condition))
p <- p + geom_violin(scale = "width")
p <- p + geom_boxplot(width = 0.1, outlier.shape = NA, notch = T)
p <- p + scale_fill_manual(values = colERC5)
p <- p + ylab("expression variation explained by miRNAs")
p <- p + guides(fill=FALSE)
p <- p + theme_bw()

pdf(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/with_transcription_withoutSVA/", "expression_variation_explaind_by_negative_and_binding_miRNAs_only_genes_with_miRNA.pdf", sep=""), width = 4, height =4)
print(p)
dev.off()



##Percent of genes regulated in different conditions
to_plot = all_CARSCORES[, .(percent_of_genes_with_miRNA = sum(variance_explained_by_miRNAs>0)/.N,
                            percent_of_genes_with_negative_miRNA = sum(variance_explained_by_negative_miRNA>0)/.N,
                          percent_of_genes_with_negative_with_BS_miRNA = sum(variance_explained_by_negative_miRNA_and_miRNABS>0)/.N), by = condition]

p <- ggplot(to_plot, aes(x = condition, y = percent_of_genes_with_miRNA, fill = condition))
p <- p + geom_bar(stat="identity")
p <- p + scale_fill_manual(values = colERC5)
p <- p + ylab("Percent genes regulated by miRNAs")
p <- p + guides(fill=FALSE)
p <- p + theme_bw()

pdf(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/with_transcription_withoutSVA/", "percent_genes_regulated_by_miRNAs.pdf", sep=""), width = 4, height =4)
print(p)
dev.off()


p <- ggplot(to_plot, aes(x = condition, y = percent_of_genes_with_negative_miRNA, fill = condition))
p <- p + geom_bar(stat="identity")
p <- p + scale_fill_manual(values = colERC5)
p <- p + ylab("Percent genes regulated by negative miRNAs")
p <- p + guides(fill=FALSE)
p <- p + theme_bw()


pdf(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/with_transcription_withoutSVA/", "percent_genes_regulated_by_negative_miRNAs.pdf", sep=""), width = 4, height =4)
print(p)
dev.off()



p <- ggplot(to_plot, aes(x = condition, y = percent_of_genes_with_negative_with_BS_miRNA, fill = condition))
p <- p + geom_bar(stat="identity")
p <- p + scale_fill_manual(values = colERC5)
p <- p + ylab("Percent genes regulated by negative binding miRNAs")
p <- p + guides(fill=FALSE)
p <- p + theme_bw()

pdf(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/figures/12.correlation_gene_miRNAs/with_transcription_withoutSVA/", "percent_genes_regulated_by_negative_binding_miRNAs.pdf", sep=""), width = 4, height =4)
print(p)
dev.off()
