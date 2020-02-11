#######################################################################
##Aims : get eQTLs and compute FDR limitations in CIS and TRANS cases##
#######################################################################

#############
##Libraries##
#############
require(ggplot2)

#############
##Load data##
#############

get_cis_eqtls <- function(chromosome, condition){
  temp = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/cis_isomiRQTL/isomirqtls", chromosome, "_", condition, "_0_matrixeQTL_output.tsv", sep=""))
  temp[, condition := condition]
  return(temp)
}
cis_eqtls = list()
for (chr in 1:22){
  for (cond in 1:5){
    print(paste("cis isomiRQTLs :", chr, cond))
    cis_eqtls[[paste(chr, cond)]] = get_cis_eqtls(chr, cond)
  }
}
cis_eqtls = rbindlist(cis_eqtls)


get_cis_permutations <- function(chromosome, condition, permutation_number){
  temp = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/cis_isomiRQTL/permutations/isomirqtls", chromosome, "_", condition, "_",permutation_number,"_matrixeQTL_output.tsv", sep=""))
  temp[, type := paste("permutation", permutation_number)]
  temp[, condition := condition]
  return(temp)
}
cis_perm = list()
for (chr in 1:22){
  for (cond in 1:5){
    for (perm in 1:100){
      print(paste("cis permutations :", chr, cond, perm))
      cis_perm[[paste(chr, cond, perm)]] = get_cis_permutations(chr, cond, perm)
    }
  }
}
cis_perm = rbindlist(cis_perm)

compute_fdr_for_pvalue_limit <-function(pvalue_limit){
  number_observation = cis_eqtls[pvalue<pvalue_limit, length(unique(gene))]
  estimate_of_false_observation = cis_perm[, .( N = sum(complete.cases(unique(ifelse(pvalue<pvalue_limit, gene, NA))))), by = type]
  # estimate_of_false_observation = cis_perm[pvalue<pvalue_limit, .(N = sum(pvalue<pvalue_limit)), by = type]
  estimate_of_false_observation[, estimate_fdr := N/number_observation]
  mean_fdr = mean(estimate_of_false_observation$estimate_fdr)
  median_fdr = median(estimate_of_false_observation$estimate_fdr)
  sup_fdr = quantile(estimate_of_false_observation$estimate_fdr, probs=0.975)
  inf_fdr = quantile(estimate_of_false_observation$estimate_fdr, probs=0.025)
  return(data.table(pvalue_limit, mean_fdr, median_fdr, sup_fdr, inf_fdr))
}
to_plot = list()
for (p in c(1e-2, 1e-3, 5e-3, 1e-4, 5e-4, 1e-5, 5e-5, 1e-6, 5e-6,1e-7,5e-7,1e-8,5e-8)){
  to_plot[[as.character(p)]] = compute_fdr_for_pvalue_limit(p)
}
to_plot = rbindlist(to_plot)
to_plot[, logpval := -log10(pvalue_limit)]

#Plotting FDR limits :
p <- ggplot(to_plot, aes(x = logpval, y = mean_fdr))
p <- p + geom_line() + geom_point()
p <- p + geom_ribbon(aes(ymax = sup_fdr, ymin = inf_fdr), fill = "darkgreen", alpha = 0.5)
p <- p + geom_hline(yintercept = 0.05, linetype="dotted", color = "Red")
p <- p + theme_bw()
p <- p + ylab("FDR Cis")
p <- p + xlab("pvalue cutoff")
pdf(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/figures/10.isomiRQTL/FDR_pvalue_relationship_in_cis_isomirQTL.pdf", sep=""))
print(p)
dev.off()

#Get fdr_permutation
cis_eqtls[, FDR_Perm := 1]
for (p in sort(c(1e-2, 1e-3, 5e-3, 1e-4, 5e-4, 1e-5, 5e-5, 1e-6, 5e-6,1e-7,5e-7,1e-8,5e-8), decreasing = T)){
  cis_eqtls[pvalue < p, FDR_Perm := to_plot[pvalue_limit == p, mean_fdr]]
}

write.table(cis_eqtls, paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/cis_isomiRQTLs_with_FDR.tsv", sep=""), quote = F, row.names = F, sep="\t")
cis_eqtls[, code := paste(gene, snps)]

filtered_cis_eqtl = cis_eqtls[1:.N]
filtered_cis_eqtl = filtered_cis_eqtl[order(pvalue)]
filtered_cis_eqtl = filtered_cis_eqtl[FDR_Perm<0.05]
ids_to_keep = filtered_cis_eqtl[, .(code = code[1]), by = gene]
filtered_cis_eqtl = cis_eqtls[code %in% ids_to_keep$code]


write.table(filtered_cis_eqtl, paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/10.isomiRQTL/cis_isomiRQTLs_with_FDR_filtered_best_isomiRNA_snp_association.tsv", sep=""), quote = F, row.names = F, sep="\t")

# ######################
# ##Now for Trans eqtl##
# ######################
#
#
# get_trans_eqtls <- function(chromosome, condition){
#   temp = fread(paste(EVO_IMMUNO_POP, "Martin/miRNA/12.isomiRQTL/results/trans/eqtls", chromosome, "_", condition, "_0.tsv", sep=""))
#   temp[, condition := condition]
#   return(temp)
# }
# trans_eqtls = list()
# for (chr in 1:22){
#   for (cond in 1:5){
#     print(paste("trans eqtls :", chr, cond))
#     trans_eqtls[[paste(chr, cond)]] = get_trans_eqtls(chr, cond)
#   }
# }
# trans_eqtls = rbindlist(trans_eqtls)
#
#
# get_trans_permutations <- function(chromosome, condition, permutation_number){
#   temp = fread(paste(EVO_IMMUNO_POP, "Martin/miRNA/12.isomiRQTL/results/trans/permutations/eqtls", chromosome, "_", condition, "_",permutation_number,".tsv", sep=""))
#   temp[, type := paste("permutation", permutation_number)]
#   temp[, condition := condition]
#   return(temp)
# }
# trans_perm = list()
# for (chr in 1:22){
#   for (cond in 1:5){
#     for (perm in 1:100){
#       print(paste("trans permutations :", chr, cond, perm))
#       trans_perm[[paste(chr, cond, perm)]] = get_trans_permutations(chr, cond, perm)
#     }
#   }
# }
# trans_perm = rbindlist(trans_perm)
#
# compute_fdr_for_pvalue_limit <-function(pvalue_limit){
#   number_observation = trans_eqtls[pvalue<pvalue_limit, length(unique(gene))]
#   estimate_of_false_observation = trans_perm[, .( N = sum(complete.cases(unique(ifelse(pvalue<pvalue_limit, gene, NA))))), by = type]
#   # estimate_of_false_observation = trans_perm[pvalue<pvalue_limit, .(N = sum(pvalue<pvalue_limit)), by = type]
#   estimate_of_false_observation[, estimate_fdr := N/number_observation]
#   mean_fdr = mean(estimate_of_false_observation$estimate_fdr)
#   median_fdr = median(estimate_of_false_observation$estimate_fdr)
#   sup_fdr = quantile(estimate_of_false_observation$estimate_fdr, probs=0.975)
#   inf_fdr = quantile(estimate_of_false_observation$estimate_fdr, probs=0.025)
#   return(data.table(pvalue_limit, mean_fdr, median_fdr, sup_fdr, inf_fdr))
# }
# to_plot = list()
# for (p in c(1e-2, 1e-3, 5e-3, 1e-4, 5e-4, 1e-5, 5e-5, 1e-6, 5e-6,1e-7,5e-7,1e-8,5e-8, 1e-9, 5e-9, 1e-10, 5e-10, 1e-11, 5e-11, 1e-12, 5e-12)){
#   to_plot[[as.character(p)]] = compute_fdr_for_pvalue_limit(p)
# }
# to_plot = rbindlist(to_plot)
# to_plot[, logpval := -log10(pvalue_limit)]
#
# #Plotting FDR limits :
# q <- ggplot(to_plot, aes(x = logpval, y = mean_fdr))
# q <- q + geom_line() + geom_point()
# q <- q + geom_ribbon(aes(ymax = sup_fdr, ymin = inf_fdr), fill = "darkgreen", alpha = 0.5)
# q <- q + geom_hline(yintercept = 0.05, linetype="dotted", color = "Red")
# q <- q + theme_bw()
# q <- q+ ylab("FDR trans")
# q <- q + xlab("pvalue cutoff")
#
#
# #Get fdr_permutation
# trans_eqtls[, FDR_Perm := 1]
# for (p in sort(c(1e-2, 1e-3, 5e-3, 1e-4, 5e-4, 1e-5, 5e-5, 1e-6, 5e-6,1e-7,5e-7,1e-8,5e-8, 1e-9, 5e-9, 1e-10, 5e-10, 1e-11, 5e-11, 1e-12, 5e-12), decreasing = T)){
#   if(min(trans_eqtls$pvalue)< p){
#     print(p)
#     trans_eqtls[pvalue < p, FDR_Perm := to_plot[pvalue_limit == p, mean_fdr]]
#   }
# }
#
# write.table(trans_eqtls, "../results/trans_isomirqtls.tsv", quote = F, row.names = F, sep="\t")
