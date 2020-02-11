#############
##Libraries##
#############
suppressMessages(require(GenomicRanges))

#############
##Load data##
#############
snps = fread(sprintf("%s/Maxime/miRNA_V2/data/08.snp_data/general_informations/all_snps.tsv",EVO_IMMUNO_POP))
hairpin_coordinates = fread(paste(EVO_IMMUNO_POP, "ERCPilot_SharedDBs/mirbase20/miRNA_hairpin_coordinates_strandinfo.bed", sep=""))
names(hairpin_coordinates) = c("chromosome", "hairpin_start", "hairpin_end", "mir1", "V5", "hairpin_strand", "hairpin")
hairpin_coordinates[, code := paste(hairpin, mir1)]
#Be carefull, the mi1 in the hairpin is ONLY the 5p (or 3p ?) version of the miRNA


mature_miRNA_coordinates = fread(paste(EVO_IMMUNO_POP, "ERCPilot_SharedDBs/mirbase20/miRNA_mature_coordinates_strandinfo.bed", sep=""))
names(mature_miRNA_coordinates) = c("chromosome", "start", "end", "mir1", "V5", "strand", "mir2", "hairpin")
mature_miRNA_coordinates[, code := paste(hairpin, mir1)]


############################
##Compute drosha cut sites##
############################
drosha_cut_sites = hairpin_coordinates[, .(start = c((hairpin_start - 22), hairpin_end+1),
                                           end = c(hairpin_start-1, hairpin_end + 22),
                                          chromosome = chromosome[1]), by = hairpin]

#############################
##Compute 5p and 3p regions##
#############################
hairpin_coordinates[, hairpin_middle := (hairpin_start + hairpin_end)/2]

hairpin_coordinates_5p_GR = GRanges(hairpin_coordinates$chromosome,
                                    IRanges(start = hairpin_coordinates[, ifelse(hairpin_strand == "+", hairpin_start, hairpin_middle)],
                                              end = hairpin_coordinates[, ifelse(hairpin_strand == "+", hairpin_middle, hairpin_end)]),
                                    strand = hairpin_coordinates$hairpin_strand)

hairpin_coordinates_3p_GR = GRanges(hairpin_coordinates$chromosome,
                                    IRanges(start = hairpin_coordinates[, ifelse(hairpin_strand == "+", hairpin_middle, hairpin_start)],
                                              end = hairpin_coordinates[, ifelse(hairpin_strand == "+", hairpin_end, hairpin_middle)]),
                                    strand = hairpin_coordinates$hairpin_strand)

mature_miRNA_coordinates_GR = GRanges(mature_miRNA_coordinates$chromosome,
                              IRanges(start = mature_miRNA_coordinates$start,
                                      end = mature_miRNA_coordinates$end), strand = mature_miRNA_coordinates$strand)

ov_5p = as.data.table(findOverlaps(mature_miRNA_coordinates_GR, hairpin_coordinates_5p_GR, minoverlap = 10))

arm_5p = mature_miRNA_coordinates[ov_5p$queryHits, .(hairpin, mir1, start, end, chromosome)]

ov_3p = as.data.table(findOverlaps(mature_miRNA_coordinates_GR, hairpin_coordinates_3p_GR, minoverlap = 10))

arm_3p = mature_miRNA_coordinates[ov_3p$queryHits, .(hairpin, mir1, start, end, chromosome)]

################################
##Give the information to snps##
################################
snps_GR = GRanges(paste("chr", snps$chromosome, sep=""), IRanges(start = snps$position, width = 1))

hairpin_GR = GRanges(hairpin_coordinates$chromosome,
                                    IRanges(start = hairpin_coordinates[, hairpin_start],
                                              end = hairpin_coordinates[, hairpin_end]),
                                    strand = hairpin_coordinates$hairpin_strand)

# wUpstream_5p = which(ifelse(hairpin_annot[, hairpin_strand]=='+',
# 							 !is.na(hairpin_annot[, mat_start_5p]),
# 							 !is.na(hairpin_annot[, mat_end_5p])))
# 
# hairpin_upstream5p_GR = GRanges(hairpin_annot$hairpin_chr[wUpstream_5p],
#                                     IRanges(start = ifelse(hairpin_annot[wUpstream_5p, hairpin_strand]=='+',hairpin_annot[wUpstream_5p, hairpin_start],hairpin_annot[wUpstream_5p, mat_end_5p]),
#                                               end = ifelse(hairpin_annot[wUpstream_5p, hairpin_strand]=='+',hairpin_annot[wUpstream_5p, mat_start_5p],hairpin_annot[wUpstream_5p,hairpin_end])),
#                                     strand = hairpin_annot$hairpin_strand[wUpstream_5p])
# 
# wLoop=which( ifelse(hairpin_annot[, hairpin_strand]=='+',
# 				!is.na(hairpin_annot[, mat_end_5p]) & !is.na(hairpin_annot[, mat_start_3p]),
# 				!is.na(hairpin_annot[, mat_start_5p]) & !is.na(hairpin_annot[, mat_end_3p])))
# 
# hairpin_loop_GR = GRanges(hairpin_annot$hairpin_chr[wLoop],
#                                     IRanges(start = ifelse(hairpin_annot[wLoop, hairpin_strand]=='+',hairpin_annot[wLoop, mat_end_5p]+1,hairpin_annot[wLoop, mat_start_3p]+1),
#                                               end = ifelse(hairpin_annot[wLoop, hairpin_strand]=='+',hairpin_annot[wLoop, mat_start_3p]-1,hairpin_annot[wLoop, mat_end_5p]-1)),
#                                     strand = hairpin_annot$hairpin_strand[wLoop])
# 
# wDownstream_3p = which(ifelse(hairpin_annot[, hairpin_strand]=='+',
# 							 !is.na(hairpin_annot[, mat_end_3p]),
# 							 !is.na(hairpin_annot[, mat_start_3p])))
# hairpin_downstream3p_GR = GRanges(hairpin_annot$hairpin_chr[wDownstream_3p],
#                                     IRanges(start = ifelse(hairpin_annot[wDownstream_3p, hairpin_strand]=='+',hairpin_annot[wDownstream_3p, mat_end_3p]+1,hairpin_annot[wDownstream_3p, hairpin_start]),
#                                               end = ifelse(hairpin_annot[wDownstream_3p, hairpin_strand]=='+',hairpin_annot[wDownstream_3p, hairpin_end],hairpin_annot[wDownstream_3p,mat_start_3p])),
#                                     strand = hairpin_annot$hairpin_strand[wDownstream_3p])

drosha_cut_sites_GR = GRanges(drosha_cut_sites$chromosome, IRanges(start = drosha_cut_sites$start,
                                                                    end = drosha_cut_sites$end))
arm_5p_GR = GRanges(arm_5p$chromosome, IRanges(start = arm_5p$start,
                                                end = arm_5p$end))

arm_3p_GR = GRanges(arm_3p$chromosome, IRanges(start = arm_3p$start,
                                                end = arm_3p$end))

ov = as.data.table(findOverlaps(hairpin_GR, snps_GR))
snps[, in_hairpin := 1:.N %in% ov[, subjectHits]]

ov = as.data.table(findOverlaps(drosha_cut_sites_GR, snps_GR))
snps[, in_drosha_cut_region := 1:.N %in% ov[, subjectHits]]

ov = as.data.table(findOverlaps(arm_5p_GR, snps_GR))
snps[, in_5p_miRNA := 1:.N %in% ov[, subjectHits]]

ov = as.data.table(findOverlaps(arm_3p_GR,snps_GR ))
snps[, in_3p_miRNA := 1:.N %in% ov[, subjectHits]]

miRNA_TSS=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_TSS_DeRie_natureBioTech_2017_TableS15.txt',EVO_IMMUNO_POP))
colnames(miRNA_TSS)=make.names(colnames(miRNA_TSS))

miRNA_TSS_GR=GRanges(miRNA_TSS$Chromosome, IRanges(start=miRNA_TSS$TSS-500,
end=miRNA_TSS$TSS+500),strand=miRNA_TSS$Strand)


ov = as.data.table(findOverlaps(miRNA_TSS_GR,snps_GR ))
snps[, in_TSS_500 := 1:.N %in% ov[, subjectHits]]


snps[,mir_cat:=ifelse(in_5p_miRNA,'5p arm',
					ifelse(in_3p_miRNA,'3p arm',
						ifelse(in_hairpin,'loop',
							ifelse(in_drosha_cut_region,'drosha',
								ifelse(in_TSS_500,'TSS_500','')))))]

# ov = as.data.table(findOverlaps(hairpin_downstream3p_GR,snps_GR ))
# snps[, down_3p_miRNA := 1:.N %in% ov[, subjectHits]]
# 
# ov = as.data.table(findOverlaps(hairpin_upstream5p_GR,snps_GR ))
# snps[, up_5p_miRNA := 1:.N %in% ov[, subjectHits]]
# 
# ov = as.data.table(findOverlaps(hairpin_loop_GR,snps_GR ))
# snps[, in_loop := 1:.N %in% ov[, subjectHits]]

write.table(snps,sprintf("%s/Maxime/miRNA_V2/data/08.snp_data/general_informations/all_snps.tsv",EVO_IMMUNO_POP), quote = F, row.names = F, sep="\t" )

snps=snps[which(MAF_EUB+MAF_AFB>.1),]
fwrite(snps,file=sprintf("%s/Maxime/miRNA_V2/data/08.snp_data/general_informations/all_snps_MAFover5.tsv",EVO_IMMUNO_POP), sep="\t" )



