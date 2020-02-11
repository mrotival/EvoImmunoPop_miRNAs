#############
##Libaririe##
#############
suppressMessages(require(GenomicRanges))

##################
##Get miRNA data##
##################

#We only look at miRNA that are detected by previous analysis
miRNAs = fread(sprintf("%s/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.RPM.GCRL_Batch_corrected_V2.0_MR.tsv",EVO_IMMUNO_POP))
names(miRNAs)[1] = "miRNA"
miRNAs = miRNAs[, .(miRNA)]
#miRNAs coordinates
miRNA_coordinate = fread(paste(EVO_IMMUNO_POP, "ERCPilot_SharedDBs/mirbase20/miRNA_mature_coordinates_strandinfo.bed", sep=""))
names(miRNA_coordinate) = c("chromosome", "start", "end", "mir1", "V5", "strand", "mir2", "hairpin")

#We compute miRNA that are duplicated (have several sets of coordinates)
duplicated_coordinates = miRNA_coordinate[, .(number_coordinates = nrow(unique(.SD[, .(chromosome, start, end)]))), by = c("mir1")]

miRNAs = miRNAs[miRNA %in% duplicated_coordinates[number_coordinates == 1, mir1]]


miRNAs[, chromosome := miRNA_coordinate[match(miRNAs$miRNA, mir1), chromosome]]
miRNAs[, start := miRNA_coordinate[match(miRNAs$miRNA, mir1), start]]
miRNAs[, end := miRNA_coordinate[match(miRNAs$miRNA, mir1), end]]
miRNAs[, strand := miRNA_coordinate[match(miRNAs$miRNA, mir1), strand]]

#######################################
##We had the gene data on top of that##
#######################################
genes_coordinates = fread(paste(EVO_IMMUNO_POP, "Maxime/Evo_Immuno_pop_data/GeneAnnotation_hg37_ens70.txt", sep=""))
genes_coordinates[get("Gene Biotype") != "miRNA"]
genes_coordinates_GR = GRanges(paste("chr", genes_coordinates$Chromosome.Name, sep=""), IRanges(start = genes_coordinates[, get("Gene Start (bp)")], end = genes_coordinates[, get("Gene End (bp)")]),  strand = ifelse(genes_coordinates$Strand == 1, "+", "-"))

miRNAs[chromosome == "chrX", chromosome := "chr23"]
miRNAs_GR = GRanges(miRNAs$chromosome, IRanges(start = miRNAs$start, end = miRNAs$end), strand = miRNAs$strand)

ov = as.data.table(findOverlaps(miRNAs_GR, genes_coordinates_GR))

miRNAs[, in_a_gene := 1:.N %in% ov$queryHits]

Genic_miRNAs = cbind(miRNAs[ov$queryHits], genes_coordinates[ov$subjectHits])

Genic_miRNAs[, intersection_start := pmax(start, get("Gene Start (bp)"))]
Genic_miRNAs[, intersection_end  := pmin(end, get("Gene End (bp)"))]
Genic_miRNAs[(intersection_start == start) & (intersection_end == end) ]
Genic_miRNAs_expressed = Genic_miRNAs[Expressed == T]
Genic_miRNAs_expressed = Genic_miRNAs_expressed[, .(ensemblGenes = paste(Ensembl.Gene.ID, collapse = "//"),
                                                    geneNames = paste(get("Associated Gene Name"), collapse = "//")), by = miRNA]

miRNAs[, in_an_expressed_gene := miRNA %in% Genic_miRNAs_expressed$miRNA]
miRNAs[, ensemblGenes := Genic_miRNAs_expressed[match(miRNAs$miRNA, miRNA), ensemblGenes]]
miRNAs[, geneNames := Genic_miRNAs_expressed[match(miRNAs$miRNA, miRNA), geneNames]]
write.table(miRNAs, sprintf("%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_intersection_with_genes.tsv",EVO_IMMUNO_POP), sep="\t", quote = F, row.names = F)
