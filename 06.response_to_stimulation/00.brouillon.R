##############
##librairies##
##############
suppressMessages(require(DESeq2))

################
##Loading data##
################
corrected_count = read.delim(paste(EVO_IMMUNO_POP, "Martin/miRNA/03.alignment_and_count_correction/data/miRNA_counts.RPM.GC_Batch_corrected.tsv", sep=""), sep="\t")
colnames(corrected_count) = gsub("\\.", "-", colnames(corrected_count))


colData = data.table(samples = colnames(corrected_count))
colData[, pop := as.factor(substr(samples, 1, 3))]
colData[, condition := as.factor(substr(samples, 8,8), levels = 1:5)]
colData[, individual :=  as.factor(substr(samples, 1,6))]
rownames(colData) = colData$samples



design <- ~ pop + condition
dds = DESeqDataSetFromMatrix(countData = round(corrected_count), colData = colData, design = design)
sizeFactors(dds)<-1

#############
##Brouillon##
#############
dds = dds[,dds$condition==5| dds$condition==1]
dds$condition<-droplevels(dds$condition,except=unique(dds$condition))
dds$individual<-droplevels(dds$individual,except=unique(dds$individual))

dds = estimateDispersions(dds) ##Carefull, it's pretty long
dds = nbinomWaldTest(dds)
pouloup = results(dds)
