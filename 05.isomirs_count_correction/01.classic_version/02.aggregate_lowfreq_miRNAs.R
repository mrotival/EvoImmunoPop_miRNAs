
###########
##Imports##
###########
suppressMessages(require("reshape"))
suppressMessages(require("DESeq2"))
suppressMessages(require(ggplot2))
suppressMessages(require("sva"))
library(stringr)

isomir_counts_file_norm=fread(sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/isomiR_counts.RPM.raw_23447_V2.0_MR.tsv",EVO_IMMUNO_POP))
isomir_annot=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_FULL_V2.0.tsv',EVO_IMMUNO_POP))

samples_names=unlist(fread(sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/sample_names_977_highQuality.tsv",EVO_IMMUNO_POP)))

all(isomir_annot$ID==isomir_counts_file_norm$isomir_ID)
isomir_counts_file_norm=merge(isomir_counts_file_norm,isomir_annot,by.x='isomir_ID',by.y='ID')

Data=str_split_fixed(isomir_counts_file_norm$isomir_ID,'_',11)
# add annotations
isomir_counts_file_norm$chrom=Data[,1]
isomir_counts_file_norm$strand=Data[,4]
isomir_counts_file_norm$subs=Data[,5]

############################
## melt isomiR data       ##
############################

isomir_counts_melted=melt(isomir_counts_file_norm, id.vars=c("isomir_ID","mirID","isomir_sequence","shift_5p","shift_3p","subs","newID","hsa_ID","MIMAT","MI_ID","miRNA_arm","mean_isoMiR","ratio_isoMiR",'chrom','strand'), measure.vars=samples_names,  variable.name = "library_ID", value.name = "counts")

# aggregate all reads from the same isomiR
isomir_grouped_counts_melted=isomir_counts_melted[,.(total_counts=sum(counts)),by=.(newID,library_ID)]

## keep only the info from the main isomiR
isomir_counts_melted = isomir_counts_melted[order(mirID, -mean_isoMiR, library_ID),]
isomir_counts_melted = unique(isomir_counts_melted, by =c('newID','library_ID'))

colnames(isomir_grouped_counts_melted)=c("isomir_ID","library_ID","count")
colnames(isomir_counts_melted)=c("dominant_sequence_ID","mirID","isomir_sequence","start","end","subs","isomir_ID","hsa_ID","MIMAT","MI_ID","miRNA_arm","mean_dominant_sequence","ratio_mean_dominant_sequence","chromosome","strand","library_ID","count_dominant_sequence")
isomir_grouped_counts_melted=merge(isomir_counts_melted,isomir_grouped_counts_melted,by=c('isomir_ID','library_ID'))

mir_grouped_counts_melted=isomir_grouped_counts_melted[,.(total_miR_counts=sum(count)),by=.(mirID,library_ID)]
isomir_grouped_counts_melted=merge(isomir_grouped_counts_melted,mir_grouped_counts_melted,by=c('mirID','library_ID'))
isomir_grouped_counts_melted[,ratio:=count/total_miR_counts]

fwrite(isomir_grouped_counts_melted,file=paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/05.isomirs_count_correction/miRNA_isomiR_raw_counts_aggregated_melded.tsv", sep=""),sep='\t')

first <- function(x){
  if (is.na(x[1])){
    return(NA)
  }else{
    return(x[1])
  }
}
raw_isomiR_counts = dcast(isomir_grouped_counts_melted,chromosome + start + end + strand + subs + hsa_ID + MIMAT + MI_ID + isomir_ID ~ library_ID, value.var = "count", first, fill=0)
raw_isomiR_ratios = dcast(isomir_grouped_counts_melted,chromosome + start + end + strand + subs + hsa_ID + MIMAT + MI_ID + isomir_ID ~ library_ID, value.var = "ratio", first, fill=0)

fwrite(raw_isomiR_counts,file=paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/05.isomirs_count_correction/miRNA_isomiR_raw_counts_aggregated_cast.tsv", sep=""),sep='\t')
fwrite(raw_isomiR_ratios,file=paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/05.isomirs_count_correction/miRNA_isomiR_raw_ratios_aggregated_cast.tsv", sep=""),sep='\t')
