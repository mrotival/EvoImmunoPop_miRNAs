######################################
##Getting reads in a workable format##
######################################
# raw_isomiR_counts_melted = fread(paste(EVO_IMMUNO_POP, "Martin/miRNA/09.isomirs/data/miRNA_transcripts_raw_counts_melded.tsv", sep=""))
# raw_isomiR_counts_melted = raw_isomiR_counts_melted[, .(count = sum(count)), by = c("chromosome", "start", "end", "strand","subsitutions", "mir1", "mir2","hairpin",  "library_ID")]
# raw_isomiR_counts_melted[, isomir_ID := paste(chromosome , start , end , strand ,  mir1 , mir2 , hairpin, sep="_")]

library(stringr)


#############################
## get Raw Read Count data ##
#############################
isomir_counts_file_norm=fread(sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/isomiR_counts.RPM.raw_23447_V2.0_MR.tsv",EVO_IMMUNO_POP))
isomir_annot=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_FULL_V2.0.tsv',EVO_IMMUNO_POP))
samples_names=unlist(fread(sprintf("%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/sample_names_977_highQuality.tsv",EVO_IMMUNO_POP)))

# merge tables
all(isomir_annot$ID==isomir_counts_file_norm$isomir_ID)
isomir_counts_file_norm=merge(isomir_counts_file_norm,isomir_annot,by.x='isomir_ID',by.y='ID')

############################
## annotate substitutions ##
############################

Data=str_split_fixed(isomir_counts_file_norm$isomir_ID,'_',11)
# add annotations
isomir_counts_file_norm$chrom=Data[,1]
isomir_counts_file_norm$strand=Data[,4]
# extract substitutions 
isomiR_subs=substr(Data[,6],3,100)
# 1 line per substitution
subs=str_match_all(isomir_counts_file_norm$isomir_ID,'([0-9]+)([ATGC])->([ATGC])')
# get position and type of substituions
subs_position=unlist(sapply(subs,function(x){as.numeric(x[,2])}))
subs_position_3p=rep(nchar(Data[,5]),sapply(subs,nrow))-subs_position+1
subs_type=unlist(sapply(subs,function(x){paste(x[,3],x[,4],sep='->')}))
# get corresponding miR
subs_miR_id=rep(isomir_counts_file_norm$isomir_ID,sapply(subs,nrow))

# identify terminal uridylation and adenylation
is_termunalUorA=substr(subs_type,4,4)%in%c('A','T') & subs_position_3p<4

# create a variable with only these susbtitutions
subs_id=by(paste(subs_position[is_termunalUorA],subs_type[is_termunalUorA],sep=''),subs_miR_id[is_termunalUorA],paste,collapse=';')
subs_id=paste(subs_id,';',sep='')[match(isomir_counts_file_norm$isomir_ID,names(subs_id))]
subs_id[is.na(subs_id)]=''


# change isomiR sequence to contain only these substitution
isomir_withSubs=which(sapply(subs,length)>0)

tic=Sys.time()
isomir_sequence_nosubs=isomir_counts_file_norm$isomir_sequence
mySeq=strsplit(isomir_counts_file_norm$isomir_sequence,"")
for( i in isomir_withSubs){
    mysubs=subs[[i]]
    Seq=mySeq[[i]]
    for (j in 1:nrow(mysubs)){
        pos=as.numeric(mysubs[j,2])
        from=mysubs[j,3]
        to= mysubs[j,4]
        if((length(Seq)-pos+1)>3 | !(to%in%c('A','T'))){
            Seq[pos]=from
        }
    }
    Seq=paste(Seq,collapse='')
    isomir_sequence_nosubs[i]=Seq
}
print(Sys.time()-tic)

# add these new info into isomir_counts_file_norm
nsubs=sapply(strsplit(subs_id,';'),length)

isomir_counts_file_norm$subs_id=paste(nsubs,subs_id,sep=';')
isomir_counts_file_norm$isomir_sequence_nosubs=isomir_sequence_nosubs


# create a new isomiR ID
isomir_counts_file_norm[,newID_nosubs:=ifelse(mean_isoMiR<1 | ratio_isoMiR<1,paste(mirID,'other',sep='_'),paste(chrom,shift_5p,shift_3p,strand,isomir_sequence_nosubs,subs_id,mirID,sep='_'))]

############################
## melt isomiR data       ##
############################
isomir_counts_melted=melt(isomir_counts_file_norm, id.vars=c("isomir_ID","mirID","isomir_sequence","shift_5p","shift_3p","subs_id","newID_nosubs","hsa_ID","MIMAT","MI_ID","miRNA_arm","mean_isoMiR","ratio_isoMiR",'chrom','strand'), measure.vars=samples_names,  variable.name = "library_ID", value.name = "counts")

# aggregate all reads from the same isomiR
isomir_grouped_counts_melted=isomir_counts_melted[,.(total_counts=sum(counts)),by=.(newID_nosubs,library_ID)]

## keep only the info from the main isomiR
isomir_counts_melted = isomir_counts_melted[order(mirID, -mean_isoMiR, library_ID),]
isomir_counts_melted = unique(isomir_counts_melted, by =c('newID_nosubs','library_ID'))

colnames(isomir_grouped_counts_melted)=c("isomir_ID","library_ID","count")
colnames(isomir_counts_melted)=c("dominant_sequence_ID","mirID","isomir_sequence","start","end","subs","isomir_ID","hsa_ID","MIMAT","MI_ID","miRNA_arm","mean_dominant_sequence","ratio_mean_dominant_sequence","chromosome","strand","library_ID","count_dominant_sequence")
isomir_grouped_counts_melted=merge(isomir_counts_melted,isomir_grouped_counts_melted,by=c('isomir_ID','library_ID'))

mir_grouped_counts_melted=isomir_grouped_counts_melted[,.(total_miR_counts=sum(count)),by=.(mirID,library_ID)]
isomir_grouped_counts_melted=merge(isomir_grouped_counts_melted,mir_grouped_counts_melted,by=c('mirID','library_ID'))
isomir_grouped_counts_melted[,ratio:=count/total_miR_counts]

fwrite(isomir_grouped_counts_melted,file=paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/05.isomirs_count_correction/miRNA_isomiR_raw_counts_aggregated_nosubs_melded.tsv", sep=""),sep='\t')
# isomir_grouped_counts_melted=fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/05.isomirs_count_correction/miRNA_isomiR_raw_counts_aggregated_nosubs_melded.tsv", sep=""))

first <- function(x){
  if (is.na(x[1])){
    return(NA)
  }else{
    return(x[1])
  }
}
raw_isomiR_counts = dcast(isomir_grouped_counts_melted,chromosome + start + end + subs +strand + hsa_ID + MIMAT + MI_ID + isomir_ID ~ library_ID, value.var = "count", first, fill=0)
raw_isomiR_ratios = dcast(isomir_grouped_counts_melted,chromosome + start + end + subs +strand + hsa_ID + MIMAT + MI_ID + isomir_ID ~ library_ID, value.var = "ratio", first, fill=0)

fwrite(raw_isomiR_counts,file=paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/05.isomirs_count_correction/miRNA_isomiR_raw_counts_aggregated_nosubs_cast.tsv", sep=""),sep='\t')
fwrite(raw_isomiR_ratios,file=paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/05.isomirs_count_correction/miRNA_isomiR_raw_ratios_aggregated_nosubs_cast.tsv", sep=""),sep='\t')


##we keep the notation the most diverse possible
#isomirs_informations = unique(raw_isomiR_counts_melted[, .(chromosome , start , end , terminal_substitutions, strand ,  mir1 , mir2 , hairpin,isomir_ID)])
#

