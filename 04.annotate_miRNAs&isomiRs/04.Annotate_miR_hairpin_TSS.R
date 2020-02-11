
###################################################
####  annotate miRNA position based on miRbase ####
###################################################

mirbase20_mature_coords = fread(sprintf("%s/ERCPilot_SharedDBs/mirbase20/miRNA_mature_coordinates_strandinfo.bed",EVO_IMMUNO_POP))
colnames(mirbase20_mature_coords)=c('mat_chr','mat_start','mat_end','hsa_ID','dot','mat_strand','MIMAT','MI_ID')

mirbase20_hairpin_coords = fread(sprintf("%s/ERCPilot_SharedDBs/mirbase20/miRNA_hairpin_coordinates_strandinfo.bed",EVO_IMMUNO_POP))
colnames(mirbase20_hairpin_coords)=c('hairpin_chr','hairpin_start','hairpin_end','hairpin_name','dot','hairpin_strand','MI_ID')

miRNA_coords=merge(mirbase20_mature_coords, mirbase20_hairpin_coords,all=TRUE)
miRNA_coords[,Nb_loci:=length(unique(MI_ID)),by=hsa_ID]

###################################################
####  annotate miRNA TSS based on Fantom data  ####
###################################################

require(splitstackshape)
miRNA_TSS=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_TSS_DeRie_natureBioTech_2017_TableS15.txt',EVO_IMMUNO_POP))
colnames(miRNA_TSS)=make.names(colnames(miRNA_TSS))
library(splitstackshape)
miRNA_TSS=cSplit(miRNA_TSS, splitCols = "Pre.microRNAs", sep = ",", direction = "long", drop = FALSE)
#miRNA_TSS$Pre.microRNAs_simplified=sapply(strsplit(miRNA_TSS[,as.character(Pre.microRNAs)],'-'),function(x){paste(x[1:3],collapse='-')})

mirID_658_expressed = fread(sprintf('%s/Maxime/miRNA_V2/data/03b.isomirs_alignment/mirID_658_expressed.tsv',EVO_IMMUNO_POP))[[1]]
#table(sapply(strsplit(miRNA_TSS[Pre.microRNAs,'-'),function(x){x[5]}))
mir_TSS_annot=miRNA_coords[hsa_ID%in%mirID_658_expressed,]
miR_TSS_expressed=miRNA_TSS[Pre.microRNAs%in%mir_TSS_annot$hairpin_name,]
#miR_TSS_expressed=miRNA_TSS[Pre.microRNAs_simplified%in%pre_mirID_658_expressed,]

cols = c('Chromosome','Strand','TSS','Promoter.name','Primary.microRNA','Intronic.intergenic','Pre.microRNAs','Conservation','Maximum.expression.level..tpm.','Highest.expressing.sample')
mir_TSS_annot=merge(mir_TSS_annot,miR_TSS_expressed[,mget(cols)],by.x='hairpin_name',all.x=TRUE,by.y='Pre.microRNAs')
mir_TSS_annot[,Nb_TSS:=sum(!is.na(TSS)),by= hsa_ID]
luq(mir_TSS_annot[Nb_TSS==1, hsa_ID]) # 476 miRNA with a single annotated TSS
mir_TSS_annot[,arm:=ifelse(grepl('3p',hsa_ID),'3p',ifelse(grepl('5p',hsa_ID),'5p',NA))]
mir_TSS_annot[,DistToHairpinStart:=ifelse( hairpin_strand=='+',mat_start-hairpin_start,hairpin_end-mat_end)]
mir_TSS_annot[,DistToHairpinEnd:=ifelse( hairpin_strand=='+',hairpin_end-mat_end,mat_start-hairpin_start)]
mir_TSS_annot[,assigned_arm:=ifelse(DistToHairpinStart<DistToHairpinEnd,"5p","3p")]

fwrite(mir_TSS_annot,file=sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_TSS_annot_DeRie2017.txt',EVO_IMMUNO_POP),sep='\t')

pdf(sprintf('%s/Maxime/miRNA_V2/figures/04.annotate_miRNAs&isomiRs/Dist_from_hairpin_end_VS_arm.pdf',EVO_IMMUNO_POP),width=4,height=2.5)
p <- ggplot(mir_TSS_annot,aes(DistToHairpinEnd,col=arm,fill=arm)) + geom_density(alpha=0.4)
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p <- p + scale_colour_manual(values=brewer.pal(3,'Set3')[-2],na.value='darkgrey') + scale_fill_manual(values=brewer.pal(3,'Set3')[-2],na.value=grey(.6))
p + xlab('Distance to hairpin end') + ylab('Density')
print(p)
dev.off()    

pdf(sprintf('%s/Maxime/miRNA_V2/figures/04.annotate_miRNAs&isomiRs/Dist_from_hairpin_start_VS_arm.pdf',EVO_IMMUNO_POP),width=4,height=2.5)
p <- ggplot(mir_TSS_annot,aes(DistToHairpinStart,col=arm,fill=arm)) + geom_density(alpha=0.4)
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p <- p + scale_colour_manual(values=brewer.pal(3,'Set3')[-2],na.value='darkgrey') + scale_fill_manual(values=brewer.pal(3,'Set3')[-2],na.value=grey(.6))
p + xlab('Distance to hairpin start') + ylab('Density')
print(p)
dev.off()    

pdf(sprintf('%s/Maxime/miRNA_V2/figures/04.annotate_miRNAs&isomiRs/Dist_from_hairpin_limits_VS_arm.pdf',EVO_IMMUNO_POP),width=4,height=3)
p <- ggplot(mir_TSS_annot,aes(x=DistToHairpinStart,DistToHairpinEnd,col=arm)) + geom_point(alpha=.8)
    p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
    p <- p + scale_color_manual(values=brewer.pal(3,'Set3')[-2],na.value='darkgrey')
    p <- p + geom_abline(aes(intercept=0, slope=1)) + xlab('Distance to hairpin start') + ylab('Distance to hairpin end') +xlim(0,100)+ylim(0,100)
print(p)
dev.off()    

###################################################
####   merge all annotations with Cis mirQTLs  ####
###################################################
hairpin_annot=dcast(mir_TSS_annot, hairpin_name + MI_ID + hairpin_chr + hairpin_strand + TSS + Conservation + hairpin_start + hairpin_end~ assigned_arm, value.var=c('hsa_ID','mat_start','mat_end'))
fwrite(hairpin_annot,file=sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_hairpin_annot.txt',EVO_IMMUNO_POP),sep='\t')
