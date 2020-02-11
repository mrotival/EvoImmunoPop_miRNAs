
#################################################################################
######################        load libraries            #########################
#################################################################################

require(ggplot2)



########################################################################################
######################             load the data               #########################
########################################################################################

miRqtl_cis_annot= fread(sprintf("%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable4A_mirQTLs_Annoted.tsv",EVO_IMMUNO_POP))
miRqtl_cis_annot[,conserved_TSS:=Conservation>.2]

mir_TSS_annot = fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_TSS_annot_DeRie2017.txt',EVO_IMMUNO_POP))
mir_TSS_annot[,has_mirQTL:=hsa_ID%in%miRqtl_cis_annot$miRNA]
mir_TSS_annot[,conserved_TSS:=Conservation>.2]

 length(mir_TSS_annot[Nb_loci==1, unique(hsa_ID)])
 
########################################################################################
################# depletion of miRQTL among conserved TSS miRNA ########################
########################################################################################


odds.ratio <- function(tab, alpha = 0.05){
	test=fisher.test(tab,conf.level=1-alpha)
	oframe = data.frame(LowerCI = test$conf.int[1], OR = test$est, UpperCI = test$conf.int[2], alpha = alpha, P = test$p.value)
	oframe
}

table(ifelse(mir_TSS_annot$has_mirQTL,'mirQTL','no mirQTL'),ifelse(mir_TSS_annot$Conservation>.2,'conserved','not conserved'))
#            conserved not conserved
#  mirQTL           45            46
#  no mirQTL       286           157

odds.ratio(table(mir_TSS_annot$has_mirQTL,mir_TSS_annot$Conservation>.2))
#             LowerCI        OR   UpperCI alpha           P
#odds ratio 0.3321654 0.5376576 0.8694213  0.05 0.008888615

########################################################################################
###################### Dist of mirQTL from TSS VS conservation #########################
########################################################################################

pdf(sprintf('%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/Dist_of_mirQTL_from_TSS.pdf',EVO_IMMUNO_POP),width=4,height=2)
colConserv=c('FALSE'="#969696AA",'TRUE'="#E31A1CAA")
p <- ggplot(miRqtl_cis_annot,aes(CisDist_TSS,col=conserved_TSS,fill=conserved_TSS)) + geom_density(alpha=0.4,bw=0.025)
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p <- p + scale_colour_manual(values=colConserv) + scale_fill_manual(values=colConserv)
p <- p + xlab('Distance to TSS') + ylab('Density')
print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/Dist_of_mirQTL_from_TSS_absolute.pdf',EVO_IMMUNO_POP),width=4,height=2)
p <- ggplot(miRqtl_cis_annot,aes(abs(CisDist_TSS),col=conserved_TSS,fill=conserved_TSS)) + geom_density(alpha=0.4,bw=0.025)
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p <- p + scale_colour_manual(values=colConserv) + scale_fill_manual(values=colConserv)
p <- p + xlab('|Distance to TSS|') + ylab('Density')
print(p)
dev.off()

wilcox.test(abs(miRqtl_cis_annot$CisDist_TSS)~miRqtl_cis_annot$conserved_TSS)
#W = 742, p-value = 0.02024

miRqtl_cis_annot[,mean(abs(CisDist_TSS)),by=conserved_TSS]
#          TRUE 0.12503202
#         FALSE 0.09044772
# + 3.5 kb 
#########################################################################################
###########################  frequency of mirQTL mechanisms ############################
########################################################################################

tab_mirQTL_type=table(miRqtl_cis_annot$typeDetail[miRqtl_cis_annot$typeDetail!=''])
tab_mirQTL_type


#                          Distant Hairpin-altering Hairpin-flanking   mirNA-altering     TSS-flanking 
#              31               50                7               17                3               14 

colorRegion=c("Hairpin-flanking"="#FFD92F","Hairpin-altering"="#FC8D62","miRNA-altering"="#E78AC3","TSS-flanking"="#66C2A5","Distant"="#8DA0CB", "Undecided"="#888888")

pdf(sprintf("%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/pie_miRNA_QTL_Type.pdf",EVO_IMMUNO_POP),width=4,height=4)
par(mar=c(3,3,3,3))
tabpct=paste(round(100*tab_mirQTL_type/sum(tab_mirQTL_type), 1),'%')
pie(tab_mirQTL_type,col=colorRegion[names(tab_mirQTL_type)],init.angle=90,labels=tabpct) # 5.5 x 5.5 inches
dev.off()


############################################################################################################
########################### Dist of mirQTL from TSS VS Dist of mirQTL from miRNA ###########################
############################################################################################################


asinh_trans=function() 
{
    trans <- function(x) asinh(x)
    inv <- function(x) sinh(x)
    scales:::trans_new(paste0("asinh"), trans, inv, scales::extended_breaks(),
        domain = c(-Inf, Inf))
}

library(ggplot2)
pdf(sprintf('%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/Dist_of_mirQTL_from_TSS_and_miRNA_byType.pdf',EVO_IMMUNO_POP),width=3.6,height=2)
p <- ggplot(miRqtl_cis_annot,aes(CisDist_miRNA*1000, CisDist_TSS*1000,col=typeDetail))+geom_point(alpha=.6) 
p <- p + theme_bw()+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p <- p + scale_colour_manual(values=colorRegion)
p <- p + xlab('Distance to miRNA (kb)') + ylab('Distance to TSS (kb)') + xlim(-600,600) + ylim(-600,600)
print(p)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/11.response_miRQTL_and_other_analyses/Dist_of_mirQTL_from_TSS_and_miRNA_byType_rescaled.pdf',EVO_IMMUNO_POP),width=2,height=2)
#p <- p + scale_y_continuous(trans = asinh_trans(),breaks=c(-1000,-100,-10,-1,0,1,10,100,1000)) + scale_x_continuous(trans = asinh_trans(),breaks=c(-1000,-100,-10,-1,0,1,10,100,1000))
p <- p + scale_y_continuous(trans = asinh_trans(),breaks=c(-50,-5,0,5,50),limits=c(-50,50))+ scale_x_continuous(trans = asinh_trans(),breaks=c(-50,-5,0,5,50),limits=c(-50,50))
p <- p + xlab('') + ylab('') + theme(legend.position = "none")
print(p) 
dev.off()



############################################################################################################
########################### enrichment of mirQTLs in specific classess of SNPs  ############################
############################################################################################################


#### add here:  enrichment in Enhancer/promoters defined by ChromHMM
#### add here:  enrichment in known TFBSs (Encode)

# for all miRNAs (122)
# for TSS miRNAs (12)
# for Distant miRNAs (42)
# for hairpin miRNAs (26)

snps=fread(sprintf("%s/Maxime/miRNA_V2/data/08.snp_data/general_informations/all_snps_MAFover5.tsv",EVO_IMMUNO_POP))
snps$mirQTL= snps$snp.name%in% miRqtl_cis_annot$snps
OR=list()
for(i in unique(snps[,mir_cat])){
	if(i!=''){
		OR[[i]]=odds.ratio(table(snps$mirQTL, snps$mir_cat==i))
		OR[[i]]$mir_cat=i
	}	
}
OR=rbindlist(OR)
fwrite(OR,sprintf("%s/Maxime/miRNA_V2/data//11.response_miRQTL_and_other_analyses/miRNA_location.tsv",EVO_IMMUNO_POP),sep='\t')

Map_imputed=fread('/Volumes/evo_immuno_pop/Maxime/SNP_annotation/Map_imputed_essential_informations.txt')
Map_imputed[,daf_EU:=as.numeric(gsub('\\((.*)\\)','\\1',daf_char_EUB))]
Map_imputed[,daf_AF:=as.numeric(gsub('\\((.*)\\)','\\1',daf_char_AFB))]
Map_imputed[,maf_EUB:=pmin(daf_EU,1-daf_EU)]
Map_imputed[,maf_AFB:=pmin(daf_AF,1-daf_AF)]
Map_imputed[,sum(maf_EUB>=0.05 | maf_AFB>=0.05)]
#10261270
Map_imputed=Map_imputed[which(maf_EUB>=0.05 | maf_AFB>=0.05),]
quantile(abs(Map_imputed$iHS_EUB[which(maf_EUB>=0.05]),0.99,na.rm=T) # 2.56
mean(abs(Map_imputed$iHS_EUB[which(maf_EUB>=0.05])>2.5,na.rm=T) # 0.011



require(GenomicRanges)

snps_GR=makeGRangesFromDataFrame(Map_imputed,start.field='position',end.field='position')
miR_GR=makeGRangesFromDataFrame(mir_TSS_annot[Nb_loci==1, ],seqnames='mat_chr',start.field='mat_start',end.field='mat_end',strand.field='mat_strand')
seqlevels(miR_GR)=gsub('chr','', seqlevels(miR_GR))
ooCis=findOverlaps(snps_GR,union(flank(miR_GR,1e6,both=T),flank(miR_GR,1e6,both=T,start=F)))
Map_imputed$Cis_1Mb_miR=FALSE
Map_imputed$Cis_1Mb_miR[unique(from(ooCis))]=TRUE
sum(Map_imputed$Cis_1Mb_miR)

