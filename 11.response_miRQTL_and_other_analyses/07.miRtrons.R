#miRtrons.R
condIndex=c("NS", "LPS", "PAM3CSK4", "R848", "IAV")

miR_QTL=fread(sprintf('%s/Maxime/miRNA_V2/data/09.mirQTL/cis_mirQTLs_with_FDR_filtered_best_miRNA_snp_association.tsv',EVO_IMMUNO_POP))
miRNA_gene_overlap=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_intersection_with_genes.tsv',EVO_IMMUNO_POP))

miRtron_QTL=merge(miRNA_gene_overlap[which(in_an_expressed_gene),],miR_QTL,by.x='miRNA',by.y='gene')

luq(miRtron_QTL$miRNA)
#[1] 64
luq(miRNA_gene_overlap[which(in_an_expressed_gene),miRNA])
#352
	
tab=table(ifelse(miRNA_gene_overlap$in_an_expressed_gene,'1_miRtron','0_other'),ifelse(miRNA_gene_overlap$miRNA%in%miR_QTL$gene,'1_hasQTL','0_noQTL'))
# odds.ratio(tab)
#              LowerCI        OR  UpperCI alpha        P
# odds ratio 0.4735984 0.7207251 1.097892  0.05 0.121939


eQTL=fread(sprintf('%s/Maxime/evo_immuno_pop_QTLs/eQTL/ALL/eQTL_ALL_bestSNP_FDR5_5cond.txt',EVO_IMMUNO_POP))

trQTL=fread(sprintf('%s/Maxime/evo_immuno_pop_QTLs/Transcription_rate_QTL/Transcription_rate_QTL_bestSNP_FDR5_5cond.txt',EVO_IMMUNO_POP)) 
miRtron_QTL[,condition:=condIndex[condition]]


miRtron_eQTL=merge(eQTL,miRtron_QTL,by.x=c('gene','condition'),by.y=c('ensemblGenes','condition'),suffix=c('.gene','.miRNA'))
miRtron_trQTL=merge(trQTL,miRtron_QTL,by.x=c('gene','condition'),by.y=c('ensemblGenes','condition'),suffix=c('.gene','.miRNA'))

luq(miRtron_eQTL$miRNA) # 26
luq(miRtron_trQTL$miRNA)# 18

length(intersect(miRtron_trQTL$miRNA,miRtron_eQTL$miRNA)) # 13
length(union(miRtron_trQTL$miRNA,miRtron_eQTL$miRNA)) # 31


miRtron_eQTL_bestSNP=miRtron_eQTL[order(miRtron_eQTL$miRNA,pvalue.miRNA),]

LD_pair=mapply(function(snp.gene,snp.miRNA){
								SNP1=getSNP(snp.gene);
								SNP2=getSNP(snp.miRNA)
								cor(SNP1,SNP2,use='p')^2},miRtron_eQTL_bestSNP[,snps.gene],miRtron_eQTL_bestSNP[,snps.miRNA])
							
miRtron_eQTL_bestSNP[,LD_pair:=LD_pair]								
LD_miRtron_eQTL_bestSNP=miRtron_eQTL_bestSNP[,.(snps.miRNA,snps.gene,miRNA,geneNames,gene,pvalue.miRNA,pvalue.gene,condition,LD_pair)]

miR_QTL_table=fread(sprintf('%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable3A_mirQTLs_Annoted_V2.tsv',EVO_IMMUNO_POP))

miR_QTL_table[,r2_with_eQTL:=0]
miR_QTL_table[!miRNA %chin% miRtron_QTL$miRNA,r2_with_eQTL:=NA]
miR_QTL_table[snps%in% miRtron_eQTL_bestSNP$snps.miRNA,r2_with_eQTL:=LD_pair[match(snps,miRtron_eQTL_bestSNP[,snps.miRNA])]]


########### same for TrQTLs

miRtron_trQTL_bestSNP=miRtron_trQTL[order(miRtron_trQTL$miRNA,pvalue.miRNA),]

LD_pair=mapply(function(snp.gene,snp.miRNA){
								SNP1=getSNP(snp.gene);
								SNP2=getSNP(snp.miRNA)
								cor(SNP1,SNP2,use='p')^2},miRtron_trQTL_bestSNP[,snps.gene],miRtron_trQTL_bestSNP[,snps.miRNA])
							
miRtron_trQTL_bestSNP[,LD_pair:=LD_pair]								



# annotate miR_QTLs


miR_QTL_table[,r2_with_trQTL:=0]
miR_QTL_table[!miRNA %chin% miRtron_QTL$miRNA,r2_with_trQTL:=NA]
miR_QTL_table[snps%in% miRtron_trQTL_bestSNP$snps.miRNA,r2_with_trQTL:=LD_pair[match(snps,miRtron_trQTL_bestSNP[,snps.miRNA])]]


fwrite(miR_QTL_table,sprintf('%s/Maxime/miRNA_V2/data/00_tables_publication/SupTable3A_mirQTLs_Annoted_V2.1_withR2_eQTL.tsv',EVO_IMMUNO_POP))


###########################################
##############@ Load data #################
###########################################

# mRNA expression
genes_expression = fread(paste(EVO_IMMUNO_POP, "Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/FPKM_matrix.txt", sep=""))
genes_expression = melt(genes_expression, id = "ID", variable.name = "sample", value.name = "expression")
genes_expression[, individual := substr(sample, 1,6)]
genes_expression[, condition := condIndex[as.numeric(substr(sample, 8,8))]]
genes_expression[, population := substr(sample, 1,3)]
names(genes_expression)[1] = "gene"

# miRNAs
miRNA_expression= fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_counts.log2RPM.GCRL_Batch_corrected_V2.0_MR.tsv", sep=""))
names(miRNA_expression)[1] = "miRNA"
miRNA_expression = melt(miRNA_expression, id = "miRNA", variable.name = "sample", value.name = "expression")
miRNA_expression[, individual := substr(sample, 1,6)]
miRNA_expression[, condition := condIndex[as.numeric(substr(sample, 8,8))]]
miRNA_expression[, population := substr(sample, 1,3)]
miRNA_expression[, sample:=substr(sample,1,8)]
# SNP genotype
SNP_genos=as.data.table(getGenos(unique(c(LD_miRtron_eQTL_bestSNP$snps.miRNA,LD_miRtron_eQTL_bestSNP$snps.gene))))
SNP_genos=melt(SNP_genos, id.vars=c("snp.name","chromosome","position","locus","posID"),variable.name='individual',value.name='genotype')
SNP_genos[,genotype:=2-genotype]

Adjust_miRNA_QTL=function(miR,geneID,snp,cond){
	miR_gene_expression=merge(miRNA_expression[miRNA==miR & condition==cond,],
							  genes_expression[gene==geneID & condition==cond,],
				 			  by=c('sample','individual','condition','population'),
					 		  suffix=c('.miRNA','.gene'))
	data_to_use = merge(miR_gene_expression,SNP_genos[snp.name==snp,],by='individual')
	mod_miR_snp=lm(expression.miRNA~genotype,data=data_to_use)
	mod_miR_snp_aj=lm(expression.miRNA~expression.gene+genotype,data=data_to_use)
	mod_miR_gene=lm(expression.miRNA~expression.gene,data=data_to_use)
	#
	mod_gene_snp=lm(expression.gene~genotype,data=data_to_use)
	mod_gene_snp_aj=lm(expression.gene~expression.miRNA+genotype,data=data_to_use)
	mod_gene_miR=lm(expression.gene~expression.miRNA,data=data_to_use)
	#
	L_miR_causal=logLik(mod_gene_miR)+logLik(mod_miR_snp)
	L_gene_causal=logLik(mod_miR_gene)+logLik(mod_gene_snp)
	L_Independant=logLik(mod_gene_snp)+logLik(mod_miR_snp)
	#
	P_miR_snp_adj=summary(mod_miR_snp_aj)$coef[3,4]
	P_gene_snp_adj=summary(mod_gene_snp_aj)$coef[3,4]
	P_gene_miR=summary(mod_miR_gene)$coef[2,4]
	P_gene_miR_adj=summary(mod_miR_snp_aj)$coef[2,4]
	data.table(miR,geneID,snp,cond,L_miR_causal,L_gene_causal,L_Independant,P_gene_snp_adj,P_miR_snp_adj,P_gene_miR, P_gene_miR_adj)
	}


res1=list()
res2=list()
for (i in 1:nrow(LD_miRtron_eQTL_bestSNP)){
	cat(i,'')
	snp.miRNA=LD_miRtron_eQTL_bestSNP[i,snps.miRNA]
	snp.gene=LD_miRtron_eQTL_bestSNP[i,snps.gene]
	cond=LD_miRtron_eQTL_bestSNP[i,condition]
	gene=LD_miRtron_eQTL_bestSNP[i,gene]
	miRNA=LD_miRtron_eQTL_bestSNP[i,miRNA]
	res1[[i]]=Adjust_miRNA_QTL(miRNA,gene,snp.miRNA,cond)
	res2[[i]]=Adjust_miRNA_QTL(miRNA,gene,snp.gene,cond)
}

res1=as.data.table(rbindlist(res1))
res2=as.data.table(rbindlist(res2))

res1=merge(res1,LD_miRtron_eQTL_bestSNP,by.x=c('snp','geneID','miR','cond'),by.y=c('snps.miRNA','gene','miRNA','condition'))
res2=merge(res2,LD_miRtron_eQTL_bestSNP,by.x=c('snp','geneID','miR','cond'),by.y=c('snps.gene','gene','miRNA','condition'))
res1[,bestModel:=c('Indep','GeneCausal','miRNACausal')[which.max(c(L_Independant,L_gene_causal,L_miR_causal))],by= .(snp, geneID, miR, cond)]
res2[,bestModel:=c('Indep','GeneCausal','miRNACausal')[which.max(c(L_Independant,L_gene_causal,L_miR_causal))],by= .(snp, geneID, miR, cond)]

fwrite(res2,file=sprintf('%s/Maxime/miRNA_V2/data/15.eQTL_comparisons/miR_QTL_host_gene_eQTL_colocalization.tsv',EVO_IMMUNO_POP),sep='\t')


#################### compute distance between gene TSS and miRNA TSS ##########################

Gene_annot=fread(sprintf('%s/Maxime/Evo_Immuno_pop_data/01_GeneFPKM_cufflinks/GeneAnnot_expressed.txt',EVO_IMMUNO_POP))
colnames(Gene_annot)[1:5]=c('gene_id','seqnames','strand','start','end')
Gene_annot[,strand:=ifelse(strand>0,'+','-')]
Gene_annot[,TSS_annot:=ifelse(strand=='+',start,end)]
colnames(Gene_annot)=make.names(colnames(Gene_annot))

miRTSS=fread('/Volumes/evo_immuno_pop/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_TSS_annot_DeRie2017.txt')


x=load(sprintf("%s/Annotation/Fantom5/Promoters_all/Promoter_Mono.Rdata",HOME))

PromAnnot=fread(sprintf("%s/Annotation/Fantom5/Promoters_all/hg19.cage_peak_phase1and2combined_ann.txt",HOME),skip='chr10:100013403')
colnames(PromAnnot)=c('id','short_description','description',' association_with_transcript','entrezgene_id','hgnc_id','uniprot_id')
PromAnnot=PromAnnot[short_description%in%paste('p1@',Gene_annot$Associated.Gene.Name,sep=''),]	
PromAnnot=merge(PromAnnot,Promoters_list,by='id')
PromAnnot[,Associated.Gene.Name:=gsub('p1@','',short_description)]
PromAnnot=merge(PromAnnot,Gene_annot[!is.na(seqnames),mget(c('Associated.Gene.Name','Gene.Biotype','TSS_annot'))],by='Associated.Gene.Name')

miRNA_gene_overlap=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_intersection_with_genes.tsv',EVO_IMMUNO_POP))
miRNA_gene_overlap=merge(miRNA_gene_overlap, miRTSS[,mget(c('TSS','hsa_ID','Nb_TSS'))],by.x='miRNA',by.y='hsa_ID')
miRNA_gene_overlap=merge(miRNA_gene_overlap,PromAnnot[,mget(c('Associated.Gene.Name','TSS_annot','seqnames','start','end'))],by.x='geneNames',by.y='Associated.Gene.Name',suffix=c('miRNA','.TSSgene'))


miRNA_gene_overlap[,DistTSS_annot:=abs(TSS-TSS_annot)]
miRNA_gene_overlap[,DistTSS_CAGE:=ifelse(TSS<start.TSSgene,start.TSSgene-TSS,ifelse(TSS>end.TSSgene,TSS-end.TSSgene,0))]

mean(miRNA_gene_overlap$DistTSS_annot>1e3,na.rm=T)
# [1] 0.2142857
mean(miRNA_gene_overlap$DistTSS_CAGE>1e3,na.rm=T)
# [1] 0.07453416
# The TSS are largely the same between host gene and miRNAs

hist(log10(abs(miRNA_gene_overlap[,.(DistTSS=TSS_annot-TSS)]$DistTSS)),br=100)


################################################################
##	marginal association of host eQTL with miRNAs	  		  ##
################################################################

miRNA_gene_overlap=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/miRNA_intersection_with_genes.tsv',EVO_IMMUNO_POP))

eQTL=fread(sprintf('%s/Maxime/evo_immuno_pop_QTLs/eQTL/ALL/eQTL_ALL_bestSNP_FDR5_5cond.txt',EVO_IMMUNO_POP))


miRNA_gene_overlap_eQTL = merge(miRNA_gene_overlap,eQTL,by.x='ensemblGenes',by.y='gene')
unique(miRNA_gene_overlap_eQTL,by=c('miRNA','ensemblGenes'))
# 81 miRNA in genes with an eQTL

SNP_genos=as.data.table(getGenos(unique(c(miRNA_gene_overlap_eQTL$snps))))
SNP_genos=melt(SNP_genos, id.vars=c("snp.name","chromosome","position","locus","posID"),variable.name='individual',value.name='genotype')
SNP_genos[,genotype:=2-genotype]


test_miR_SNP=function(i){
	miR=miRNA_gene_overlap_eQTL[i,miRNA]
	snps=miRNA_gene_overlap_eQTL[i,snps]
	cond=miRNA_gene_overlap_eQTL[i,condition]
	data_to_use = merge(miRNA_expression[miRNA==miR & condition==cond,],SNP_genos[snp.name==snps,],by='individual')
	mod=lm(expression~genotype+population,data=data_to_use)
	Pval=summary(mod)$coeff[2,4]
	beta=summary(mod)$coeff[2,1]
	tstat=summary(mod)$coeff[2,3]
	data.table(miR,cond,snps,Pval, beta, tstat)
}
miRNA_gene_overlap_eQTL=miRNA_gene_overlap_eQTL[miRNA%chin%miRNA_expression$miRNA,]
test_miR_Hostgene_eQTL=lapply(1:nrow(miRNA_gene_overlap_eQTL),test_miR_SNP)
test_miR_Hostgene_eQTL=rbindlist(test_miR_Hostgene_eQTL)
Best_miR_Hostgene_eQTL=unique(test_miR_Hostgene_eQTL[order(miR,Pval),],by=c('miR'))
Best_miR_Hostgene_eQTL[fdr<.01,]

