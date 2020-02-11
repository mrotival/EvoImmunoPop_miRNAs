library(data.table)
library(stringr)

mirAnnotFile=sprintf('%s/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_raw_counts.tsv',EVO_IMMUNO_POP)

isoMir=fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/05.isomirs_count_correction/miRNA_isomiR_raw_counts_aggregated_nosubs_cast.tsv", sep=""))
# 2722 isomiRs inlcuding 2068 merged isomiRs (frequent) and 654 other categories 
# 19 frequent isomiRs have no other category and are excluded. 

samples=colnames(isoMir)[grep('AFB|EUB',colnames(isoMir))]
isoMir_mat=as.matrix(as.data.frame(isoMir[,mget(samples)]))
rownames(isoMir_mat)=isoMir$isomir_ID

all_test_var_cat = fread(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/covariates/all_categorical_variables.tsv", sep=""))
all_test_var_cat = all_test_var_cat[match(colnames(isoMir_mat),V1),]


condIndex=c("NS", "LPS", "PAM3CSK4", "R848", "IAV")
for (i in 1:5){
    isoMir[,paste('mean_isoMiR_',condIndex[i],sep='')]=apply(isoMir_mat[,all_test_var_cat$condition==i],1,mean)
}

isoMir[,mirID:=paste(hsa_ID,MIMAT,MI_ID)]
isoMir$mean_isoMiR=apply(isoMir_mat,1,mean)
isoMir$ratio_isoMiR=as.vector(isoMir$mean_isoMiR/by(isoMir$mean_isoMiR,isoMir$mirID,sum,na.rm=T)[isoMir$mirID]*100)

for (i in 1:5){
    expr_isomiR=isoMir[,paste('mean_isoMiR_',condIndex[i],sep='')]
    isoMir[,paste('ratio_isoMiR_',condIndex[i],sep='')]=as.vector(expr_isomiR/by(expr_isomiR,isoMir$mirID,sum,na.rm=T)[isoMir$mirID]*100)
}

isoMir = isoMir[-grep('other',isoMir$isomir_ID),]
#isoMir$V1=isoMir$isomir_ID
#isoMir$isomir_ID=NULL

Data=str_split_fixed(isoMir$isomir_ID,'_',11)
shift_5p=as.numeric(Data[,2])
shift_3p=as.numeric(Data[,3])
mirID=gsub('_*$','',apply(Data[,7:11],1,paste,collapse='_'))
MIMAT=gsub('.*(MIMAT[0-9]+(_[0-9])?).*','\\1',mirID)
MI_ID=gsub('.*(MIMAT[0-9]+(_[0-9])?)_(MI[0-9]+)','\\3',mirID)
hsa_ID=gsub('(.*)_(MIMAT[0-9]+(_[0-9])?)_(MI[0-9]+)','\\1',mirID)
strand=Data[,4]
isomir_sequence=Data[,5]
isomir_seed=substr(isomir_sequence,2,7)
isomiR_subs=substr(Data[,6],3,100)

## load miRNA annot
miRNA_Annot=fread(mirAnnotFile)
miRNA_Annot = miRNA_Annot[,.(mir_1=miRNA_name, mir_2=miRNA_name2,hairpin= miRNA_hairpin)]
miRNA_Annot=miRNA_Annot[!duplicated(miRNA_Annot)]

# all substitutions, one sub by isoMIR
subs=str_match_all(isoMir$isomir_ID,'([0-9]+)([ATGC])->([ATGC])')
nsubs=sapply(subs,nrow)
subs_position=unlist(sapply(subs,function(x){as.numeric(x[,2])}))
subs_position_3p=rep(nchar(Data[,5]),sapply(subs,nrow))-subs_position+1
subs_type=unlist(sapply(subs,function(x){paste(x[,3],x[,4],sep='->')}))
strand_rep=rep(strand,sapply(subs,nrow))
shift_5P_subs=rep(shift_5p,sapply(subs,nrow))
shift_3P_subs=rep(shift_3p,sapply(subs,nrow))

subs_miR_id=rep(isoMir$isomir_ID,sapply(subs,nrow))
nta1=(subs_position_3p==1 & shift_3P_subs==1 & shift_5P_subs==0)
nta22=(subs_position_3p==2 & shift_3P_subs==2 & shift_5P_subs==0)
nta12=(subs_position_3p==1 & shift_3P_subs==2 & shift_5P_subs==0)

extractSubs=function(subsCriteria){
    subs_miR=by(subsCriteria,subs_miR_id,any)
    subs_miR=subs_miR[match(isoMir$isomir_ID,names(subs_miR))]
    subs_miR[is.na(subs_miR)]=FALSE
    as.vector(subs_miR)
    }

any_to_any_LastPos=extractSubs(subs_position_3p==1)
any_to_any_BeforeLastPos=extractSubs(subs_position_3p==2)
any_to_any_BeforeBeforeLastPos=extractSubs(subs_position_3p==3)

######## annotate miRNAs

# non template additions (strict)
is_nta_miR=by(nta1,subs_miR_id,any) | (by(nta12,subs_miR_id,any) & by(nta12,subs_miR_id,any))
is_nta_miR=is_nta_miR[match(isoMir$isomir_ID,names(is_nta_miR))]
is_nta_miR[is.na(is_nta_miR)]=FALSE

# classify isomiRs
isomiR_type = ifelse(nsubs==0 & shift_5p==0 & shift_3p==0, 'CAN',
                 ifelse(nsubs==0 & shift_5p!=0 & shift_3p==0, '5PC',
                     ifelse(nsubs==0 & shift_5p==0 & shift_3p!=0, '3PC',
                         ifelse(nsubs==0 & shift_5p!=0 & shift_3p!=0 & shift_3p==shift_5p, 'SFT',
                             ifelse(nsubs!=0 & shift_5p==0 & shift_3p==0, 'SUBS',
                                 ifelse(nsubs!=0 & shift_5p==0 & shift_3p>0 & is_nta_miR, 'NTA','OTH'))))))
miR_seed=isomir_seed[isomiR_type=='CAN']
names(miR_seed)=mirID[isomiR_type=='CAN']


# classify isomiRs (complex)
isomiR_type_cmplx = ifelse(nsubs==0 & shift_5p==0 & shift_3p==0, 'CAN',
                 ifelse(nsubs==0 & shift_5p<0 & shift_3p==0, '5P_EXT',
                 ifelse(nsubs==0 & shift_5p>0 & shift_3p==0, '5P_RED',
                     ifelse(nsubs==0 & shift_5p==0 & shift_3p>0, '3P_EXT',
                     ifelse(nsubs==0 & shift_5p==0 & shift_3p<0, '3P_RED',
                         ifelse(nsubs==0 & shift_5p!=0 & shift_3p!=0 & shift_3p==shift_5p, 'SFT',
                         ifelse(nsubs==0 & shift_5p!=0 & shift_3p!=0 & shift_3p!=shift_5p, 'CPX_SFT',
                             ifelse(nsubs!=0 & shift_5p==0 & shift_3p==0, 'SUBS',
                                 ifelse(nsubs!=0 & shift_5p==0 & shift_3p>0 & is_nta_miR, 'NTA','OTH')))))))))

# miRNA arm
miRNA_arm=ifelse(grepl('3p_MIMAT',isoMir$isomir_ID),'3p',ifelse(grepl('5p_MIMAT',isoMir$isomir_ID),'5p',NA))
miRNA_arm_usage=apply(2^isoMir_mat-1,2,function(x){by(x,miRNA_arm,sum)/sum(x)})
miRNA_arm_num=ifelse(miRNA_arm=='3p',1,0)
# non template additions in 3p (any) 
nta_3p1=by(subs_position_3p==1,subs_miR_id,any)
nta_3p2=by(subs_position_3p==1,subs_miR_id,any) & by(subs_position_3p==2,subs_miR_id,any)
nta_3p_miR=nta_3p1 + nta_3p2
nta_3p_miR=nta_3p_miR[match(isoMir$isomir_ID,names(nta_3p_miR))]
nta_3p_miR[is.na(nta_3p_miR)]=0

# non template additions in 3p (A nucleotide)
nta_3p1_a=by(subs_position_3p==1 & substr(subs_type,4,4)=='A',subs_miR_id,any)
nta_3p2_a=by(subs_position_3p==1 & substr(subs_type,4,4)=='A',subs_miR_id,any) & by(subs_position_3p==2 & substr(subs_type,4,4)=='A',subs_miR_id,any)
nta_3p_a_miR=nta_3p1_a + nta_3p2_a
nta_3p_a_miR=nta_3p_a_miR[match(isoMir$isomir_ID,names(nta_3p_a_miR))]
nta_3p_a_miR[is.na(nta_3p_a_miR)]=0

# non template additions in 3p (U nucleotide)
nta_3p1_u=by(subs_position_3p==1 & substr(subs_type,4,4)=='T',subs_miR_id,any)
nta_3p2_u=by(subs_position_3p==1 & substr(subs_type,4,4)=='T',subs_miR_id,any) & by(subs_position_3p==2 & substr(subs_type,4,4)=='T',subs_miR_id,any)
nta_3p_u_miR=nta_3p1_u + nta_3p2_u
nta_3p_u_miR=nta_3p_u_miR[match(isoMir$isomir_ID,names(nta_3p_u_miR))]
nta_3p_u_miR[is.na(nta_3p_u_miR)]=0

shift_3p_template=shift_3p-as.vector(nta_3p_miR)

isomiR_annot=data.frame(ID=isoMir$isomir_ID,mirID,hsa_ID,MIMAT,MI_ID,miRNA_arm,miRNA_arm_num,
                        shift_3p,
                        is_3p_extended=shift_3p>0,
                        is_3p_reduced=shift_3p<0,
                        is_3p_cannonical=shift_3p==0,
                        shift_5p,
                        is_5p_extended=shift_5p<0,
                        is_5p_reduced=shift_5p>0,
                        is_5p_canonical=shift_5p==0,
                        isomiR_type,isomiR_type_cmplx,
                        is_nta_miR=as.vector(is_nta_miR),
                        nta_3p_miR=as.vector(nta_3p_miR),
                        nta_3p_a_miR=as.vector(nta_3p_a_miR),
                        nta_3p_u_miR=as.vector(nta_3p_u_miR),
                        shift_3p_template=shift_3p_template,
                        shift_3p_template_ext=pmax(0,shift_3p_template),
                        shift_3p_template_red=pmin(0,shift_3p_template),
                        is_3p_template_extended=shift_3p_template>0,
                        is_3p_template_reduced=shift_3p_template<0,
                        is_3p_template_canonical=shift_3p_template==0,
                        nsubs,isomiR_subs,
                        any_to_any_BeforeBeforeLastPos,any_to_any_BeforeLastPos,any_to_any_LastPos,
                        isomir_sequence, isomir_seed,miR_seed=miR_seed[mirID],
                        has_canonical_seed=(isomir_seed==miR_seed[mirID]),
                        is_cannonical=(isomiR_type=='CAN'))

isomiR_annot$mean_isoMiR=apply(isoMir_mat,1,mean)

Start_site_shift=isomiR_annot$shift_5p
Start_site_shift[isomiR_annot$shift_5p < -2]='< -2'
Start_site_shift[isomiR_annot$shift_5p>2]='> 2'
isomiR_annot$Start_site_shift=factor(Start_site_shift,levels=c('< -2','-2','-1','0','1','2','> 2'))

End_site_shift=isomiR_annot$shift_3p
End_site_shift[isomiR_annot$shift_3p < -2]='< -2'
End_site_shift[isomiR_annot$shift_3p>2]='> 2'
isomiR_annot$End_site_shift=factor(End_site_shift,levels=c('< -2','-2','-1','0','1','2','> 2'))

End_site_shift_template=isomiR_annot$shift_3p_template
End_site_shift_template[isomiR_annot$shift_3p_template < -2]='< -2'
End_site_shift_template[isomiR_annot$shift_3p_template>2]='> 2'
isomiR_annot$End_site_shift_template=factor(End_site_shift_template,levels=c('< -2','-2','-1','0','1','2','> 2'))

fwrite(isomiR_annot,file=sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_nosubs_FULL.tsv',EVO_IMMUNO_POP),sep='\t')
#fwrite(isomiR_annot,file=sprintf('%s/Martin/miRNA/17.annotate_miRNAs.R/results/isomiR_annotation.tsv',EVO_IMMUNO_POP),sep='\t')
