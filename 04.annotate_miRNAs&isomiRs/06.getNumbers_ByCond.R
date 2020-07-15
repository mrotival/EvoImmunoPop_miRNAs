isomiR_annot=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_FULL_V2.0.tsv',EVO_IMMUNO_POP))
isomiR_annot_nosubs=fread(sprintf('%s/Maxime/miRNA_V2/data/04.annotate_miRNAs&isomiRs/isomiR_annotation_nosubs_FULL.tsv',EVO_IMMUNO_POP))



##### counting in ratio from miRNAs (hsa_ID)
Pct_can_NS=isomiR_annot[,ratio_hsa:=mean_isoMiR_NS/sum(mean_isoMiR_NS)*100,by=hsa_ID]
Pct_can_NS=Pct_can_NS[isomiR_type=='CAN',.(Pct_can=sum(ratio_hsa)),by=hsa_ID]

Pct_can_seed_NS=isomiR_annot[,ratio_hsa := mean_isoMiR_NS/sum(mean_isoMiR_NS)*100,by=hsa_ID]
Pct_can_seed_NS=Pct_can_seed_NS[has_canonical_seed==TRUE, .(Pct_can_seed=sum(ratio_hsa)),by=hsa_ID]

Pct_can_IAV=isomiR_annot[,ratio_hsa:=mean_isoMiR_IAV/sum(mean_isoMiR_IAV)*100,by=hsa_ID]
Pct_can_IAV=Pct_can_IAV[isomiR_type=='CAN',.(Pct_can=sum(ratio_hsa)),by=hsa_ID]

Pct_can_seed_IAV=isomiR_annot[,ratio_hsa := mean_isoMiR_IAV/sum(mean_isoMiR_IAV)*100,by=hsa_ID]
Pct_can_seed_IAV=Pct_can_seed_IAV[has_canonical_seed==TRUE, .(Pct_can_seed=sum(ratio_hsa)),by=hsa_ID]

DF=merge(Pct_can_seed_NS,Pct_can_NS)
mean(DF$Pct_can<50,na.rm=T) # 0.5826446
sum(DF$Pct_can<50,na.rm=T) # 282
sum(DF$Pct_can<50 & DF$Pct_can_seed<80,na.rm=T) # 75
sum(DF$Pct_can<50 & DF$Pct_can_seed<80,na.rm=T)/sum(DF$Pct_can<50,na.rm=T) # 0.2659574
DF$cond='NS'
DF_cond=DF
library(ggplot2)
pdf(sprintf('%s/Maxime/miRNA_V2/figures/Revisions/PercentageCanonicalSeed_NS.pdf',EVO_IMMUNO_POP),height=3,width=4)
p1 = ggplot(DF,aes(x=Pct_can),size = 2.5)+geom_density(bw=3,alpha=0.6,fill=colERC[1])+coord_flip()+xlab('% of canonical isomiR')+theme_classic()
p2 = ggplot(DF,aes(x=Pct_can_seed))+geom_density(bw=3,alpha=0.6,fill=colERC[1])+coord_flip()+xlab('% of isomiR with the canonical seed')+theme_classic()
library(gridExtra)
grid.arrange(p1, p2, nrow = 1)
dev.off()

DF=merge(Pct_can_seed_IAV,Pct_can_IAV)
mean(DF$Pct_can<50,na.rm=T) # 0.5958763
sum(DF$Pct_can<50,na.rm=T) # 289
sum(DF$Pct_can<50 & DF$Pct_can_seed<80,na.rm=T) # 76
sum(DF$Pct_can<50 & DF$Pct_can_seed<80,na.rm=T)/sum(DF$Pct_can<50,na.rm=T) # 0.2629758
DF$cond='IAV'
DF_cond=rbind(DF_cond,DF)


pdf(sprintf('%s/Maxime/miRNA_V2/figures/Revisions/PercentageCanonicalSeed_IAV.pdf',EVO_IMMUNO_POP),height=3,width=4)
p1 = ggplot(DF,aes(x=Pct_can),size = 2.5)+geom_density(bw=3,alpha=0.6,fill=colERC[1])+coord_flip()+xlab('% of canonical isomiR')+theme_classic()
p2 = ggplot(DF,aes(x=Pct_can_seed))+geom_density(bw=3,alpha=0.6,fill=colERC[1])+coord_flip()+xlab('% of isomiR with the canonical seed')+theme_classic()
library(gridExtra)
grid.arrange(p1, p2, nrow = 1)
dev.off()

pdf(sprintf('%s/Maxime/miRNA_V2/figures/Revisions/PercentageCanonicalSeed_IAV.pdf',EVO_IMMUNO_POP),height=3,width=4)
p1 = ggplot(DF_cond,aes(x=Pct_can,fill=cond),size = 2.5)+geom_density(bw=3,alpha=0.6)+coord_flip()+xlab('% of canonical isomiR')+theme_classic()
p2 = ggplot(DF_cond,aes(x=Pct_can_seed,fill=cond))+geom_density(bw=3,alpha=0.6)+coord_flip()+xlab('% of isomiR with the canonical seed')+theme_classic()
library(gridExtra)
grid.arrange(p1, p2, nrow = 1)
dev.off()
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


isomiR_annot_nosubs_2=merge(isoMir[,mget(colnames(isoMir)[!grepl('EUB|AFB',colnames(isoMir))])],isomiR_annot_nosubs,by.x='isomir_ID',by.y='ID',suffix=c('','.2'))


miR_expr=By(isomiR_annot$mean_isoMiR, isomiR_annot$mirID,sum)

luq(isomiR_annot[mean_isoMiR>1 & ratio_isoMiR>1 & isomiR_type!='CAN',mirID])
#[1] 424
luq(isomiR_annot[mean_isoMiR>1 & ratio_isoMiR>1,mirID])
#[1] 492
luq(isomiR_annot[mean_isoMiR>1 & ratio_isoMiR>1 & isomiR_type=='CAN',mirID])
# [1] 390
> 424/492
# [1] 0.8617886
> 390/492
# [1] 0.7926829
NumberOfIsomiR_permiRNA_1pct=table(pmin(15,table(isomiR_annot[mean_isoMiR>1 & ratio_isoMiR>1 & ratio_isoMiR<100,mirID])))
NumberOfIsomiR_permiRNA_5pct=table(table(isomiR_annot[mean_isoMiR>1 & ratio_isoMiR>5 & ratio_isoMiR<100,mirID]))
NumberOfIsomiR_permiRNA=cbind(	NumberOfIsomiR_permiRNA_1pct[as.character(1:15)],
								NumberOfIsomiR_permiRNA_5pct[as.character(1:15)])
NumberOfIsomiR_permiRNA_before=NumberOfIsomiR_permiRNA
pdf(sprintf('%s/03_Analysis/miRNA_project/figures/Nb_of_isomiRs_per_miRNA_before.pdf',HOME),width=5,height=4)
x=barplot(rev(cumsum(rev(NumberOfIsomiR_permiRNA_1pct))),col=grey(0.8),las=1,xlab='numbers of isomiRs',ylab='Number of miRNA with >K isomiRs', axisnames=FALSE)
barplot(rev(cumsum(rev(NumberOfIsomiR_permiRNA_5pct))),add=T,col="#E31A1CAA",axes = FALSE, axisnames=FALSE)
legend('topright',fill=c(grey(0.8),"#E31A1CAA"),legend=c('at frequency >1%','at frequency >5% '),bty='n')
mtext(c(1:14,'15+'),1,at=x, line=.5)
dev.off()


#### same figure for isomiR_annot_no_subs
miR_expr=By(isomiR_annot_nosubs_2$mean_isoMiR, isomiR_annot_nosubs_2$mirID,sum)
luq(isomiR_annot_nosubs_2[mean_isoMiR>1 & ratio_isoMiR>1 & isomiR_type!='CAN',mirID])
#[1] 421
luq(isomiR_annot_nosubs_2[mean_isoMiR>1 & ratio_isoMiR>1,mirID])
#[1] 492
luq(isomiR_annot[mean_isoMiR>1 & ratio_isoMiR>1 & isomiR_type=='CAN',mirID])
# [1] 390
NumberOfIsomiR_permiRNA_1pct=table(pmin(15,table(isomiR_annot_nosubs_2[mean_isoMiR>1 & ratio_isoMiR>1 & ratio_isoMiR<100,mirID])))
NumberOfIsomiR_permiRNA_5pct=table(table(isomiR_annot_nosubs_2[mean_isoMiR>1 & ratio_isoMiR>5 & ratio_isoMiR<100,mirID]))
NumberOfIsomiR_permiRNA=cbind(	NumberOfIsomiR_permiRNA_1pct[as.character(1:15)],
								NumberOfIsomiR_permiRNA_5pct[as.character(1:15)])
NumberOfIsomiR_permiRNA_after=NumberOfIsomiR_permiRNA
pdf(sprintf('%s/03_Analysis/miRNA_project/figures/Nb_of_isomiRs_per_miRNA_nosubs.pdf',HOME),width=5,height=4)
x=barplot(rev(cumsum(rev(NumberOfIsomiR_permiRNA_1pct))),col=grey(0.8),las=1,xlab='numbers of isomiRs',ylab='Number of miRNA with >K isomiRs', axisnames=FALSE)
barplot(rev(cumsum(rev(NumberOfIsomiR_permiRNA_5pct))),add=T,col="#E31A1CAA",axes = FALSE, axisnames=FALSE)
legend('topright',fill=c(grey(0.8),"#E31A1CAA"),legend=c('at frequency >1%','at frequency >5% '),bty='n')
mtext(c(1:14,'15+'),1,at=x, line=.5)
dev.off()

DF=cbind(NumberOfIsomiR_permiRNA_before,NumberOfIsomiR_permiRNA_after)
DF=as.data.table(DF)
rownames(DF)[15]='15+'
colnames(DF)=c('Nb IsomiR 1pct','Nb IsomiR 5pct','Nb IsomiR 1pt excluding subs','Nb IsomiR 5pt excluding subs')

fwrite(DF,file=sprintf('%s/03_Analysis/miRNA_project/figures/Nb_of_isomiRs_per_miRNA_w_or_wo_subs.txt',HOME),sep='\t')

##### counting in ratio from miRNAs derived from each locus
Pct_can=isomiR_annot[isomiR_type=='CAN',.(Pct_can=sum(ratio_isoMiR)),by=mirID]
Pct_can_seed=isomiR_annot[has_canonical_seed==TRUE,.(Pct_can_seed=sum(ratio_isoMiR)),by=mirID]

DF=merge(Pct_can_seed,Pct_can)
mean(DF$Pct_can<50) # 0.5769944
sum(DF$Pct_can<50) # 311
sum(DF$Pct_can<50 & DF$Pct_can_seed<80) # 79
sum(DF$Pct_can<50 & DF$Pct_can_seed<80)/sum(DF$Pct_can<50) # 0.2540193

##### counting in ratio from miRNAs (hsa_ID)
Pct_can=isomiR_annot[,ratio_hsa:=mean_isoMiR/sum(mean_isoMiR)*100,by=hsa_ID]
Pct_can=Pct_can[isomiR_type=='CAN',.(Pct_can=sum(ratio_hsa)),by=hsa_ID]

Pct_can_seed=isomiR_annot[,ratio_hsa := mean_isoMiR/sum(mean_isoMiR)*100,by=hsa_ID]
Pct_can_seed=Pct_can_seed[has_canonical_seed==TRUE, .(Pct_can_seed=sum(ratio_hsa)),by=hsa_ID]
DF=merge(Pct_can_seed,Pct_can)
mean(DF$Pct_can<50) # 0.5896907
sum(DF$Pct_can<50) # 286
sum(DF$Pct_can<50 & DF$Pct_can_seed<80) # 75
sum(DF$Pct_can<50 & DF$Pct_can_seed<80)/sum(DF$Pct_can<50) # 0.2622378

library(ggplot2)
pdf(sprintf('%s/03_Analysis/miRNA_project/figures/PercentageCanonicalSeed.pdf',HOME),height=3,width=4)
p1 = ggplot(DF,aes(x=Pct_can),size = 2.5)+geom_density(bw=3,alpha=0.6,fill=colERC[1])+coord_flip()+xlab('% of canonical isomiR')+theme_classic()
p2 = ggplot(DF,aes(x=Pct_can_seed))+geom_density(bw=3,alpha=0.6,fill=colERC[1])+coord_flip()+xlab('% of isomiR with the canonical seed')+theme_classic()
library(gridExtra)
grid.arrange(p1, p2, nrow = 1)
dev.off()


hist(Pct_can$Pct_can,br=30)
dd=density(Pct_can$Pct_can)
nNeg=sum(dd$x<0)
n100=sum(dd$x>100)
dd$y[nNeg+1:nNeg]=dd$y[nNeg+1:nNeg]+rev(dd$y[1:nNeg])
dd$y[1:nNeg]=0
ny=length(dd$y)
dd$y[ny:1][n100+1:n100]=dd$y[ny:1][n100+1:n100]+rev(dd$y[ny:1][1:n100])
dd$y[ny:1][1:n100]=0



layout(1:2)
plot(dd$y,dd$x,col='#00000000')
polygon(c(0,dd$y[dd$y!=0],0),c(min(dd$x[dd$y!=0]),dd$x[dd$y!=0],max(dd$x[dd$y!=0])))
polygon(c(0,dd$y[dd$y!=0 & dd$x>50],0),c(min(dd$x[dd$y!=0 & dd$x>50]),dd$x[dd$y!=0 & dd$x>50],max(dd$x[dd$y!=0 & dd$x>50])),col='grey')


dds=density(Pct_can_seed$Pct_can)
nNeg=sum(dds$x<0)
n100=sum(dds$x>100)
dds$y[nNeg+1:nNeg]=dds$y[nNeg+1:nNeg]+rev(dds$y[1:nNeg])
dds$y[1:nNeg]=0
ny=length(dds$y)
dds$y[ny:1][n100+1:n100]=dds$y[ny:1][n100+1:n100]+rev(dds$y[ny:1][1:n100])
dds$y[ny:1][1:n100]=0

plot(dds$y,dds$x,col='#00000000')

polygon(c(0,dds$y[dds$y!=0],0),c(min(dds$x[dds$y!=0]),dds$x[dds$y!=0],max(dds$x[dds$y!=0])))
polygon(c(0,dds$y[dds$y!=0 & dds$x>50],0),c(min(dds$x[dds$y!=0 & dds$x>50]),dds$x[dds$y!=0 & dds$x>50],max(dds$x[dds$y!=0 & dds$x>50])),col='grey')


hh=hist(Pct_can_seed$Pct_can,breaks=seq(0,100,5))



p1 = ggplot(DF,aes(x=Pct_can)) +
    geom_histogram(data=subset(DF,Pct_can>20),fill = "darkgrey", alpha = 0.5,breaks=seq(0,100,5)) +
    geom_histogram(data=subset(DF,Pct_can<=20),fill = "red", alpha = 0.3,breaks=seq(0,100,5)) +theme_minimal()+coord_flip()+scale_y_continuous(trans='sqrt')


p2 = ggplot(aes(x=Pct_can_seed, fill=Pct_can>20, color=Pct_can>20), data=DF) + geom_histogram(alpha=0.6,breaks=seq(0,100,5)) +theme_minimal()+coord_flip()+scale_y_continuous(trans='sqrt')

ggplot(aes(x=Pct_can, fill=Pct_can>50, color=Pct_can>50), data=DF) + geom_histogram(alpha=0.6,breaks=seq(0,100,5)) +theme_minimal()+coord_flip()

ggplot(DF,aes(x=Pct_can_seed)) +
    geom_histogram(data=subset(DF,Pct_can>50),fill = "darkgrey", alpha = 0.5,breaks=seq(0,100,5)) +
    geom_histogram(data=subset(DF,Pct_can<=50),fill = "red", alpha = 0.3,breaks=seq(0,100,5)) +theme_minimal()+coord_flip()



ggplot(DF, aes(x=Pct_can, y=Pct_can_seed) )+stat_density_2d(h=3, geom = "polygon", colour="white")
 

p1 = ggplot(DF,aes(x=Pct_can,y=Pct_can_seed))+geom_contour(bw=3,alpha=0.6)+theme_minimal()+coord_flip()

p2 = ggplot(DF,aes(x=Pct_can_seed,fill=Pct_can>50, color=Pct_can>50))+geom_density(bw=3,alpha=0.6)+theme_minimal()+coord_flip()


p2 = ggplot(aes(x=Pct_can_seed, fill=Pct_can_seed>20, color=Pct_can>20), data=DF) + geom_histogram(alpha=0.6,breaks=seq(0,100,5)) +theme_minimal()+coord_flip()+scale_y_continuous(trans='sqrt')

