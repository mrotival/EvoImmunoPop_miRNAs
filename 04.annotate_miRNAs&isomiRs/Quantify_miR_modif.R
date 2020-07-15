#########################################################
###  Quantify miRNA modifications                     ###
#########################################################
isomiR_annot = as.data.frame(fread(sprintf('%s/Martin/miRNA/17.annotate_miRNAs.R/results/isomiR_annotation_FULL_V2.0.tsv',EVO_IMMUNO_POP)))
mirID=isomiR_annot$mirID

miRNA_arm_miR=ifelse(grepl('3p_MIMAT',sort(unique(mirID))),'3p',ifelse(grepl('5p_MIMAT',sort(unique(mirID))),'5p',NA))
#isomiR_annot$mean_isoMiR=apply(2^isoMir_mat-1,1,mean)

getModif_frequency=function(miR_modif,cutoff=c(.01,.05,.1,.5,.80,.95,.99),ratio_Var='ratio_isoMiR'){
    isomiR_annot=as.data.table(isomiR_annot)
    isomiR_annot$miRNA_modifs=miR_modif
    Expr_modif=dcast(isomiR_annot[,.(ratio=sum(get(eval(ratio_Var)))/100),by=.(mirID,miRNA_modifs)], mirID~miRNA_modifs,fill=0)
    Expr_total=isomiR_annot[,.(ratio_total=sum(get(eval(ratio_Var)))),by=.(mirID)]
    Expr_modif = merge(Expr_total, Expr_modif,by='mirID')
    Expr_modif=melt(Expr_modif,id.vars=c('mirID','ratio_total'),measure.vars=levels(isomiR_annot$miRNA_modifs),variable.name='miRNA_modifs',value.name='ratio')

    modification_Frequency=sapply(cutoff,function(x){table(Expr_modif[ratio>x,'miRNA_modifs'])/length(unique(Expr_modif$mirID))*100})
    colnames(modification_Frequency)=paste(cutoff*100,'%',sep='')
    modification_Frequency
    }

getModifVects_frequency=function(miR_modif_mat,cumulative=TRUE){
    modification_Frequency=list()
    for (i in 1:ncol(miR_modif_mat)){
        isomiR_annot$miRNA_modifs=miR_modif_mat[,i]
        isomiR_annot=as.data.table(isomiR_annot)
        Expr_modif=dcast(isomiR_annot[,.(ratio=sum(ratio_isoMiR)/100),by=.(mirID,miRNA_modifs)], mirID~miRNA_modifs,fill=0)
        Expr_total=isomiR_annot[,.(ratio_total=sum(ratio_isoMiR)),by=.(mirID)]
        Expr_modif = merge(Expr_total, Expr_modif,by='mirID')
        Expr_modif=melt(Expr_modif,id.vars=c('mirID','ratio_total'),measure.vars=levels(isomiR_annot$miRNA_modifs),variable.name='miRNA_modifs',value.name='ratio')
        modification_Frequency_i=sapply(c(.01,.05,.1,.5),function(x){table(Expr_modif[ratio>x,'miRNA_modifs'])/length(unique(Expr_modif$mirID))*100})
        colnames(modification_Frequency_i)=c("1%","5%","10%","50%")
        modification_Frequency[[i]]=as.data.frame(modification_Frequency_i)['TRUE',]
        }
    modification_Frequency=rbindlist(modification_Frequency)
    modification_Frequency=as.matrix(modification_Frequency)
    rownames(modification_Frequency)=colnames(miR_modif_mat)
    if(!cumulative){
        modification_Frequency[,1]=modification_Frequency[,1]-modification_Frequency[,2]
        modification_Frequency[,2]=modification_Frequency[,2]-modification_Frequency[,3]
#        modification_Frequency[,3]=modification_Frequency[,3]-modification_Frequency[,4]
    }
    modification_Frequency
    
    }

getModifVects_frequency_custom=function(miR_modif_mat,cutoff=c(.01,.05,.1,.5,.95,.99)){
    modification_Frequency=list()
    for (i in 1:ncol(miR_modif_mat)){
        isomiR_annot$miRNA_modifs=miR_modif_mat[,i]
        isomiR_annot=as.data.table(isomiR_annot)
        Expr_modif=dcast(isomiR_annot[,.(ratio=sum(ratio_isoMiR)/100),by=.(mirID,miRNA_modifs)], mirID~miRNA_modifs,fill=0)
        Expr_total=isomiR_annot[,.(ratio_total=sum(ratio_isoMiR)),by=.(mirID)]
        Expr_modif = merge(Expr_total, Expr_modif,by='mirID')
        Expr_modif=melt(Expr_modif,id.vars=c('mirID','ratio_total'),measure.vars=levels(isomiR_annot$miRNA_modifs),variable.name='miRNA_modifs',value.name='ratio')
        modification_Frequency_i=sapply(cutoff,function(x){table(Expr_modif[ratio>x,'miRNA_modifs'])/length(unique(Expr_modif$mirID))*100})
        colnames(modification_Frequency_i)=paste(cutoff*100,'%')
        modification_Frequency[[i]]=as.data.frame(modification_Frequency_i)['TRUE',]
        }
    modification_Frequency=rbindlist(modification_Frequency)
    modification_Frequency=as.matrix(modification_Frequency)
    rownames(modification_Frequency)=colnames(miR_modif_mat)
    modification_Frequency
   }

getModifVects_lowfrequency=function(miR_modif_mat){
    modification_Frequency=list()
    for (i in 1:ncol(miR_modif_mat)){
        isomiR_annot$miRNA_modifs=miR_modif_mat[,i]
        isomiR_annot=as.data.table(isomiR_annot)
        Expr_modif=dcast(isomiR_annot[,.(ratio=sum(ratio_isoMiR)/100),by=.(mirID,miRNA_modifs)], mirID~miRNA_modifs,fill=0)
        Expr_total=isomiR_annot[,.(ratio_total=sum(ratio_isoMiR)),by=.(mirID)]
        Expr_modif = merge(Expr_total, Expr_modif,by='mirID')
        Expr_modif=melt(Expr_modif,id.vars=c('mirID','ratio_total'),measure.vars=levels(isomiR_annot$miRNA_modifs),variable.name='miRNA_modifs',value.name='ratio')

        modification_Frequency_i=sapply(c(0.001,.01),function(x){table(Expr_modif[ratio>0 & ratio<x,'miRNA_modifs'])/length(unique(Expr_modif$mirID))*100})
        colnames(modification_Frequency_i)=c('<.1%',"<1%")
        modification_Frequency[[i]]=as.data.frame(modification_Frequency_i)['TRUE',]
        }
    modification_Frequency=rbindlist(modification_Frequency)
    modification_Frequency=as.matrix(modification_Frequency)
    rownames(modification_Frequency)=colnames(miR_modif_mat)
    modification_Frequency
    }
    
barplot_modif=function(modif_var,modif_name,modif_short_name,Palette='BuPu',ratio_Var='ratio_isoMiR',...){
    if(!is.matrix(modif_var)){
        modification_Frequency=getModif_frequency(modif_var,cutoff=c(.01,.05,.1,.5),ratio_Var=ratio_Var)
        }else{
        modification_Frequency=getModifVects_frequency(modif_var,ratio_Var=ratio_Var)
        }
    require(RColorBrewer)
    mypalette <- brewer.pal(4,Palette)
    pdf(sprintf('%s/Maxime/miRNA_V2/figures/Revisions/miR_modifs/frequency_of_isomiRs_%s.pdf',EVO_IMMUNO_POP,modif_short_name),...)
    par(mar=c(4,7,4,1))
    barplot(modification_Frequency[,1],ylim=c(0,100),las=1,col=mypalette[1],border=NA,ylab='Percentage of miRNAs where the \n modification affects > X % of the reads',xlab=modif_name)
    barplot(modification_Frequency[,2],add=T,col=mypalette[2],axes=F,names.arg='',border=NA)
    barplot(modification_Frequency[,3],add=T,col=mypalette[3],axes=F,names.arg='',border=NA)
    barplot(modification_Frequency[,4],add=T,col=mypalette[4],axes=F,names.arg='',border=NA)
    barplot(modification_Frequency[,1],add=T,col="#00000000",axes=F,names.arg='')

    abline(h=100,col='lightgrey',lty=3)
    ncol_Legend=ifelse(nrow(modification_Frequency)<4,2,4)
    ncol_height=ifelse(nrow(modification_Frequency)<4,-0.4,-0.2)
    legend('top',fill=mypalette,legend=c("1%","5%","10%","50%"),xpd=T,bty='n',ncol=ncol_Legend,inset=-0.2)
    dev.off()
    DF=data.frame(ID=rownames(modification_Frequency),modification_Frequency)
    colnames(DF)[1]=modif_short_name
    fwrite(DF,file=sprintf('%s/Maxime/miRNA_V2/figures/Revisions/miR_modifs/frequency_of_isomiRs_%s.txt',EVO_IMMUNO_POP,modif_short_name),sep='\t')
}

barplot_modif_subs=function(modif_var,modif_name,modif_short_name,ratio_Var='ratio_isoMiR',...){
    if(!is.matrix(modif_var)){
        modification_Frequency=getModif_frequency(modif_var,ratio_Var=ratio_Var)
        }else{
        modification_Frequency=getModifVects_frequency(modif_var,ratio_Var=ratio_Var)
        }
    require(RColorBrewer)
    mypalette <- brewer.pal(4,"BuPu")
    pdf(sprintf('%s/Maxime/miRNA_V2/figures/Revisions/miR_modifs/frequency_of_isomiRs_%s.pdf',EVO_IMMUNO_POP,modif_short_name),...)
    par(mar=c(4,7,4,1))
    rownames(modification_Frequency)=rep('',nrow(modification_Frequency))
    myspace=rep(c(.2,1,.2),c(10,1,9))
    x=barplot(modification_Frequency[,1],ylim=c(0,100),col=mypalette[1],border=NA,las=1,ylab='Percentage of miRNAs where the \n modification affects > X % of the reads',xlab=modif_name,space=myspace,cex.names=0.8)
    barplot(modification_Frequency[,2],add=T,col=mypalette[2],axes=F,names.arg='',border=NA,space=myspace)
    barplot(modification_Frequency[,3],add=T,col=mypalette[3],axes=F,names.arg='',border=NA,space=myspace)
    barplot(modification_Frequency[,4],add=T,col=mypalette[4],axes=F,names.arg='',border=NA,space=myspace)
    barplot(modification_Frequency[,1],add=T,col="#00000000",axes=F,names.arg='',space=myspace)
    text(x,-8,labels=c(1:10,-(10:1)),xpd=T)
    arrows(x[1]-0.5,-14,x[10]+0.7,-14,xpd=T,length=0.08)
    arrows(x[20]+0.5,-14,x[11]-0.7,-14,xpd=T,length=0.08)
    text(x[1],-16,labels="5' start",xpd=T,adj=c(0,1))
    text(x[20],-16,labels="3' end",xpd=T,adj=c(1,1))
    abline(h=100,col='lightgrey',lty=3)
    ncol_Legend=ifelse(nrow(modification_Frequency)<4,2,4)
    ncol_height=ifelse(nrow(modification_Frequency)<4,-0.4,-0.2)
    legend('top',fill=mypalette,legend=c("1%","5%","10%","50%"),xpd=T,bty='n',ncol=ncol_Legend,inset=-0.2)
    dev.off()
    DF=data.frame(ID=rownames(modification_Frequency),modification_Frequency)
    colnames(DF)[1]=modif_short_name
    fwrite(DF,file=sprintf('%s/Maxime/miRNA_V2/figures/Revisions/miR_modifs/frequency_of_isomiRs_%s.txt',EVO_IMMUNO_POP,modif_short_name),sep='\t')
}


Subs_pos_var=colnames(isomiR_annot)[grep('subs_Pos_',colnames(isomiR_annot))]
Subs_Pos=as.matrix(as.data.frame(isomiR_annot)[,Subs_pos_var])
rownames(Subs_Pos)=isomiR_annot$ID

barplot_modif_subs(Subs_Pos, modif_name="position of substitution", modif_short_name="subs_pos", height=4, width=8.5)

Subs_count=getModifVects_frequency(Subs_Pos)*length(unique(isomiR_annot$mirID))/100
# test for enrichment of susbtitution at sepcific positions
pbinom(t(Subs_count),apply(Subs_count,2,sum),low=F,0.05)#
#    subs_Pos_5p_1 subs_Pos_5p_2 subs_Pos_5p_3 subs_Pos_5p_4 subs_Pos_5p_5 subs_Pos_5p_6 subs_Pos_5p_7 subs_Pos_5p_8 subs_Pos_5p_9 subs_Pos_5p_10 subs_Pos_3p_10 subs_Pos_3p_9 subs_Pos_3p_8 subs_Pos_3p_7
#1%   2.163288e-05  1.915712e-07     1.0000000  1.306264e-37     1.0000000     1.0000000     1.0000000     1.0000000     1.0000000      1.0000000      1.0000000     1.0000000     1.0000000     1.0000000
#5%   1.515465e-01  9.999997e-01     1.0000000  3.763859e-33     1.0000000     1.0000000     1.0000000     1.0000000     1.0000000      1.0000000      1.0000000     1.0000000     1.0000000     1.0000000
#10%  2.002609e-01  9.999936e-01     1.0000000  2.484143e-17     1.0000000     1.0000000     1.0000000     1.0000000     0.9999998      1.0000000      0.9999998     1.0000000     1.0000000     1.0000000
#50%  8.972354e-01  7.692049e-01     0.9951777  9.951777e-01     0.9951777     0.9951777     0.9951777     0.9951777     0.7692049      0.9687819      0.9687819     0.9951777     0.9951777     0.9951777
#    subs_Pos_3p_6 subs_Pos_3p_5 subs_Pos_3p_4 subs_Pos_3p_3 subs_Pos_3p_2 subs_Pos_3p_1
#1%      1.0000000     1.0000000     1.0000000     0.9994535  3.570277e-38 6.064733e-295
#5%      1.0000000     1.0000000     1.0000000     0.9999843  1.909833e-07 7.326662e-292
#10%     1.0000000     1.0000000     1.0000000     0.9999064  8.165217e-04 1.379736e-257
#50%     0.9951777     0.9687819     0.9951777     0.8972354  3.562364e-02  1.823060e-87


isomiR_annot$Start_site_shift=factor(isomiR_annot$Start_site_shift,levels=c('< -2','-2','-1','0','1','2','> 2'))
isomiR_annot$End_site_shift=factor(isomiR_annot$End_site_shift,levels=c('< -2','-2','-1','0','1','2','> 2'))
isomiR_annot$End_site_shift_template=factor(isomiR_annot$End_site_shift_template,levels=c('< -2','-2','-1','0','1','2','> 2'))
isomiR_annot$Start_site_shift_template=factor(isomiR_annot$Start_site_shift_template,levels=c('< -2','-2','-1','0','1','2','> 2'))

DF=list()
modification_Frequency=getModif_frequency(isomiR_annot$Start_site_shift)
DF[['Start_site_shift']]=data.frame(type='Start_site_shift',ID=rownames(modification_Frequency),modification_Frequency,ratio_Var='ratio_isoMiR')
modification_Frequency=getModif_frequency(isomiR_annot$End_site_shift)
DF[['End_site_shift']]=data.frame(type='End_site_shift',ID=rownames(modification_Frequency),modification_Frequency,ratio_Var='ratio_isoMiR')
modification_Frequency=getModif_frequency(isomiR_annot$End_site_shift_template)
DF[['End_site_shift_template']]=data.frame(type='End_site_shift_template',ID=rownames(modification_Frequency),modification_Frequency,ratio_Var='ratio_isoMiR')
modification_Frequency=getModif_frequency(isomiR_annot$Start_site_shift_template)
DF[['Start_site_shift_template']]=data.frame(type='Start_site_shift_template',ID=rownames(modification_Frequency),modification_Frequency,ratio_Var='ratio_isoMiR')

modification_Frequency=getModif_frequency(isomiR_annot$Start_site_shift,ratio_Var='ratio_isoMiR_NS')
DF[['Start_site_shift_NS']]=data.frame(type='Start_site_shift',ID=rownames(modification_Frequency),modification_Frequency,ratio_Var='ratio_isoMiR_NS')
modification_Frequency=getModif_frequency(isomiR_annot$End_site_shift,ratio_Var='ratio_isoMiR_NS')
DF[['End_site_shift_NS']]=data.frame(type='End_site_shift',ID=rownames(modification_Frequency),modification_Frequency,ratio_Var='ratio_isoMiR_NS')
modification_Frequency=getModif_frequency(isomiR_annot$End_site_shift_template,ratio_Var='ratio_isoMiR_NS')
DF[['End_site_shift_template_NS']]=data.frame(type='End_site_shift_template',ID=rownames(modification_Frequency),modification_Frequency,ratio_Var='ratio_isoMiR_NS')
modification_Frequency=getModif_frequency(isomiR_annot$Start_site_shift_template,ratio_Var='ratio_isoMiR_NS')
DF[['Start_site_shift_template_NS']]=data.frame(type='Start_site_shift_template',ID=rownames(modification_Frequency),modification_Frequency,ratio_Var='ratio_isoMiR_NS')

modification_Frequency=getModif_frequency(isomiR_annot$Start_site_shift,ratio_Var='ratio_isoMiR_IAV')
DF[['Start_site_shift_IAV']]=data.frame(type='Start_site_shift',ID=rownames(modification_Frequency),modification_Frequency,ratio_Var='ratio_isoMiR_IAV')
modification_Frequency=getModif_frequency(isomiR_annot$End_site_shift,ratio_Var='ratio_isoMiR_IAV')
DF[['End_site_shift_IAV']]=data.frame(type='End_site_shift',ID=rownames(modification_Frequency),modification_Frequency,ratio_Var='ratio_isoMiR_IAV')
modification_Frequency=getModif_frequency(isomiR_annot$End_site_shift_template,ratio_Var='ratio_isoMiR_IAV')
DF[['End_site_shift_template_IAV']]=data.frame(type='End_site_shift_template',ID=rownames(modification_Frequency),modification_Frequency,ratio_Var='ratio_isoMiR_IAV')
modification_Frequency=getModif_frequency(isomiR_annot$Start_site_shift_template,ratio_Var='ratio_isoMiR_IAV')
DF[['Start_site_shift_template_IAV']]=data.frame(type='Start_site_shift_template',ID=rownames(modification_Frequency),modification_Frequency,ratio_Var='ratio_isoMiR_IAV')


fwrite(rbindlist(DF),file=sprintf('%s/Maxime/miRNA_V2/figures/Revisions/miR_modifs/frequency_of_isomiRs_AllSHIFTS_new.txt',EVO_IMMUNO_POP),sep='\t')

barplot_modif(isomiR_annot$Start_site_shift, modif_name="shift of 5' start site", modif_short_name="5p_shift",height=4,width=5)
barplot_modif(isomiR_annot$End_site_shift, modif_name="shift of 3' end site", modif_short_name="3p_shift",height=4,width=5)
barplot_modif(isomiR_annot$End_site_shift_template, modif_name="shift of 3' template end site ", modif_short_name="3p_shift_template",height=4,width=5)
# barplot_modif(isomiR_annot$Start_site_shift_template, modif_name="shift of 5' template start site ", modif_short_name="5p_shift_template",height=4,width=5)

barplot_modif(isomiR_annot$nta_3p_a_miR, modif_name="end site adenylation", modif_short_name="3p_adenylation",height=4,width=3)
barplot_modif(isomiR_annot$nta_3p_u_miR, modif_name="end site uridylation", modif_short_name="3p_uridylation",height=4,width=3)
# barplot_modif(isomiR_annot$nta_3p_o_miR, modif_name="non template additions - other", modif_short_name="3p_nta_other",height=4,width=3)
# barplot_modif(isomiR_annot$nta_3p_miR, modif_name="non template additions", modif_short_name="3p_nta",height=4,width=3)


barplot_modif(isomiR_annot$Start_site_shift, modif_name="shift of 5' start site (NS)", modif_short_name="5p_shift_NS",height=4,width=5,Palette='Greys',ratio_Var='ratio_isoMiR_NS')
barplot_modif(isomiR_annot$End_site_shift, modif_name="shift of 3' end site (NS)", modif_short_name="3p_shift_NS",height=4,width=5,Palette='Greys',ratio_Var='ratio_isoMiR_NS')
barplot_modif(isomiR_annot$End_site_shift_template, modif_name="shift of 3' template end site (NS)", modif_short_name="3p_shift_template_NS",height=4,width=5,Palette='Greys',ratio_Var='ratio_isoMiR_NS')
barplot_modif(isomiR_annot$nta_3p_a_miR, modif_name="end site adenylation (NS)", modif_short_name="3p_adenylation_NS",height=4,width=3,Palette='Greys',ratio_Var='ratio_isoMiR_NS')
barplot_modif(isomiR_annot$nta_3p_u_miR, modif_name="end site uridylation (NS)", modif_short_name="3p_uridylation_NS",height=4,width=3,Palette='Greys',ratio_Var='ratio_isoMiR_NS')

barplot_modif(isomiR_annot$Start_site_shift, modif_name="shift of 5' start site (IAV)", modif_short_name="5p_shift_IAV",height=4,width=5,Palette='OrRd',ratio_Var='ratio_isoMiR_IAV')
barplot_modif(isomiR_annot$End_site_shift, modif_name="shift of 3' end site (IAV)", modif_short_name="3p_shift_IAV",height=4,width=5,Palette='OrRd',ratio_Var='ratio_isoMiR_IAV')
barplot_modif(isomiR_annot$End_site_shift_template, modif_name="shift of 3' template end site (IAV)", modif_short_name="3p_shift_template_IAV",height=4,width=5,Palette='OrRd',ratio_Var='ratio_isoMiR_IAV')
barplot_modif(isomiR_annot$nta_3p_a_miR, modif_name="end site adenylation (IAV)", modif_short_name="3p_adenylation_IAV",height=4,width=3,Palette='OrRd',ratio_Var='ratio_isoMiR_IAV')
barplot_modif(isomiR_annot$nta_3p_u_miR, modif_name="end site uridylation (IAV)", modif_short_name="3p_uridylation_IAV",height=4,width=3,Palette='OrRd',ratio_Var='ratio_isoMiR_IAV')


###############################################################################################
################################# Quantify substitution types ################################# 
###############################################################################################

library(stringr)
subs=str_match_all(isoMir$V1,'([0-9]+)([ATGC])->([ATGC])')
subs_position_5p=unlist(sapply(subs,function(x){as.numeric(x[,2])}))
subs_position_3p=rep(nchar(isomiR_annot$isomir_sequence),sapply(subs,nrow))-subs_position_5p+1
subs_type=unlist(sapply(subs,function(x){paste(x[,3],x[,4],sep='->')}))

subs_miR_id=rep(isoMir$V1,sapply(subs,nrow))

extractSubs=function(subsCriteria){
   subs_miR=by(subsCriteria,subs_miR_id,any)
   subs_miR=subs_miR[match(isoMir$V1,names(subs_miR))]
   subs_miR[is.na(subs_miR)]=FALSE
   subs_miR
   }

subs_type_all=sort(unique(subs_type))
subs_type_all=subs_type_all[order(substr(subs_type_all,4,4),substr(subs_type_all,1,1))]
hasMutation1=sapply(subs_type_all,function(mysub){extractSubs(subs_type==mysub & subs_position_5p==1)})
hasMutation2=sapply(subs_type_all,function(mysub){extractSubs(subs_type==mysub & subs_position_5p==2)})
hasMutation4=sapply(subs_type_all,function(mysub){extractSubs(subs_type==mysub & subs_position_5p==4)})
hasMutation5=sapply(subs_type_all,function(mysub){extractSubs(subs_type==mysub & subs_position_5p==5)})

hasMutationN2=sapply(subs_type_all,function(mysub){extractSubs(subs_type==mysub & subs_position_3p==3)})
hasMutationN1=sapply(subs_type_all,function(mysub){extractSubs(subs_type==mysub & subs_position_3p==2)})
hasMutationN=sapply(subs_type_all,function(mysub){extractSubs(subs_type==mysub & subs_position_3p==1)})

colors=c(brewer.pal(3,'Reds'),brewer.pal(3,'Blues'),brewer.pal(3,'Greens'),brewer.pal(3,'Purples'))

Freq_Mut1=getModifVects_frequency(hasMutation1,cumul=TRUE)
Count_Mut1=Freq_Mut1*length(unique(isomiR_annot$mirID))/100

Freq_Mut2=getModifVects_frequency(hasMutation2,cumul=TRUE)
Count_Mut2=Freq_Mut2*length(unique(isomiR_annot$mirID))/100

Freq_Mut4=getModifVects_frequency(hasMutation4,cumul=TRUE)
Count_Mut4=Freq_Mut4*length(unique(isomiR_annot$mirID))/100
pbinom(t(Count_Mut1),apply(Count_Mut1,2,sum),low=F,1/12)#
#          C->A       G->A       T->A       A->C      G->C       T->C       A->G       C->G         T->G       A->T       C->T         G->T
#1%  0.98854583 0.95494612 0.99853493 0.99853493 0.7589941 0.99853493 0.98854583 0.99853493 2.416087e-08 0.99853493 0.99853493 9.869278e-31
#5%  0.81365533 0.58658155 0.95638742 0.95638742 0.3526267 0.95638742 0.95638742 0.95638742 9.563874e-01 0.95638742 0.95638742 8.688880e-29
#10% 0.64980781 0.64980781 0.89588881 0.89588881 0.3701703 0.89588881 0.89588881 0.89588881 8.958888e-01 0.89588881 0.89588881 3.057609e-22
#50% 0.08333333 0.08333333 0.08333333 0.08333333 0.0000000 0.08333333 0.08333333 0.08333333 8.333333e-02 0.08333333 0.08333333 8.333333e-02
 
pbinom(t(Count_Mut2),apply(Count_Mut2,2,sum),low=F,1/12)
#1%  0.9989656 0.965196780 0.9989656 0.99153643 0.96519678 0.9989656 0.99153643 0.9989656 0.06299712 0.9989656 0.9989656 1.277865e-54
#5%  0.4067078 0.009545262 0.4067078 0.08309388 0.08309388 0.4067078 0.08309388 0.4067078 0.08309388 0.4067078 0.4067078 4.067078e-01
#10% 0.3527722 0.005087770 0.3527722 0.05857767 0.05857767 0.3527722 0.05857767 0.3527722 0.35277215 0.3527722 0.3527722 3.527722e-01
#50% 0.2297454 0.019675926 0.2297454 0.01967593 0.01967593 0.2297454 0.22974537 0.2297454 0.22974537 0.2297454 0.2297454 2.297454e-01
pbinom(t(Count_Mut4),apply(Count_Mut4,2,sum),low=F,1/12)
#         C->A      G->A      T->A      A->C      G->C      T->C      A->G      C->G         T->G      A->T      C->T          G->T
#1%  0.9999996 0.9999927 0.9999996 0.9999996 0.8357579 0.9999996 0.9999996 0.9999927 0.0003337379 0.9999996 0.9999927 5.247767e-103
#5%  0.9999303 0.9999303 0.9999303 0.9999303 0.9999303 0.9999303 0.9999303 0.9992333 0.9999302959 0.9999303 0.9999303 1.950127e-119
#10% 0.9965026 0.9965026 0.9965026 0.9965026 0.9965026 0.9965026 0.9965026 0.9758362 0.9965026142 0.9965026 0.9965026  7.132126e-71
#50% 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000000 0.0000000 0.0000000  0.000000e+00

pbinom(t(Freq_Mut1*675/100)/
create_pie_plots=function(cumul){
    Freq_Mut1=getModifVects_frequency(hasMutation1,cumul=cumul)
    lowFreq_Mut1=getModifVects_lowfrequency(hasMutation1)
    Freq_Mut2=getModifVects_frequency(hasMutation2,cumul=cumul)
    lowFreq_Mut2=getModifVects_lowfrequency(hasMutation2)

    Freq_Mut4=getModifVects_frequency(hasMutation4,cumul=cumul)
    lowFreq_Mut4=getModifVects_lowfrequency(hasMutation4)

    Freq_MutN2=getModifVects_frequency(hasMutationN2,cumul=cumul)
    lowFreq_MutN2=getModifVects_lowfrequency(hasMutationN2)

    Freq_MutN1=getModifVects_frequency(hasMutationN1,cumul=cumul)
    lowFreq_MutN1=getModifVects_lowfrequency(hasMutationN1)
    Freq_MutN=getModifVects_frequency(hasMutationN,cumul=cumul)
    lowFreq_MutN=getModifVects_lowfrequency(hasMutationN)

    library(mapplots)
    par(mar=c(0.1,0.1,0.1,0.1))
    r=0.3
    d1=1
    d2=3
    plot(c(0,1+6*d1+d2),c(0,5),col='#00000000',axes=F,xlab='',ylab='',mar=c(0,0,0,0))
    add.pie(lowFreq_Mut1[,2],d1,4,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_Mut1[,1],d1,3,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_Mut1[,2],d1,2,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_Mut1[,3],d1,1,col=colors,radius=r,labels='',lty = 0)

    add.pie(lowFreq_Mut2[,2],2*d1,4,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_Mut2[,1],2*d1,3,col=colors,radius=r,labels='',lty = 0)

    add.pie(lowFreq_Mut4[,2],3*d1,4,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_Mut4[,1],3*d1,3,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_Mut4[,2],3*d1,2,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_Mut4[,3],3*d1,1,col=colors,radius=r,labels='',lty = 0)

    add.pie(lowFreq_MutN2[,2],d2+4*d1,4,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_MutN2[,1],d2+4*d1,3,col=colors,radius=r,labels='',lty = 0)

    add.pie(lowFreq_MutN1[,2],d2+5*d1,4,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_MutN1[,1],d2+5*d1,3,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_MutN1[,2],d2+5*d1,2,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_MutN1[,3],d2+5*d1,1,col=colors,radius=r,labels='',lty = 0)
    add.pie(lowFreq_MutN[,2],d2+6*d1,4,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_MutN[,1],d2+6*d1,3,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_MutN[,2],d2+6*d1,2,col=colors,radius=r,labels='',lty = 0)
    add.pie(Freq_MutN[,3],d2+6*d1,1,col=colors,radius=r,labels='',lty = 0)
    legend(3*d1+d2*0.45,4.4,ncol=1,fill=colors[c(4:12,1:3)],legend=gsub('T','U',subs_type_all[c(4:12,1:3)]),bty='n',cex=0.8)
    rect(0.5,3.5,3*d1+.5,4.5,border=NA,col='#DDDDDD88')
    rect(0.5+d2+3*d1,3.5,d2+6*d1+.5,4.5,border=NA,col='#DDDDDD88')
}

pdf(sprintf('%s/03_Analysis/miRNA_project/figures/pie_plots_substitution_events_cumulative.pdf',HOME), height=4, width=8.5)
create_pie_plots(cumul=TRUE)
text(0,4:1,labels=c('<1%','>1%','>5%','>10%'),cex=0.8)
dev.off()

pdf(sprintf('%s/03_Analysis/miRNA_project/figures/pie_plots_substitution_events_noncumulative.pdf',HOME), height=4, width=8.5)
create_pie_plots(cumul=FALSE)
text(0,4:1,labels=c('<1%','1-5%','5-10%','>10%'),cex=0.8)
dev.off()

Freq_Mut1=getModifVects_frequency(hasMutation1,cumul=TRUE)

###############################################################################################
############################# END Quantify substitution types ################################# 
###############################################################################################



########################################################################################################################## 
############################# Quantify difference in substitution types upon stimulation ################################# 
##########################################################################################################################
barplot_modif(isomiR_annot$Start_site_shift, modif_name="shift of 5' start site", modif_short_name="5p_shift",height=4,width=5)
barplot_modif(isomiR_annot$End_site_shift, modif_name="shift of 3' end site", modif_short_name="3p_shift",height=4,width=5)
barplot_modif(isomiR_annot$End_site_shift_template, modif_name="shift of 3' template end site ", modif_short_name="3p_shift_template",height=4,width=5)
barplot_modif(isomiR_annot$nta_3p_a_miR, modif_name="end site adenylation", modif_short_name="3p_adenylation",height=4,width=3)
barplot_modif(isomiR_annot$nta_3p_u_miR, modif_name="end site uridylation", modif_short_name="3p_uridylation",height=4,width=3)

modifs=c('shift_3p','shift_5p','shift_3p_template','nta_3p_miR','nta_3p_a_miR','nta_3p_u_miR')

DT=list()
for (modif in modifs){
	for (cond in condIndex[-1]){
	modif_byMiR_NS=isomiR_annot[,sum(ratio_isoMiR_NS* get(eval(modif)),na.rm=T), by='mirID']$V1
	modif_byMiR_cond=isomiR_annot[,sum(get(eval(paste("ratio_isoMiR",cond,sep='_')))* get(eval(modif))), by='mirID']$V1

	DT[[paste(cond,modif)]]=data.table(condition=cond,Modif=modif,arm='any',
						Pval_nonpaired=wilcox.test(modif_byMiR_NS,modif_byMiR_cond,paired=F)$p.value,
						mean_Delta_modif=(mean(modif_byMiR_cond,na.rm=T)-mean(modif_byMiR_NS,na.rm=T))/100,
						Pval_paired=wilcox.test(modif_byMiR_NS,modif_byMiR_cond,paired=T)$p.value)
	}
	}
for (myarm in c('3p','5p')){
for (modif in modifs){
	for (cond in condIndex[-1]){
	modif_byMiR_NS=isomiR_annot[miRNA_arm==myarm,sum(ratio_isoMiR_NS* get(eval(modif)),na.rm=T), by='mirID']$V1
	modif_byMiR_cond=isomiR_annot[miRNA_arm==myarm,sum(get(eval(paste("ratio_isoMiR",cond,sep='_')))* get(eval(modif))), by='mirID']$V1

	DT[[paste(cond,modif,myarm)]]=data.table(condition=cond,Modif=modif,arm=myarm,
						Pval_nonpaired=wilcox.test(modif_byMiR_NS,modif_byMiR_cond,paired=F)$p.value,
						mean_Delta_modif=(mean(modif_byMiR_cond,na.rm=T)-mean(modif_byMiR_NS,na.rm=T))/100,
						Pval_paired=wilcox.test(modif_byMiR_NS,modif_byMiR_cond,paired=T)$p.value)
	}
	}

}

DT=rbindlist(DT)
DT[,Fdr_nonpaired:=p.adjust(Pval_nonpaired,'fdr')]
DT[,Fdr_paired:=p.adjust(Pval_paired,'fdr')]
for (myarm in c('3p','5p','any')){
	DT[arm==myarm,global:= DT[Modif=='shift_3p' & arm==myarm, mean_Delta_modif][match(condition ,DT[Modif=='shift_3p' & arm==myarm, condition])]]
	}
DT[Modif=='shift_5p',global:=mean_Delta_modif]
DT[,pct:= mean_Delta_modif/global]

fwrite(DT,file=sprintf('%s/Maxime/miRNA_V2/figures/Revisions/miR_modifs/Change_in_miRsmodifs.txt',EVO_IMMUNO_POP),sep='\t')


