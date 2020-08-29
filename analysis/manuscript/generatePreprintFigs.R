library(ggpubr)

swabseq.dir='/data/Covid/swabseq/'
outdir=paste0(swabseq.dir, 'analysis/manuscript/')

#get table munging function
source(paste0(swabseq.dir, 'code/helper_functions.R')) 

#Figure 1 ========================================================================================================= 
rundir=paste0(swabseq.dir, 'runs/v11/')
dfs=mungeTables(paste0(rundir, 'countTable.RDS'))
levels(dfs$virus_copy)=as.character(round(as.numeric(levels(dfs$virus_copy))))

fig1e1=
dfs%>%filter(lysate=='MTS_TE-1_to_2') %>% filter(Stotal>1000) %>%  
ggplot(aes(x=virus_copy, y=S2+1))+geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)')+xlab('GCE per reaction')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()
fig1e2=
dfs%>%filter(lysate=='MTS_TE-1_to_2') %>% filter(Stotal>1000) %>%  
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike))+geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 Spike + 1)')+xlab('GCE per reaction')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()

ggarrange(fig1e1,fig1e2)
ggsave('/data/Covid/swabseq/analysis/manuscript/fig1e.png')
ggsave('/data/Covid/swabseq/analysis/manuscript/fig1e.svg')
#=========================================================================================================


#Figure 2 ================================================================================================
# note can reuse fig1e2 as figure 2c
#LoD Purified 
rundir=paste0(swabseq.dir, 'runs/v19/')
dfsR=mungeTables(paste0(rundir, 'countTable.RDS'))
dfsR$Plate_384=droplevels(dfsR$Plate_384)
dfsR$virus_copy[dfsR$virus_copy=='Negative Patient']='0'
dfsR$virus_copy=as.factor(round(as.numeric(as.character(dfsR$virus_copy)))) #*28/1000))
dfsR$virus_copy=droplevels(dfsR$virus_copy) #Plate_384=droplevels(dfsR$Plate_384)
dfsR$virus_copy=factor(dfsR$virus_copy, levels(dfsR$virus_copy)[order(as.numeric(levels(dfsR$virus_copy)))])
#dfsR$virus_copy=as.character(dfsR$virus_copy)
dfsR$virus_ident2=dfsR$virus_identity
dfsR$virus_ident2[grepl('^N|^P', dfsR$virus_identity)]='NegPatient'
dfsR$virus_ident2[grepl('^U', dfsR$virus_identity)]='PosPatient'

fig2a=
dfsR=dfsR %>% filter(virus_ident2!='PosPatient') %>% filter(Stotal>2000)

fig2a=ggplot(dfsR, aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+ #log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike+1)')+xlab('GCE per mL')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+ ggtitle('Nasal Swab, Purified RNA')+
    theme_bw()

# combine in data from EUA confirmation
#EUA confirmation of LoD
swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v20/')
outdir=paste0(swabseq.dir, 'analysis/EUA/')
dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfsC=dfL$dfs

#filters and reformat for plots 
dfsC= dfsC %>% 
  filter(Plate_ID=='Plate8') %>%filter(Stotal>2000)
dfsC$virus_copy=droplevels(dfsC$virus_copy) #Plate_384=droplevels(dfsC$Plate_384)
dfsC$virus_copy=factor(dfsC$virus_copy, levels(dfsC$virus_copy)[order(as.numeric(levels(dfsC$virus_copy)))])
#dfsC=dfsC[dfsC$virus_copy!='125',]

ggplot(dfsC, aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+ #log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2_spike+1)')+xlab('copies/mL')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
    theme_bw()+ggtitle('LoD Confirmation')+
    geom_hline(yintercept=.003, color='red')
dfsB=rbind(dfsR[,c('virus_copy', 'S2_normalized_to_S2_spike')],dfsC[,c('virus_copy', 'S2_normalized_to_S2_spike')])


fig2a2=ggplot(dfsB, aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+ #log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike+1)')+xlab('GCE per mL')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+ ggtitle('Nasal Swab, Purified RNA')+
    theme_bw()



#LoD Saliva 
rundir=paste0(swabseq.dir, 'runs/v13/')
dfs=mungeTables(paste0(rundir, 'countTable.RDS'))
dfsc = dfs %>% filter(lysate=='NS_Saliva_1_to_2') %>% filter(MasterMixVol=='20uL') %>% filter(Stotal>2000) 

#hooray for dynamic typing
dfsc$virus_copy=as.character(round(as.numeric(as.character(dfsc$virus_copy))))
dfsc$virus_copy=factor(dfsc$virus_copy, levels=sort(unique(as.numeric(dfsc$virus_copy))))

fig2b=
    ggplot(dfsc, aes(x=virus_copy, y=S2_normalized_to_S2_spike))+geom_quasirandom(alpha=.75)+
     scale_y_log10() + annotation_logticks(sides='l') + ylab('(S2+1)/(S2 spike+1)')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+
    xlab('GCE per reaction')+
    ggtitle('LoD - Saliva 1:1 in Normal Saline')

fig2c=fig1e2+ggtitle('LoD - Mid Turbinate 1:1 in TE')

# also contains taqpath vs luna 
rundir=paste0(swabseq.dir, 'runs/v10/')
dfs=mungeTables(paste0(rundir, 'countTable.RDS'))
#Plates1 and Plate2 are taqpath
dfs$RT=as.factor(dfs$Plate_ID %in% c("Plate1", "Plate2"))
dfs$PositivePatient=as.factor(dfs$lysate %in% c('R003', 'O-003', 'R-005','C-009', 'R-014', 'R006', 'O-009', 'R-004', 'O-014', 'C-008'))
levels(dfs$RT)=c('Luna', 'Taqpath')
levels(dfs$PositivePatient)=c('Negative', 'Positive')
dfsc = dfs %>% filter(RT=='Taqpath') %>% filter(Stotal>1000) %>% filter(grepl('^R|^O|^C', lysate)) #%>% filter( !(Col %in% c('09','10','11','12')))



fig2d=
    ggplot(dfsc, aes(x=PositivePatient, y=S2_normalized_to_S2_spike, group=virus_copy))+
    geom_quasirandom(alpha=.75)+scale_colour_manual(values=c('blue', 'red'))+
       scale_y_log10() + annotation_logticks() + xlab("")+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+ylab('(S2+1)/(S2 spike+1)')+
    ggtitle('Nasal Swab in Normal Saline, extraction-free') #+geom_hline(yintercept=.003)
#------------------------------------------


rundir=paste0(swabseq.dir, 'runs/v19/')
dfsR=mungeTables(paste0(rundir, 'countTable.RDS'))
dfsR$Plate_384=droplevels(dfsR$Plate_384)
dfsR$virus_copy[dfsR$virus_copy=='Negative Patient']='0'
dfsR$virus_copy=droplevels(dfsR$virus_copy) #Plate_384=droplevels(dfsR$Plate_384)
dfsR$virus_copy=factor(dfsR$virus_copy, levels(dfsR$virus_copy)[order(as.numeric(levels(dfsR$virus_copy)))])
#dfsR$virus_copy=as.character(dfsR$virus_copy)
dfsR$virus_ident2=dfsR$virus_identity
dfsR$virus_ident2[grepl('^N|^P', dfsR$virus_identity)]='Negative'
dfsR$virus_ident2[grepl('^U', dfsR$virus_identity)]='Positive'
dfsR=dfsR[dfsR$SARS_COV_2_Detected!='Inconclusive',]
fig2e=dfsR %>% filter(virus_ident2=='Negative' | virus_ident2=='Positive') %>% filter(SARS_COV_2_Detected!='Inconclusive') %>% 
ggplot(aes(x=virus_ident2, y=S2_normalized_to_S2_spike))+ #, color=SARS_COV_2_Detected))+ #spike+S2, group=virus_copy))+ #log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike+1)')+ #xlab('patient SARs-CoV-2 status')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+ xlab("")+
    theme_bw()+ggtitle('Nasal Swab, Purified RNA') 
    #facet_wrap(~virus_ident2)+ 
    #scale_color_viridis(option = 'plasma')+
    #+    geom_hline(yintercept=.003, color='red') 


#fig 2b version 2, NP in VTM 
rundir=paste0(swabseq.dir, 'runs/v24/')
dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfs=dfL$dfs
fig2b2=dfs %>%  filter(grepl('VTM', Description)) %>%  filter(Stotal>2000) %>% filter(lysate=='VTM 1:4') %>%
    ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+
    geom_quasirandom(alpha=.75, size=2)+
    #facet_wrap(~lysate)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    #xlab('copies per mL of lysate')+
    xlab('GCE per mL')+
    theme_bw()+ggtitle("NP in VTM diluted 1:4, extraction-free")
    #$ggtitle('NP into VTM, diluted 1:4')


# fig 2c version 2, saliva into water 
#rundir=paste0(swabseq.dir, 'runs/v25/')
#dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
#df=dfL$df
#dfs=dfL$dfs
#fig2c2=dfs %>%  filter(Plate_ID=='Plate9') %>%  filter(Stotal>1000) %>% filter(Row != 'A') %>% filter(lysate=='Saliva 1:2') %>% filter( treatment=='TBE+0.5% Tween') %>%
#ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+
#    geom_quasirandom(alpha=.75, size=2)+
#    #facet_wrap(~lysate+treatment)+
#    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+xlab('copies per mL of lysate')+
#    theme_bw()+ggtitle('Saliva, 95C for 15 min; TBE+0.5% Tween; diluted 1:2')
rundir=paste0(swabseq.dir, 'runs/v26/')
dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfs2c=dfL$dfs
dfs2c=   dfs2c %>%  filter(Plate_ID=='Plate4') %>%  filter(Stotal>1000) %>% filter(virus_identity!='ED') %>% filter(virus_copy!='NA')
#%>%
   #filter(Row != 'A') %>%
fig2c2=
ggplot(dfs2c,aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+
    geom_quasirandom(alpha=.75, size=2)+
    #facet_wrap(~lysate+treatment)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    #xlab('copies per mL of lysate')+
    xlab('GCE per mL')+
    theme_bw()+ggtitle('Saliva, extraction-free')

swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v34/')
outdir=paste0(swabseq.dir, 'analysis/v34/')
dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfsM=dfL$dfs
munge=dfsM %>% filter(Plate_ID=='Plate11' | Plate_ID=='Plate2') %>% filter(Stotal>1000)
munge$lysate=gsub( '^3.*' , 'EDneg',munge$lysate)
munge$lysate=droplevels(factor(munge$lysate))
#filter(Col!='06'|Col!='07'|Col!='08') %>%
#'368309416'
munge=munge %>%  filter(lysate=='EDneg' | lysate == 'saliva in 1XTBE+0.5%tw') %>% filter(virus_identity!='368309416')
munge$virus_copy[is.na(munge$virus_copy)]=0

dfs2c2=rbind(dfs2c[,c('virus_copy', 'S2_normalized_to_S2_spike')],munge[,c('virus_copy', 'S2_normalized_to_S2_spike')])

dfs2c2$virus_copy=factor(as.numeric(as.character(dfs2c2$virus_copy)))
  
fig2c3=ggplot(dfs2c2,aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+
    geom_quasirandom(alpha=.75, size=2)+
    #facet_wrap(~lysate+treatment)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    #xlab('copies per mL of lysate')+
    xlab('GCE per mL')+
    theme_bw()+ggtitle('Saliva, extraction-free')









rundir=paste0(swabseq.dir, 'runs/v11/')
dfs=mungeTables(paste0(rundir, 'countTable.RDS'))
dfs=dfs%>%filter(lysate=='MTS_TE-1_to_2')
vc=(as.numeric(as.character((dfs$virus_copy))))
vc.n=2*((1000*(vc/7))/2)
dfs$vc.n=factor(vc.n) #factor(vc.n)
levels(dfs$vc.n)=as.character(round(as.numeric(levels(dfs$vc.n))))

fig2b3=
dfs%>% 
ggplot(aes(x=vc.n, y=S2_normalized_to_S2_spike))+geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 Spike + 1)')+xlab('GCE per mL')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+
    ggtitle('Nasal Swab in TE, extraction-free')







ggarrange(fig2a,fig2b,fig2c,fig2d,fig2e, labels=c('A', 'B', 'C', 'D', 'E'),ncol=3, nrow=2)
ggsave('/data/Covid/swabseq/analysis/manuscript/fig2.png')
ggsave('/data/Covid/swabseq/analysis/manuscript/fig2.svg')

ggarrange(fig2a,fig2b3,fig2c2,fig2e, fig2d, labels=c('A', 'B', 'C', 'D', 'E'),ncol=3, nrow=2)
ggsave('/data/Covid/swabseq/analysis/manuscript/fig2v2.png')
ggsave('/data/Covid/swabseq/analysis/manuscript/fig2v2.svg')

ggarrange(fig2a,fig2e,fig2b3,fig2d,fig2c2, labels=c('A', 'B', 'C', 'D', 'E'),ncol=3, nrow=2)
ggsave('/data/Covid/swabseq/analysis/manuscript/fig2v3.png')

ggarrange(fig2a2,fig2e,fig2b3,fig2d,fig2c3, labels=c('A', 'B', 'C', 'D', 'E'),ncol=3, nrow=2)
ggsave('/data/Covid/swabseq/analysis/manuscript/fig2v4.png')



#========================================================================================================================


# supplementary figures
#------------------ linearity given high virus input -----------------------------------
rundir=paste0(swabseq.dir, 'runs/v17/')
dfsL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
dfs=dfsL$dfs
dfs$Plate_384=droplevels(dfs$Plate_384)
dfs$virus_copy[dfs$virus_copy=='Negative Patient']='0.0'
dfs$virus_copy=as.numeric(as.character(dfs$virus_copy))*1000
dfs$virus_copy=droplevels(as.factor(dfs$virus_copy))
dfs$virus_copy=factor(dfs$virus_copy, levels(dfs$virus_copy)[order(as.numeric(levels(dfs$virus_copy)))])
dfs$virus_ident2=dfs$virus_identity
dfs$virus_ident2[grepl('^N|^P', dfs$virus_identity)]='NegPatient'
dfs$virus_copy=as.numeric(as.character(dfs$virus_copy))
dfs$zeros=as.factor(dfs$virus_copy=='0')
levels(dfs$zeros)=c('>0 Virus Added', 'No Virus Added')
dfs$zeros=relevel(dfs$zeros,'No Virus Added')

s1a=
dfs %>% filter(Plate_ID!='Plate8') %>% filter(Stotal>2000) %>% filter(virus_ident2!='D1' & virus_ident2!='D2') %>% 
ggplot( aes(x=virus_copy+1, y=S2_normalized_to_S2_spike))+geom_quasirandom(alpha=.75)+
    scale_x_log10()+facet_grid(~zeros, scales='free_x', space='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike+1)')+xlab('copies/mL +1')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()

s1b=
dfs %>% filter(Plate_ID!='Plate8') %>% filter(Stotal>2000) %>% filter(virus_ident2!='D1' & virus_ident2!='D2') %>% 
ggplot( aes(x=virus_copy+1, y=S2+1))+geom_quasirandom(alpha=.75)+
    scale_x_log10()+facet_grid(~zeros, scales='free_x', space='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('S2+1')+xlab('copies/mL +1')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()

s1c=
dfs %>% filter(Plate_ID!='Plate8') %>% filter(Stotal>2000) %>% filter(virus_ident2!='D1' & virus_ident2!='D2') %>% 
ggplot( aes(x=virus_copy+1, y=S2_spike+1))+geom_quasirandom(alpha=.75)+
    scale_x_log10()+facet_grid(~zeros, scales='free_x', space='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('S2 spike+1')+xlab('copies/mL +1')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()

s1bc=ggarrange(s1b,s1c, nrow=2, labels=c('A','B'))
s1abc=ggarrange(s1bc, s1a, ncol=2, labels=c('A', 'C'))
ggsave('/data/Covid/swabseq/analysis/manuscript/Sfig_linearity2.svg')
ggsave('/data/Covid/swabseq/analysis/manuscript/Sfig_linearity2.png')
# -------------------------------------------------------------------------------------------



# with and without bleach wash ---------------------------------------------------------------
#14 no bleach
rundir=paste0(swabseq.dir, 'runs/v14/')
df14=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)$df
df14$bleach='No Bleach Wash'
df14=df14 %>%filter(Plate_ID=='Plate9' | Plate_ID=='Plate10'| Plate_ID=='Plate11' |Plate_ID=='Plate12')
df14$Plate_384_Quadrant=''
df14$Plate_384_Quadrant[df14$Plate_ID=='Plate9']='A'
df14$Plate_384_Quadrant[df14$Plate_ID=='Plate10']='B'
df14$Plate_384_Quadrant[df14$Plate_ID=='Plate11']='C'
df14$Plate_384_Quadrant[df14$Plate_ID=='Plate12']='D'
df14$Plate_384='No Bleach Wash'




#15 bleach 
rundir=paste0(swabseq.dir, 'runs/v15/')
df15=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)$df
df15$bleach='With Bleach Wash'
df15=df15 %>%filter(Plate_ID=='Plate13' | Plate_ID=='Plate14'| Plate_ID=='Plate15' |Plate_ID=='Plate16')
df15$Plate_384_Quadrant=''
df15$Plate_384_Quadrant[df15$Plate_ID=='Plate13']='A'
df15$Plate_384_Quadrant[df15$Plate_ID=='Plate14']='B'
df15$Plate_384_Quadrant[df15$Plate_ID=='Plate15']='C'
df15$Plate_384_Quadrant[df15$Plate_ID=='Plate16']='D'
df15$Plate_384='Bleach Wash'
dfb=rbind(df14,df15)

dfb$Plate_384=factor(dfb$Plate_384, levels=c('No Bleach Wash', 'Bleach Wash'))
#relevel(dfb$Plate_384)='No Bleach Wash'
#move this to helper functions
col384=sprintf('%02d', 1:24)
row384=toupper(letters[1:16])
col384L=list(
'A'=col384[seq(1,24,2)],
'B'=col384[seq(2,24,2)],
'C'=col384[seq(1,24,2)],
'D'=col384[seq(2,24,2)])
col384L=lapply(col384L, function(x) { names(x)=sprintf('%02d', 1:12); return(x); })

row384L=list(
'A'=row384[seq(1,16,2)],
'B'=row384[seq(1,16,2)],
'C'=row384[seq(2,16,2)],
'D'=row384[seq(2,16,2)])
row384L=lapply(row384L, function(x) { names(x)=toupper(letters[1:8]); return(x); })

dfb$Row384=''
dfb$Col384=''
dfb$Row384[dfb$Plate_384_Quadrant=='A']=as.character(row384L[['A']][as.character(dfb$Row[dfb$Plate_384_Quadrant=='A'])])
dfb$Row384[dfb$Plate_384_Quadrant=='B']=as.character(row384L[['B']][as.character(dfb$Row[dfb$Plate_384_Quadrant=='B'])])
dfb$Row384[dfb$Plate_384_Quadrant=='C']=as.character(row384L[['C']][as.character(dfb$Row[dfb$Plate_384_Quadrant=='C'])])
dfb$Row384[dfb$Plate_384_Quadrant=='D']=as.character(row384L[['D']][as.character(dfb$Row[dfb$Plate_384_Quadrant=='D'])])
dfb$Col384[dfb$Plate_384_Quadrant=='A']=as.character(col384L[['A']][as.character(dfb$Col[dfb$Plate_384_Quadrant=='A'])])
dfb$Col384[dfb$Plate_384_Quadrant=='B']=as.character(col384L[['B']][as.character(dfb$Col[dfb$Plate_384_Quadrant=='B'])])
dfb$Col384[dfb$Plate_384_Quadrant=='C']=as.character(col384L[['C']][as.character(dfb$Col[dfb$Plate_384_Quadrant=='C'])])
dfb$Col384[dfb$Plate_384_Quadrant=='D']=as.character(col384L[['D']][as.character(dfb$Col[dfb$Plate_384_Quadrant=='D'])])

filled.plates=df %>% filter(Description!='' & Description!=' ') %>% filter(Plate_384!='')
#filled.plates$Row384=droplevels(factor(filled.plates$Row384))
#filled.plates$Col384=droplevels(factor(filled.plates$Col384))

dfb$Row384=factor(dfb$Row384, levels=c(rev(toupper(letters[1:16]))))



















dfb %>%
  #filter(str_detect(amplicon, "spike")) %>%
  ggplot(aes(x=Col384, y=Row384, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384) +
  scale_fill_viridis_c(option = 'plasma')+xlab('Col')+ylab('Row')
ggsave('/data/Covid/swabseq/analysis/manuscript/Sfig_bleachwash.svg')
ggsave('/data/Covid/swabseq/analysis/manuscript/Sfig_bleachwash.png')
#--------------------------------------------------------------------------------------------------------


# luna vs taqpath ------------------------------------------------
rundir=paste0(swabseq.dir, 'runs/v10/')
dfs=mungeTables(paste0(rundir, 'countTable.RDS'))
#Plates1 and Plate2 are taqpath
dfs$RT=as.factor(dfs$Plate_ID %in% c("Plate1", "Plate2"))
dfs$PositivePatient=as.factor(dfs$lysate %in% c('R003', 'O-003', 'R-005','C-009', 'R-014', 'R006', 'O-009', 'R-004', 'O-014', 'C-008'))
levels(dfs$RT)=c('Luna', 'Taqpath')
levels(dfs$PositivePatient)=c('Negative', 'Positive')
dfsc = dfs %>% filter(Stotal>1000) %>% filter(grepl('^R|^O|^C', lysate)) #%>% filter( !(Col %in% c('09','10','11','12')))

slt1=
dfsc %>% filter(Stotal>1000) %>%
    ggplot( aes(x=PositivePatient, y=S2_normalized_to_S2_spike))+geom_quasirandom(alpha=.75)+
   facet_grid(~RT)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike+1)')+xlab('Patient SARS-CoV-2 status')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()
slt2=
dfsc %>% filter(Stotal>1000) %>%
    ggplot( aes(x=PositivePatient, y=S2+1))+geom_quasirandom(alpha=.75)+
   facet_grid(~RT)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)')+xlab('Patient SARS-CoV-2 status')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()

slt=ggarrange(slt2, slt1, nrow=2)
slt=ggarrange(slt2, slt1, nrow=2, labels=c('A' ,'B'))
ggsave('/data/Covid/swabseq/analysis/manuscript/Sfig_neb_v_luna.png')

#-----------------------------------------------------------------



#demonstration of hopping normalization ---------------------------------------------------
rundir=paste0(swabseq.dir, 'runs/v15/')
dfs=mungeTables(paste0(rundir, 'countTable.RDS'))
dfs$lS2=log10(dfs$S2_normalized_to_S2_spike)
library(lme4)
#full.model=lmer(scale(log10(dfs$S2_normalized_to_S2_spike))~scale(log10(well_total))+(1|virus_copy)+(1|Sample_Well)+(1|Plate_ID)+(1|lysate), data=dfs)
#reduced.model=lmer(scale(log10(dfs$S2_normalized_to_S2_spike))~scale(log10(well_total))+(1|virus_copy)+(1|Plate_ID)+(1|lysate), data=dfs)
#summary(full.model)
#anova(full.model, reduced.model)
simpler.model=lmer(scale(log10(dfs$S2_normalized_to_S2_spike))~(1|Sample_Well)+(1|Plate_ID), data=dfs)
dfs$lS2_corrected=residuals(simpler.model) #=log10(dfs$S2_normalized_to_S2_spike)
dfn=dfs %>%
 filter(Plate_ID %in% c('Plate5', 'Plate6', 'Plate7')) %>%
 filter(lysate == 'SimulatedPatient') 

dfn$Plate_384=droplevels(dfn$Plate_384)
dfn$Plate_ID=droplevels(dfn$Plate_ID)
dfn$virus_copy=droplevels(dfn$virus_copy)
dfn$virus_copy=factor(dfn$virus_copy, levels(dfn$virus_copy)[order(as.numeric(levels(dfn$virus_copy)))])
dfn$PositivePatient=dfn$virus_identity %in% c('U0002', 'U0003', 'U0004', 'U0070', 
                                              'U0036', 'U0037', 'U0040', 'U0041', 'U0077', 'U0078')
dfn$Ct=NA
dfn$Ct[dfn$virus_identity=='U0002']=12.1
dfn$Ct[dfn$virus_identity=='U0003']=23
dfn$Ct[dfn$virus_identity=='U0004']=16.2
dfn$Ct[dfn$virus_identity=='U0070']=31
dfn$Ct[dfn$virus_identity=='U0036']=36.4
dfn$Ct[dfn$virus_identity=='U0037']=11
dfn$Ct[dfn$virus_identity=='U0040']=26.2
dfn$Ct[dfn$virus_identity=='U0041']=22.7
dfn$Ct[dfn$virus_identity=='U0077']=32.4
dfn$Ct[dfn$virus_identity=='U0078']=34.7

fn1=dfn %>% filter(Stotal>1000) %>% filter(Ct<32 | is.na(Ct) )%>%
ggplot(aes(x=PositivePatient, y=S2_normalized_to_S2_spike, color=log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75,size=2.5)+
    scale_y_log10() + 
    annotation_logticks() +
    scale_color_viridis(option = 'plasma',name='log10(S2 total per i7)')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+ylab('(S2+1)/(S2 spike+1)')

fn2=dfn %>% filter(Stotal>1000) %>% filter(Ct<32 | is.na(Ct) )%>%
ggplot(aes(x=PositivePatient, y=lS2_corrected, color=log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75,size=2)+
    scale_color_viridis(option = 'plasma',name='log10(S2 total per i7)')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+ylab('Normalized (S2+1/S2 spike+1)') 
sfigfn=ggarrange(fn1, fn2, ncol=2, labels=c('A', 'B'))

library(patchwork)
co=fn1+fn2 +plot_layout(guides='collect')
ggsave('/data/Covid/swabseq/analysis/manuscript/Sfig_mm_correct.png')

#---------------------------------------------------------------------------------





rundir=paste0(swabseq.dir, 'runs/v25/')
dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfs=dfL$dfs
#typo in sample sheet
dfs$Description=gsub( '15', '30',dfs$Description)
#plate visualization 
nh1=df %>% filter(Description!='' & Description!=' ') %>% filter(Plate_ID=='Plate3') %>% # | Plate_ID=='Plate6') %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle('A  no preheat')+theme(legend.position="bottom")

nh2=df %>% filter(Description!='' & Description!=' ') %>% filter(Plate_ID=='Plate9') %>% # | Plate_ID=='Plate16') %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle('B  preheat 95C for 30 min')+theme(legend.position="bottom")


ggarrange(nh1, nh2, ncol=2, labels=c('A' ,'B'))
ggsave(paste0(outdir,'Sfig_heat_treat_saliva.png'))

nh1+nh2 +plot_layout(guides='collect')

ggsave(paste0(outdir,'plateVis_plates_run.png'))















rundir=paste0(swabseq.dir, 'runs/v12/')
dfs=mungeTables(paste0(rundir, 'countTable.RDS'))
dfs$Plate_384=droplevels(dfs$Plate_384)
dfs %>% filter(Stotal>1000) %>% 
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy, color=Plate_ID))+geom_quasirandom(alpha=.75)+
    facet_wrap(~lysate+virus_identity+rpp30)+  scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2_spike+1)')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()


rundir=paste0(swabseq.dir, 'runs/v22/')
dfs=mungeTables(paste0(rundir, 'countTable.RDS'))
dfs$Plate_384=droplevels(dfs$Plate_384)
dfs$virus_copy=factor(dfs$virus_copy, levels(dfs$virus_copy)[order(as.numeric(levels(dfs$virus_copy)))])

dfs %>% filter(Stotal>1000) %>%  filter(Plate_ID=='Plate6') %>% 
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy, color=Plate_ID))+geom_quasirandom(alpha=.75)+
     scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2_spike+1)')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()

dfs %>% filter(Stotal>1000) %>%  filter(Plate_ID=='Plate6') %>% 
ggplot(aes(x=virus_copy, y=S2+1, group=virus_copy, color=Plate_ID))+geom_quasirandom(alpha=.75)+
     scale_y_log10() + annotation_logticks() + ylab('(Rspike)')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()





#hamming distance between S2 and S2 amplicon ( and base quality)
library(ShortRead)
library(stringdist)
rundir='/data/Covid/swabseq/runs/v18/'
fastq_dir  <- paste0(rundir, 'bcls/out/')
nbuffer=3e7
in.fileI1  <- paste0(fastq_dir, 'Undetermined_S0_I1_001.fastq.gz')
in.fileI2  <- paste0(fastq_dir, 'Undetermined_S0_I2_001.fastq.gz')
in.fileR1  <- paste0(fastq_dir, 'Undetermined_S0_R1_001.fastq.gz')
i1 <- FastqStreamer(in.fileI1, nbuffer, readerBlockSize = 1e9, verbose = T)
i2 <- FastqStreamer(in.fileI2, nbuffer, readerBlockSize = 1e9, verbose = T)
r1 <- FastqStreamer(in.fileR1, nbuffer, readerBlockSize = 1e9, verbose = T)
rfq1 <- yield(i1) 
rfq2 <- yield(i2) 
rfq3 <- yield(r1) 
ind1 <- sread(rfq1)
ind2 <- sread(rfq2)
rd1  <- sread(rfq3)

test=quality(rfq3)
test=as(quality(rfq3), "matrix") #test=quality(rfq3)
#rmean=apply(test,2, mean)
#rvar=apply(test,2, var)    
q12=apply(test, 2, function(x) sum(x<12))
#  1 in 16 chance of error  Q12 
amplicons=list(
    S2=      'TATCTTCAACCTAGGACTTTTCTATT',
    S2_spike='ATAGAACAACCTAGGACTTTTCTATT',
    RPP30   ='CGCAGAGCCTTCAGGTCAGAACCCGC'
)
png(file='/data/Covid/swabseq/analysis/manuscript/Sfig_ampdist.png',width=1024, height=512)
barplot(q12/nrow(test)*100, ylab='percentage of reads' ,xlab='read 1 position', main='percentage of bases with Q<12 (base call accuracy <92%)', names.arg=1:26)
dev.off()

d1=stringdist(amplicons[[1]], rd1, method='hamming')
d2=stringdist(amplicons[[2]], rd1, method='hamming')
d3=stringdist(amplicons[[3]], rd1, method='hamming')
distance_from_S2=d1
distance_from_S2_Spike=d2
d4=table(distance_from_S2,distance_from_S2_Spike)
hdist.table=d4[1:8,1:8]
write.table(hdist.table, file='/data/Covid/swabseq/analysis/manuscript/Sfig_ampdist_table.csv', row.names=T, col.names=T, quote=F)
#library(gridExtra)
#library(grid)
#par(mfrow=c(2,1))
#plot(0,0, type='n')
#ad2=grid.table(hdist.table)
  


#sup fig 14 VTM, Amies, Normal Saline

rundir=paste0(swabseq.dir, 'runs/v24/')
dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfs=dfL$dfs
fig14a=dfs %>%  filter(grepl('VTM', Description)) %>%  filter(Stotal>2000) %>% filter(lysate=='VTM 1:4') %>%
    ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+
    geom_quasirandom(alpha=.75, size=2)+
    #facet_wrap(~lysate)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    #xlab('copies per mL of lysate')+
    xlab('GCE per mL')+
    theme_bw()+ggtitle("NP in VTM diluted 1:4, extraction-free")
    #$ggtitle('NP into VTM, diluted 1:4')
fig14b=dfs %>%  filter(lysate == "NS (1:4) ") %>%  filter(Stotal>2000) %>%
    ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+
    geom_quasirandom(alpha=.75, size=2)+
    #facet_wrap(~lysate)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    #xlab('copies per mL of lysate')+
    xlab('GCE per mL')+
    theme_bw()+ggtitle("NP in NS diluted 1:4, extraction-free")
    #$ggtitle('NP into VTM, diluted 1:4')



library(gdata)
swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v32/')
outdir=paste0(swabseq.dir, 'analysis/v32/')

dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfs=dfL$dfs
amies=dfs %>%   filter(Stotal>1000) %>% filter(Plate_ID=='Plate1') 

amies.ct=read.xls('/data/Covid/swabseq/2020_0706 DeIdentified POS COVID Samples.xlsx')
amies.ct[,1]=gsub('VA','',amies.ct[,1])
amies.ct[,1]=gsub('^0','',amies.ct[,1])
amies.ct[,1]=paste0('Aimes ', amies.ct[,1])

names(amies.ct)[1]='virus_identity'

aa=right_join(amies, amies.ct)
aa$virus_identity=gsub('Aimes ', 'Patient ', aa$virus_identity) 


fig14c=aa %>%  filter(Row!='D') %>%  ggplot(aes(x=S.Gene.Ct.Value, y=S2_normalized_to_S2_spike, color=virus_identity))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    #facet_wrap(~Description, scales='free_x')+
    scale_x_reverse()+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+ xlab('S Gene Ct')+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
     ggtitle("Eswab into Amies, diluted 1:10, extraction-free")

ggarrange(fig14a,fig14b,fig14c, labels=c('A', 'B', 'C'),ncol=3, nrow=1)
ggsave('/data/Covid/swabseq/analysis/manuscript/figS_nswab_media.png')




# for comparison with patient data ---------------------------
swabseq.dir="/mnt/e/swabseq/"
swabseq.dir="/data/Covid/swabseq/"

source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v40/')
outdir=paste0(swabseq.dir, 'analysis/v40/')

dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T, Stotal_filt=500)
df=dfL$df
dfs=dfL$dfs
#usable swabseq reads
#sum(df$Count)
#17477902

titl='v40 - ED Study + NS 1:5 Clinical Samples'

#plate visualization 
#plate visualization 
df %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle(titl)
ggsave(paste0(outdir,'plateVis_all_indices.png'))

filled.plates=add384Mapping(df)
filled.plates %>% 
    ggplot(aes(x=Col384, y=Row384, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle(titl)
ggsave(paste0(outdir,'plateVis_384.png'))


df %>% filter(Description!='' & Description!=' ') %>% #filter(Plate_ID!='Plate1') %>%
    ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle(titl)
#ggtitle('v30 Saliva; Nasal; ED; Ashe')
ggsave(paste0(outdir,'plateVis_plates_run.png'))




ed.table=read.delim('/data/Covid/clinic/ED_patients.csv', stringsAsFactors=F, header=T, sep=',')
edo=dfs[as.character(dfs$virus_identity) %in% ed.table$Matrix.Tube.Barcode | as.character(dfs$virus_identity) %in% ed.table$Deidenfied.ID,]
mu=match(edo$lysate, ed.table$Matrix.Tube.Barcode)
edo$virus_identity[!is.na(mu)]=na.omit(ed.table$Deidenfied.ID[mu])
edo=merge(edo, ed.table, all.x=T, by.x='virus_identity', by.y='Deidenfied.ID')

#all.ed=edo[,c(1,2,3,4,5,6,8,26:36,40,52,53)]
#all.ed %>% write.csv(paste0(outdir, 'ED_rerun.csv'))

library(data.table)
edo=rbindlist(lapply(split(edo, edo$virus_identity), function(x) x[which.max(x$Stotal),]))

edo$Ct[is.na(edo$Ct)]=0
edo=edo[1:63,]

edo$NPResult[edo$NPResult=='ND']='Negative'
edo$NPResult[edo$NPResult=='DETECTED']='Positive'


gt=edo%>%filter(Stotal>500)%>%
ggplot(aes(x=Ct, y=S2_normalized_to_S2_spike))+
geom_quasirandom(size=2)+
    scale_y_log10() + annotation_logticks(sides="l")+
    geom_hline(yintercept=3e-3, color='red')+theme_bw()+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    facet_grid(~NPResult, scales='free_x')

gt1=ggplot_gtable(ggplot_build(gt))
gt1$widths[5]=.3*gt1$widths[5]
grid.draw(gt1)



swabseq.dir="/mnt/e/swabseq/"
swabseq.dir="/data/Covid/swabseq/"
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v41/')
outdir=paste0(swabseq.dir, 'analysis/v41/')
dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T, Stotal_filt=500)
df=dfL$df
dfs=dfL$dfs
#usable swabseq reads
#sum(df$Count)
#17477902
library(gdata)
pd=read.xls(paste0(swabseq.dir,'analysis/v37/', 'EUA_data.xlsx'), stringsAsFactors=F)
pd$Sample[2]='R-003'
pd$Sample[19]='R-006'
dfsF =dfs %>%  filter(Plate_ID=='Plate15')
dfsF=merge(dfsF, pd, by.x='virus_identity', by.y='Sample')
names(dfsF)[41]='CovidDetectedNP'


library(grid)

names(dfsF)[38]='Ct'
names(dfsF)[41]='NPResult'
gt=ggplot(dfsF, aes(x=Ct, y=S2_normalized_to_S2_spike))+
geom_quasirandom(size=2)+
    scale_y_log10() + annotation_logticks(sides="l")+
    geom_hline(yintercept=3e-3, color='red')+theme_bw()+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    facet_grid(~NPResult, scales='free_x')

gt2=ggplot_gtable(ggplot_build(gt))
gt2$widths[5]=.3*gt2$widths[5]
#grid.draw(gt2)

library(ggpubr)
ggarrange(gt1,gt2, nrow=2, labels=c('A','B'))
ggsave(file='/data/Covid/swabseq/analysis/manuscript/patient_comparison.png')

ggarrange(gt1)
ggsave(file='/data/Covid/swabseq/analysis/manuscript/patient_comparison_saliva.png')

ggarrange(gt2)
ggsave(file='/data/Covid/swabseq/analysis/manuscript/patient_comparison_NS.png')
