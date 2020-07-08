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
dfs%>%filter(lysate=='MTS_TE-1_to_2') %>% 
ggplot(aes(x=virus_copy, y=S2+1))+geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)')+xlab('GCE per reaction')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()
fig1e2=
dfs%>%filter(lysate=='MTS_TE-1_to_2') %>% 
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
dfsR %>% filter(virus_ident2!='PosPatient') %>% filter(Stotal>2000) %>%
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+ #log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike+1)')+xlab('GCE per mL')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+ ggtitle('LoD - Purified RNA')+
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
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+ylab('(S2+1)/(S2 spike+1)')+ggtitle('NP in Normal Saline') #+geom_hline(yintercept=.003)
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
    theme_bw()+ggtitle('Purifed Extract') 
    #facet_wrap(~virus_ident2)+ 
    #scale_color_viridis(option = 'plasma')+
    #+    geom_hline(yintercept=.003, color='red') 

ggarrange(fig2a,fig2b,fig2c,fig2d,fig2e, labels=c('A', 'B', 'C', 'D', 'E'),ncol=3, nrow=2)
ggsave('/data/Covid/swabseq/analysis/manuscript/fig2.png')
ggsave('/data/Covid/swabseq/analysis/manuscript/fig2.svg')


#========================================================================================================================


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


#v23
swabseq.dir='/data/Covid/swabseq/'

rundir=paste0(swabseq.dir, 'runs/v23/')
outdir=paste0(swabseq.dir, 'analysis/v23/')

dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfs=dfL$dfs

#plate visualization 
#plate visualization 
df %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')
ggsave(paste0(outdir,'plateVis_all_indices.png'))

df %>% filter(Description!='' & Description!=' ') %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')
ggsave(paste0(outdir,'plateVis_plates_run.png'))

