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
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike+1)')+xlab('log10(copies/mL +1)')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()

s1b=
dfs %>% filter(Plate_ID!='Plate8') %>% filter(Stotal>2000) %>% filter(virus_ident2!='D1' & virus_ident2!='D2') %>% 
ggplot( aes(x=virus_copy+1, y=S2+1))+geom_quasirandom(alpha=.75)+
    scale_x_log10()+facet_grid(~zeros, scales='free_x', space='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('S2+1')+xlab('log10(copies/mL +1)')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()

s1c=
dfs %>% filter(Plate_ID!='Plate8') %>% filter(Stotal>2000) %>% filter(virus_ident2!='D1' & virus_ident2!='D2') %>% 
ggplot( aes(x=virus_copy+1, y=S2_spike+1))+geom_quasirandom(alpha=.75)+
    scale_x_log10()+facet_grid(~zeros, scales='free_x', space='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('S2 spike+1')+xlab('log10(copies/mL +1)')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()

s1bc=ggarrange(s1b,s1c, nrow=2)
s1abc=ggarrange(s1bc, s1a, ncol=2)
ggsave('/data/Covid/swabseq/analysis/manuscript/Sfig_linearity2.svg')
ggsave('/data/Covid/swabseq/analysis/manuscript/Sfig_linearity2.png')
# -------------------------------------------------------------------------------------------



# with and without bleach wash ---------------------------------------------------------------
#14 no bleach
rundir=paste0(swabseq.dir, 'runs/v14/')
df14=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)$df
df14$bleach='No Bleach Wash'
df14=df14 %>%filter(Plate_ID=='Plate9' | Plate_ID=='Plate10'| Plate_ID=='Plate11' |Plate_ID=='Plate12')
#15 bleach 
rundir=paste0(swabseq.dir, 'runs/v15/')
df15=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)$df
df15$bleach='With Bleach Wash'
df15=df15 %>%filter(Plate_ID=='Plate13' | Plate_ID=='Plate14'| Plate_ID=='Plate15' |Plate_ID=='Plate16')
dfb=rbind(df14,df15)

dfb %>%
  #filter(str_detect(amplicon, "spike")) %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~bleach+Plate_384+Plate_ID) +
  scale_fill_viridis_c(option = 'plasma')
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

slt=ggarrange(slt1, slt2, nrow=2)
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



