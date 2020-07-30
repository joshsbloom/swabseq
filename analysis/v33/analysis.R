 #v23
swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v33/')
outdir=paste0(swabseq.dir, 'analysis/v33/')

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
  scale_fill_viridis_c(option = 'plasma')+ggtitle('v33, LoD, contam test and ED;')
ggsave(paste0(outdir,'plateVis_all_indices.png'))

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

df$Row384=''
df$Col384=''
df$Row384[df$Plate_384_Quadrant=='A']=as.character(row384L[['A']][as.character(df$Row[df$Plate_384_Quadrant=='A'])])
df$Row384[df$Plate_384_Quadrant=='B']=as.character(row384L[['B']][as.character(df$Row[df$Plate_384_Quadrant=='B'])])
df$Row384[df$Plate_384_Quadrant=='C']=as.character(row384L[['C']][as.character(df$Row[df$Plate_384_Quadrant=='C'])])
df$Row384[df$Plate_384_Quadrant=='D']=as.character(row384L[['D']][as.character(df$Row[df$Plate_384_Quadrant=='D'])])
df$Col384[df$Plate_384_Quadrant=='A']=as.character(col384L[['A']][as.character(df$Col[df$Plate_384_Quadrant=='A'])])
df$Col384[df$Plate_384_Quadrant=='B']=as.character(col384L[['B']][as.character(df$Col[df$Plate_384_Quadrant=='B'])])
df$Col384[df$Plate_384_Quadrant=='C']=as.character(col384L[['C']][as.character(df$Col[df$Plate_384_Quadrant=='C'])])
df$Col384[df$Plate_384_Quadrant=='D']=as.character(col384L[['D']][as.character(df$Col[df$Plate_384_Quadrant=='D'])])

filled.plates=df %>% filter(Description!='' & Description!=' ') %>% filter(Plate_384!='')
#filled.plates$Row384=droplevels(factor(filled.plates$Row384))
#filled.plates$Col384=droplevels(factor(filled.plates$Col384))

filled.plates$Row384=factor(filled.plates$Row384, levels=c(rev(toupper(letters[1:16]))))
#filled.plates$Col384=factor(filled.plates$Col384, levels=c(sprintf('%2d', 1:24))) #rev(toupper(letters[1:16])))


filled.plates %>% 
    ggplot(aes(x=Col384, y=Row384, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle('v33, LoD, contam test and ED;')
ggsave(paste0(outdir,'plateVis_384.png'))


df %>% filter(Description!='' & Description!=' ') %>% 
    ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle('v33, LoD, contam test and ED;')
#ggtitle('v30 Saliva; Nasal; ED; Ashe')
ggsave(paste0(outdir,'plateVis_plates_run.png'))


#saliva
dfs %>% filter(Plate_ID=='Plate10') %>% filter(Col!='06') %>% filter(Col!='07')%>% filter(Col!='08') %>%
    filter(Stotal>1000) %>%
    ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate, scales='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(paste0(outdir,'XY_saliva_confirmatory_LoD.png'))


dfs %>% filter(Plate_ID=='Plate5') %>% #filter(Col!='06'|Col!='07'|Col!='08') %>%
    filter(Stotal>1000) %>%
    ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate, scales='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(paste0(outdir,'XY_MNS_prelim_LoD.png'))




library(gdata)

amies=dfs %>%   filter(Stotal>1000) %>% filter(Plate_ID=='Plate1') 

amies.ct=read.xls('/data/Covid/swabseq/2020_0706 DeIdentified POS COVID Samples.xlsx')
amies.ct[,1]=gsub('VA','',amies.ct[,1])
amies.ct[,1]=gsub('^0','',amies.ct[,1])
amies.ct[,1]=paste0('Aimes ', amies.ct[,1])

names(amies.ct)[1]='virus_identity'

aa=right_join(amies, amies.ct)
aa$virus_identity=gsub('Aimes ', 'Patient ', aa$virus_identity) 


aa %>%  filter(Row!='D') %>%  ggplot(aes(x=S.Gene.Ct.Value, y=S2_normalized_to_S2_spike, color=virus_identity))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    #facet_wrap(~Description, scales='free_x')+
    scale_x_reverse()+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))











library(ggpubr)
df %>% filter(Description!='' & Description!=' ') %>% 
    ggplot(aes(x=Col, y=Row, fill=log10(as.numeric(as.character(virus_copy))))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')+
  theme(legend.title = element_blank())+
  ggtitle('ATCC Heat Inactivated virus added (copies/mL)')
ggsave(paste0(outdir,'ExperimentSetup_ViralDilution.png'))

df %>% filter(Description!='' & Description!=' ') %>%
  ggplot(aes(x=Col, y=Row, fill=lysate) )+ 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  ggtitle('Experiment Setup - lysate')
ggsave(paste0(outdir,'ExperimentSetup_Lysate.png'))

df %>% filter(Description!='' & Description!=' ') %>%
  ggplot(aes(x=Col, y=Row, fill=virus_identity) )+ 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  ggtitle('Experiment Setup - Identity')
ggsave(paste0(outdir,'ExperimentSetup_treatment.png'))



dfs %>%  filter(Description!='' & Description!=' ') %>%
  ggplot(aes(x=Col, y=Row,fill=SARS_COV_2_Detected)) + 
  geom_raster() +
  coord_equal() +
  scale_fill_manual(values=c("lightblue", "white", "red"))+
  facet_grid(~Plate_384+Plate_ID+Description) 
ggsave(paste0(outdir,'plateVis_CovidDetected.png'))


dfs %>%  filter(Description!='' & Description!=' ') %>%
ggplot(aes(x=Col, y=Row,fill=Stotal>1000)) + 
  geom_raster() +
  coord_equal() +
  #scale_fill_manual(values=c("lightblue", "white", "red"))+
  facet_grid(~Plate_384+Plate_ID+Description) 
ggsave(paste0(outdir,'plateVis_Stotal_gt_1000.png'))

#dfs$lysate[dfs$lysate=='Saliva 1:1'] = 'saliva'

dfs %>%   filter(Stotal>1000) %>% filter(Plate_ID=='Plate2' | Plate_ID=='Plate12') %>% 
    #$filter(virus_identity=='Contrived' & lysate=='Saliva') %>%
    ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate, scales='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(paste0(outdir,'Amies_and_MNS.png'))
 

dfs %>%   filter(Stotal>1000) %>% filter(Plate_ID=='Plate4' | Plate_ID=='Plate5' | Plate_ID=='Plate10') %>% 
    ggplot(aes(x=virus_identity, y=S2_normalized_to_S2_spike))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~Description, scales='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(paste0(outdir,'Contam.png'))

dfs %>%   filter(Stotal>1000) %>% filter(Plate_ID=='Plate4' | Plate_ID=='Plate5' | Plate_ID=='Plate10') %>% 
    ggplot(aes(x=virus_identity, y=S2+1, color=log10(S2_total_across_all_wells)))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~Description, scales='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(paste0(outdir,'Contam_S2.png'))
library(lme4)



r=residuals(lmer(log10(S2_normalized_to_S2_spike)~(1|S_index)+(1|S_index2)+(1|Plate_ID), data=dfs))
dfs$r=r
dfs %>%   filter(Stotal>1000) %>% filter(Plate_ID=='Plate4' | Plate_ID=='Plate5' | Plate_ID=='Plate10') %>% 
    ggplot(aes(x=virus_identity, y=r)) + #, color=log10(S2_total_across_all_wells)))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~Description, scales='free_x')+
    #scale_y_log10() 
    annotation_logticks() + ylab('normalized')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('log10(S2+1/Spike+1)=a*i7+b*i5+c*plate+e')

dfs %>%   filter(Stotal>1000) %>% filter(Plate_ID!='Plate5')  %>% filter(lysate=='Saliva') %>% 
    filter(!grepl('^C', virus_identity))%>% 
     filter(!grepl('^S',virus_identity)) %>%
    group_by(virus_identity)%>% 
    mutate(two_hits=sum(S2_normalized_to_S2_spike>1e-2)>1) %>% 
ggplot(aes(x=virus_identity, y=S2_normalized_to_S2_spike, color=two_hits))+
    geom_point(alpha=.75, size=2)+
     #geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate, scales='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(paste0(outdir,'Ashe_Saliva.png'))

368286765











 

dfs %>%  filter(Plate_ID=='Plate2' | Plate_ID=='Plate5') %>% 
    filter(Stotal>1000) %>% filter(lysate=='saliva') %>% 
       ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike))+
    #geom_point(alpha=.75, size=2)+
     geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate+treatment+Plate_ID, scales='free_x')+
    scale_y_log10() +  ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('Contrived for EUA')
ggsave(paste0(outdir,'Contrived.png'))

dfs %>%  filter(Plate_ID=='Plate2' | Plate_ID=='Plate5') %>% 
    filter(Stotal>2000) %>% filter(lysate=='saliva') %>% 
       ggplot(aes(x=virus_copy, y=S2+1))+
    #geom_point(alpha=.75, size=2)+
     geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate+treatment+Plate_ID, scales='free_x')+
    scale_y_log10() +  ylab('(S2+1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('Contrived for EUA')
ggsave(paste0(outdir,'Contrived_S2.png'))





 +geom_hline(yintercept=1e-2, color='red')+ggtitle('Plates 2 and 3')

 ggsave(paste0(outdir,'Ashe_ED_plates.png'))




dfs %>%  filter(Stotal>1000) %>% 
    ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike))+
    geom_point(alpha=.75, size=2)+
     #geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate+treatment+virus_identity, scales='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

 
 +geom_hline(yintercept=1e-2, color='red')+ggtitle('Plates 2 and 3')

dfs %>%  filter(Plate_ID=='Plate6' | Plate_ID=='Plate11') %>%  filter(Stotal>1000) %>% 
    ggplot(aes(x=virus_identity, y=S2))+
    geom_point(alpha=.75, size=2)+
     #geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate+treatment, scales='free_x')+
    scale_y_log10() +  ylab('(S2)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
     ggtitle('Plates 2 and 3')
ggsave(paste0(outdir,'Ashe_ED_plates_S2.png'))


dfs %>%  filter(Plate_ID=='Plate6' | Plate_ID=='Plate11') %>%  filter(Stotal>1000) %>% 
    ggplot(aes(x=virus_identity, y=S2_spike))+
    geom_point(alpha=.75, size=2)+
     #geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate+treatment, scales='free_x')+
    scale_y_log10() +  ylab('(S2_spike)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
     ggtitle('Plates 2 and 3')

dfs %>%  filter(Plate_ID=='Plate6' | Plate_ID=='Plate11') %>%  
    ggplot(aes(x=virus_identity, y=Stotal+1))+
    geom_point(alpha=.75, size=2)+

    facet_wrap(~lysate+treatment, scales='free_x')+
    scale_y_log10() + ylab('(S2 + S2 spike)')+xlab('identity')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+geom_hline(yintercept=1e-2, color='red')+ggtitle('Plates 2 and 3')
ggsave(paste0(outdir,'Ashe_ED_plates_dropout.png'))




dfs %>%  filter(Description!='' & Description!=' ') %>% filter(Stotal>1000) %>%
  ggplot(aes(x=Col, y=Row,fill=S2_normalized_to_S2_spike>1e-2)) + 
  geom_raster() +
  geom_text(aes(label=as.factor(virus_identity)))+
  coord_equal() +
  #scale_fill_manual(values=c("lightblue", "white", "red"))+
  facet_grid(~Plate_384+Plate_ID+Description) 

dfs %>%  filter(Description!='' & Description!=' ') %>%
ggplot(aes(x=Col, y=Row,fill=SARS_COV_2_Detected)) + 
  geom_raster() +
  coord_equal() +
  #scale_fill_manual(values=c("lightblue", "white", "red"))+
  facet_grid(~Plate_384+Plate_ID+Description) 



#remove row A from plate 9 (and 3) due to pipetting error 
dfs %>%  filter(Plate_ID=='Plate4') %>%  filter(Stotal>1000) %>% 
   #filter(Row != 'A') %>%
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy, col=virus_identity=='ED' & (Row=='F' | Row=='H')))+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate+treatment+virus_identity)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+xlab('copies per mL of lysate')+
    theme_bw()+ggtitle('Saliva; 2X TBE; 95C 30 min')
ggsave(paste0(outdir,'XY_Saliva_TBE.png'))



#remove row A from plate 9 (and 3) due to pipetting error 
dfs %>%  filter(Plate_ID=='Plate16') %>%  filter(Stotal>1000) %>%
    ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate+treatment)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+xlab('copies per mL of lysate')+
    theme_bw()+ggtitle('Saliva; TE,RNAsec,w/woQP; 95C 15 min')
ggsave(paste0(outdir,'XY_SalivaTE RNAsec w_woQP 95C 15 min.png'))




dfs %>%  filter(Plate_ID=='Plate9') %>%  filter(Stotal>1000) %>% filter(Row != 'A') %>%
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+xlab('copies per mL of lysate')+
    theme_bw()+ggtitle('Saliva; TBE w/wo tween; 95C 15 min')

dfs %>%  filter(Plate_ID=='Plate4') %>%  filter(Stotal>200) %>%
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+xlab('copies per mL of lysate')+
    theme_bw()+ggtitle('Saliva + RNA Secure')
ggsave(paste0(outdir,'XY_saliva.png'))


dfs %>%  filter(grepl('VTM', Description)) %>%  filter(Stotal>1000) %>%
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+xlab('copies per mL of lysate')+
    theme_bw()+ggtitle('VTM')
ggsave(paste0(outdir,'XY_VTM.png'))


dfs %>%  filter(grepl('NS', Description)) %>%  filter(Stotal>1000) %>%
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+xlab('copies per mL of lysate')+
    theme_bw()+ggtitle('NS')
ggsave(paste0(outdir,'XY_NS.png'))


dfs %>%  filter(grepl('Amies', Description)) %>%  filter(Stotal>1000) %>%
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+xlab('copies per mL of lysate')+
    theme_bw()+ggtitle('Amies')
ggsave(paste0(outdir,'XY_Amies.png'))









dfs %>%  filter(grepl('Nasal', Description)) %>%  filter(Stotal>2000) %>%
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy, color=lysate))+
    geom_quasirandom(alpha=.75)+
    facet_wrap(~lysate)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/S2 spike + 1)')+xlab('copies per mL of lysate')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
    #scale_color_viridis(option = 'plasma')+
    theme_bw()
ggsave(paste0(outdir,'XY_Nasal_stotal_gt_1000.png'))

dfs %>%  filter(grepl('Saliva', Description)) %>%  filter(Stotal>1000) %>%
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy, color=lysate))+
    geom_quasirandom(alpha=.75)+
    facet_wrap(~lysate)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/S2 spike + 1)')+xlab('copies per mL of lysate')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
    #scale_color_viridis(option = 'plasma')+
    theme_bw()
ggsave(paste0(outdir,'saliva_c2cdna_and_te_gt_1000.png'))

