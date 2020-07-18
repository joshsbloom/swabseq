 #v23
swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v28/')
outdir=paste0(swabseq.dir, 'analysis/v28/')

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
  scale_fill_viridis_c(option = 'plasma')+ggtitle('v28 Saliva; ED; Ashe')
ggsave(paste0(outdir,'plateVis_all_indices.png'))

df %>% filter(Description!='' & Description!=' ') %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle('v28 Saliva; ED; Ashe')
ggsave(paste0(outdir,'plateVis_plates_run.png'))

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

dfs$lysate[dfs$lysate=='Saliva 1:1'] = 'saliva'

dfs %>%   filter(Stotal>1000) %>% 
    ggplot(aes(x=virus_identity, y=S2_normalized_to_S2_spike))+
    geom_point(alpha=.75, size=2)+
     #geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate+treatment, scales='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(paste0(outdir,'Ashe_ED_1.png'))
 
 

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



dfs %>%   filter(Stotal>1000) %>% 
    ggplot(aes(x=virus_identity, y=S2_normalized_to_S2_spike))+
    geom_point(alpha=.75, size=2)+
     #geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate+treatment, scales='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


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

