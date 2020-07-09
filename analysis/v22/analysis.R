swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v22/')
outdir=paste0(swabseq.dir, 'analysis/v22/')

dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfs=dfL$dfs

#plate visualization 
df %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID) +
  scale_fill_viridis_c(option = 'plasma')
ggsave(paste0(outdir,'plateVis.png'))

#subset of plates run
df %>% filter(Description!= ' ') %>% 
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')
ggsave(paste0(outdir,'plateVis_NonEmpty.png'))


#assay results
dfs= df %>%filter(amplicon=='S2') %>%  
  count(Sample_Well, wt=Count, name='S2_total_across_all_wells') %>%
  right_join(df)
dfs= df %>%filter(amplicon=='S2'|amplicon=='S2_spike') %>%  
  count(Sample, wt=Count, name='Stotal') %>%
  right_join(dfs)
dfs= dfs %>% count(Sample, wt=Count, name='well_total') %>%
  right_join(dfs) %>% 
  select(-mergedIndex, -Sample_ID, -index, -index2 ) %>% 
  spread(amplicon, Count) %>% 
  mutate(S2_normalized_to_S2_spike=(S2+1)/(S2_spike+1))%>%
  mutate(RPP30_Detected=RPP30>10) %>%  
  mutate(SARS_COV_2_Detected=S2_normalized_to_S2_spike>.003)
dfs$SARS_COV_2_Detected[!dfs$RPP30_Detected]='Inconclusive'
dfs$SARS_COV_2_Detected[dfs$Stotal<2000]='Inconclusive'

dfs %>% filter(Description!= ' ')  %>% 
  ggplot(aes(x=Col, y=Row,fill=SARS_COV_2_Detected)) + 
  geom_raster() +
  coord_equal() +
  scale_fill_manual(values=c("lightblue", "white", "red"))+
  facet_grid(~Plate_384+Plate_ID+Description) 
ggsave(paste0(outdir,'plateVis_Calls.png'))

dfs %>% filter(Description!= ' ')  %>% 
  ggplot(aes(x=Col, y=Row,fill=Stotal>2000)) + 
  geom_raster() +
  coord_equal() +
  #scale_fill_manual(values=c("lightblue", "white", "red"))+
  facet_grid(~Plate_384+Plate_ID+Description) 
ggsave(paste0(outdir,'plateVis_Stotal_gt_2000.png'))

dfs %>% filter(Description!= ' ')  %>% 
  ggplot(aes(x=Col, y=Row,fill=treatment)) + 
  geom_raster() +
  coord_equal() +
  #scale_fill_manual(values=c("lightblue", "white", "red"))+
  facet_grid(~Plate_384+Plate_ID+Description) 
ggsave(paste0(outdir,'plateVis_treatment.png'))

dfs %>% filter(Description!= ' ')  %>% 
  ggplot(aes(x=Col, y=Row,fill=spin)) + 
  geom_raster() +
  coord_equal() +
  #scale_fill_manual(values=c("lightblue", "white", "red"))+
  facet_grid(~Plate_384+Plate_ID+Description)
ggsave(paste0(outdir,'plateVis_spin.png'))

dfs %>% filter(Description!= ' ')  %>% 
  ggplot(aes(x=Col, y=Row,fill=lysate)) + 
  geom_raster() +
  coord_equal() +
  #scale_fill_manual(values=c("lightblue", "white", "red"))+
  facet_grid(~Plate_384+Plate_ID+Description) 
ggsave(paste0(outdir,'plateVis_lysate.png'))



## in progress ....
#
#
#
#dfs %>%
#  write_csv(paste0(rundir, 'Calls_per_sample.csv')) 
#
#
#
#dfF=dfs %>%  filter(Plate_ID=='Plate16') %>% filter(Stotal>1000)
#dfF$virus_copy=droplevels(dfF$virus_copy)
#dfF$virus_copy=factor(dfF$virus_copy, levels(dfF$virus_copy)[order(as.numeric(levels(dfF$virus_copy)))])
#dfF$sample_type=gsub('_1', '', dfF$lysate) #sapply(strsplit(dfF$lysate, '_'), function(x) x[1])
#dfF$sample_type=gsub('_2', '', dfF$sample_type) #sapply(strsplit(dfF$lysate, '_'), function(x) x[1])
#
#
#
##plots slyaste
#ggplot(dfF, aes(x=virus_copy, y=S2+1, group=virus_copy, color=lysate))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~sample_type)+
#    scale_y_log10() + annotation_logticks() + ylab('(S2+1)')+xlab('copies per rxn')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()
#ggsave(paste0(rundir,'Slysate_S2.png'))
#ggplot(dfF, aes(x=virus_copy, y=S2_spike+1, group=virus_copy, color=lysate))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~sample_type)+
#    scale_y_log10() + annotation_logticks() + ylab('(S2_spike+1)')+xlab('copies per rxn')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()
#ggsave(paste0(rundir,'Slysate_S2_spike.png'))
#ggplot(dfF, aes(x=virus_copy, y=RPP30+1, group=virus_copy, color=lysate))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~sample_type)+
#    scale_y_log10() + annotation_logticks() + ylab('(RPP30+1)')+xlab('copies per rxn')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()
#ggsave(paste0(rundir,'Slysate_RPP30_spike.png'))
#
#ggplot(dfF, aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy, color=lysate))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~sample_type)+
#    scale_y_log10() + annotation_logticks() + ylab('(S2+1/S2_spike+1)')+xlab('copies per rxn')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()
#ggsave(paste0(rundir,'Slysate_S2nrom.png'))
#
#
#
#
#
#
#
#
#dfF=dfs %>%  filter(Plate_ID=='Plate16') %>% filter(Stotal>1000)
#dfF$virus_copy=droplevels(dfF$virus_copy)
#dfF$virus_copy=factor(dfF$virus_copy, levels(dfF$virus_copy)[order(as.numeric(levels(dfF$virus_copy)))])
#dfF$sample_type=gsub('_1', '', dfF$lysate) #sapply(strsplit(dfF$lysate, '_'), function(x) x[1])
#dfF$sample_type=gsub('_2', '', dfF$sample_type) #sapply(strsplit(dfF$lysate, '_'), function(x) x[1])
#
#
#
##plots slyaste
#ggplot(dfF, aes(x=virus_copy, y=S2+1, group=virus_copy, color=lysate))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~sample_type)+
#    scale_y_log10() + annotation_logticks() + ylab('(S2+1)')+xlab('copies per rxn')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()
#ggsave(paste0(rundir,'Slysate_S2.png'))
#ggplot(dfF, aes(x=virus_copy, y=S2_spike+1, group=virus_copy, color=lysate))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~sample_type)+
#    scale_y_log10() + annotation_logticks() + ylab('(S2_spike+1)')+xlab('copies per rxn')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()
#ggsave(paste0(rundir,'Slysate_S2_spike.png'))
#ggplot(dfF, aes(x=virus_copy, y=RPP30+1, group=virus_copy, color=lysate))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~sample_type)+
#    scale_y_log10() + annotation_logticks() + ylab('(RPP30+1)')+xlab('copies per rxn')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()
#ggsave(paste0(rundir,'Slysate_RPP30_spike.png'))
#
#ggplot(dfF, aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy, color=lysate))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~sample_type)+
#    scale_y_log10() + annotation_logticks() + ylab('(S2+1/S2_spike+1)')+xlab('copies per rxn')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()
#ggsave(paste0(rundir,'Slysate_S2nrom.png'))
#
#
#######################3
#dfF=dfs %>%  filter(Plate_ID=='Plate1') #%>%
#   # filter(Stotal>1000)
#dfF$virus_copy=droplevels(dfF$virus_copy)
#dfF$virus_copy=factor(dfF$virus_copy, levels(dfF$virus_copy)[order(as.numeric(levels(dfF$virus_copy)))])
##dfF$sample_type=gsub('_1', '', dfF$lysate) #sapply(strsplit(dfF$lysate, '_'), function(x) x[1])
##dfF$sample_type=gsub('_2', '', dfF$sample_type) #sapply(strsplit(dfF$lysate, '_'), function(x) x[1])
#dfF$treatment=gsub(' only' ,'' ,dfF$treatment)
#dfF$spin=gsub('spin down' ,'Spin down' ,dfF$spin)
#
#
##plots saliva
#ggplot(dfF, aes(x=virus_copy, y=S2+1, group=virus_copy, color=Stotal>1000))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~treatment+spin)+
#    scale_y_log10() + annotation_logticks() + ylab('(S2+1)')+xlab('copies per rxn')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()
#ggsave(paste0(rundir,'Saliva_S2.png'))
#ggplot(dfF, aes(x=virus_copy, y=S2_spike+1, group=virus_copy, color=Stotal>1000))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~treatment+spin)+
#
#    scale_y_log10() + annotation_logticks() + ylab('(S2_spike+1)')+xlab('copies per rxn')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()
#ggsave(paste0(rundir,'Saliva_S2_spike.png'))
#ggplot(dfF, aes(x=virus_copy, y=RPP30+1, group=virus_copy, color=Stotal>1000))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~treatment+spin)+
#
#    scale_y_log10() + annotation_logticks() + ylab('(RPP30+1)')+xlab('copies per rxn')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()
#ggsave(paste0(rundir,'Saliva_RPP30_spike.png'))
#
#ggplot(dfF[dfF$Stotal>1000,], aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~treatment+spin)+
#
#    scale_y_log10() + annotation_logticks() + ylab('(S2+1/S2_spike+1)')+xlab('copies per rxn')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()
#ggsave(paste0(rundir,'Saliva_S2nrom.png'))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#ggplot(dfF, aes(x=virus_copy, y=S2_spike+1, group=virus_copy, color=Plate_ID))+ #log10(S2_total_across_all_wells)))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~lysate)+ scale_y_log10() + annotation_logticks() + ylab('(S2_spike+1)')+xlab('copies pre rxn')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()+ggtitle('LoD')
#ggsave(paste0(rundir,'S2_spike.png'))
#
#ggplot(dfF, aes(x=virus_copy, y=RPP30+1, group=virus_copy, color=Plate_ID))+ #log10(S2_total_across_all_wells)))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~lysate)+ scale_y_log10() + annotation_logticks() + ylab('(RPP30+1)')+xlab('copies/mL')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()+ggtitle('LoD')
#ggsave(paste0(rundir,'RPP30_JB.png'))
#
#ggplot(dfF, aes(x=virus_copy, y=RPP30+1, group=virus_copy, color=Plate_ID))+ #log10(S2_total_across_all_wells)))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~lysate)+ scale_y_log10() + annotation_logticks() + ylab('(RPP30+1)')+xlab('copies/mL')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()+ggtitle('LoD')
#ggsave(paste0(rundir,'RPP30_spike_JB.png'))
#
#ggplot(dfF, aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy, color=Plate_ID))+ #log10(S2_total_across_all_wells)))+
#    geom_quasirandom(alpha=.75)+
#    facet_wrap(~lysate)+ scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2_spike+1)')+xlab('copies/mL')+
#    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
#    #scale_color_viridis(option = 'plasma')+
#    theme_bw()+ggtitle('LoD')+
#    geom_hline(yintercept=.003)
