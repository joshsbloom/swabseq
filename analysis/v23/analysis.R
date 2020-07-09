 #v23
swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
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

dfs %>%  filter(Description!='' & Description!=' ') %>%
  ggplot(aes(x=Col, y=Row,fill=SARS_COV_2_Detected)) + 
  geom_raster() +
  coord_equal() +
  scale_fill_manual(values=c("lightblue", "white", "red"))+
  facet_grid(~Plate_384+Plate_ID+Description) 
ggsave(paste0(outdir,'plateVis_CovidDetected.png'))
 
dfs %>%  filter(Description!='' & Description!=' ') %>%
ggplot(aes(x=Col, y=Row,fill=Stotal>2000)) + 
  geom_raster() +
  coord_equal() +
  #scale_fill_manual(values=c("lightblue", "white", "red"))+
  facet_grid(~Plate_384+Plate_ID+Description) 
ggsave(paste0(outdir,'plateVis_Stotal_gt_2000.png'))


dfs %>%  filter(Description!='' & Description!=' ') %>%
ggplot(aes(x=Col, y=Row,fill=Stotal>1000)) + 
  geom_raster() +
  coord_equal() +
  #scale_fill_manual(values=c("lightblue", "white", "red"))+
  facet_grid(~Plate_384+Plate_ID+Description) 
ggsave(paste0(outdir,'plateVis_Stotal_gt_1000.png'))


dfs %>%  filter(grepl('VTM', Description)) %>% filter(Stotal>2000) %>%
ggplot(aes(x=virus_copy, y=S2+1, group=virus_copy, color=lysate))+
    geom_quasirandom(alpha=.75)+
    facet_wrap(~lysate)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)')+xlab('copies per mL of lysate')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
    #scale_color_viridis(option = 'plasma')+
    theme_bw()

dfs %>%  filter(grepl('VTM', Description)) %>%  filter(Stotal>1000) %>%
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy, color=lysate))+
    geom_quasirandom(alpha=.75)+
    facet_wrap(~lysate)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/S2 spike + 1)')+xlab('copies per mL of lysate')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
    #scale_color_viridis(option = 'plasma')+
    theme_bw()
ggsave(paste0(outdir,'XY_VTM_stotal_gt_1000.png'))

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

