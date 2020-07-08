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

