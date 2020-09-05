swabseq.dir="/mnt/e/swabseq/"
swabseq.dir="/data/Covid/swabseq/"
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v45/')
outdir=paste0(swabseq.dir, 'analysis/v45/')

dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T, Stotal_filt=500, input=384)
df=dfL$df
dfs=dfL$dfs
#usable swabseq reads
#sum(df$Count)
#17477902

titl='v45 - rerun of v43 rapid run miniseq 1536 UDI plate 2 only'

#plate visualization 
#plate visualization 
df %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle(titl)
ggsave(paste0(outdir,'plateVis_all_indices.png'))

#filled.plates=add384Mapping(df)
#filled.plates %>% 
#    ggplot(aes(x=Col384, y=Row384, fill=log10(Count))) + 
#  geom_raster() +
#  coord_equal() +
#  facet_grid(amplicon~Plate_384) +
#  scale_fill_viridis_c(option = 'plasma')+ggtitle(titl)
#ggsave(paste0(outdir,'plateVis_384.png'))


df %>% filter(Plate_ID=='Plate2') %>% #filter(Plate_ID!='Plate1') %>%
    ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle(titl)
#ggtitle('v30 Saliva; Nasal; ED; Ashe')
ggsave(paste0(outdir,'plateVis_plates_run.png'))



