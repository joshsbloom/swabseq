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
edo%>%filter(Stotal>500)%>%
ggplot(aes(x=Ct, y=S2_normalized_to_S2_spike))+geom_quasirandom()+
    scale_y_log10() + annotation_logticks(sides="l")+geom_hline(yintercept=3e-3, color='red')+theme_bw()+
    facet_wrap(~NPResult, scales='free_x')
