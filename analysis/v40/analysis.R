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
