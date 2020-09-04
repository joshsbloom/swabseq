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

titl='v41 - Ashe Negative  + NS 1:4 Clinical Samples'

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


library(gdata)
pd=read.xls(paste0(swabseq.dir,'analysis/v37/', 'EUA_data.xlsx'), stringsAsFactors=F)
pd$Sample[2]='R-003'
pd$Sample[19]='R-006'
#dfsF =dfs %>%  filter(Plate_ID=='Plate15')
#dfsF=merge(dfsF, pd, by.x='virus_identity', by.y='Sample', all.x=T)
#for genapsys
dfsF=merge(dfs, pd, by.x='virus_identity', by.y='Sample', all.x=T)
#names(dfsF)[41]='CovidDetectedNP'
names(dfsF)[38]='Ct'
names(dfsF)[41]='NPResult'

dfsF=dfsF %>% filter(Plate_ID=='Plate4' | Plate_ID=='Plate5'| Plate_ID=='Plate10' | Plate_ID=='Plate15')
dfsF$Ct[dfsF$Ct==0]=NA
dfsF %>% filter(Plate_ID=='Plate15') %>% write.csv(paste0(outdir, 'resultsTable.csv'))
library(grid)

gt=ggplot(dfsF, aes(x=Ct, y=S2_normalized_to_S2_spike))+
geom_quasirandom(size=2)+
    scale_y_log10() + annotation_logticks(sides="l")+
    geom_hline(yintercept=3e-3, color='red')+theme_bw()+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    facet_grid(~NPResult, scales='free_x')

gt2=ggplot_gtable(ggplot_build(gt))
gt2$widths[5]=.3*gt2$widths[5]
grid.draw(gt2)
#png(paste0(outdir,'NS_patients.png'), width=512, height=400)
#grid.draw(gt2)
#dev.off()
#ggsave(paste0(outdir,'NS_patients.png'))
#ggarrange(gt2)













library(ggrepel)
gt=dfsF %>%  
gt=    ggplot(dfsF, aes(x=virus_identity, y=S2_normalized_to_S2_spike, color=Stotal>500, label=S_gene))+
    geom_text_repel()+
    geom_point(alpha=.75, size=2)+
    #geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~CovidDetectedNP, scales='free_x')+
    scale_y_log10() + 
    annotation_logticks() + ylab('(S2+1)/(S2 spike +1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('NS MNS 1:4 LoD')



