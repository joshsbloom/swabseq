swabseq.dir="/mnt/e/swabseq/"
swabseq.dir="/data/Covid/swabseq/"
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v49/')
outdir=paste0(swabseq.dir, 'analysis/v49/')

dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T, Stotal_filt=500, input=384)
df=dfL$df
dfs=dfL$dfs
#usable swabseq reads
#sum(df$Count)
#17477902

titl='v49 - Ashe / Gamma LoD Saliva'

#plate visualization 
#plate visualization 
df %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle(titl)
ggsave(paste0(outdir,'plateVis_all_indices.png'))


df %>% filter(Description!='' & Description!=' ') %>% #filter(Plate_ID!='Plate1') %>%
    ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle(titl)
#ggtitle('v30 Saliva; Nasal; ED; Ashe')
ggsave(paste0(outdir,'plateVis_plates_run.png'))

# Ashe samples 
dfs %>% filter(grepl('^\\d\\d', virus_identity)) %>% 
ggplot(aes(x=virus_identity, y=S2_normalized_to_S2_spike, color=SARS_COV_2_Detected,
           fill =RPP30_Detected))+
geom_point(size=2, alpha=.75)+
    scale_y_log10() + annotation_logticks(sides="l")+
    #geom_hline(yintercept=3e-3, color='red')+
    theme_bw()+geom_hline(yintercept=3e-3)+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    #facet_grid(~NPResult, scales='free_x')+
    theme(axis.text.x = element_text(angle = 90))+
    ggtitle('Ashe Samples')
ggsave(paste0(outdir,'Ashe.png'))

#LoD
dfs %>% filter(lysate=='C' | lysate=='D') %>% filter(Stotal>500) %>%
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike))+
geom_quasirandom(size=2, alpha=.75)+
    scale_y_log10() + annotation_logticks(sides="l")+
    #geom_hline(yintercept=3e-3, color='red')+
    theme_bw()+geom_hline(yintercept=3e-3)+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    #facet_grid(~NPResult, scales='free_x')+
    theme(axis.text.x = element_text(angle = 90))+
    ggtitle('Gamma Saliva LoD')
ggsave(paste0(outdir,'Gamma_saliva_LoD.png'))


df %>% filter(Description!='' & Description!=' ') %>% filter(lysate=='C') %>% #filter(Plate_ID!='Plate1') %>%
    ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle(titl)
#ggtitle('v30 Saliva; Nasal; ED; Ashe')


#ggsave(paste0(outdir,'plateVis_plates_run.png'))
pid=c('U166',
'U266',	
'U201',	
'U258',	
'U215',	
'U191',	
'U241',	
'U234',	
'U281',	
'U252',	
'U267',	
'U238',	
'U225',	
'U185',	
'U176',	
'U279',	
'U228',	
'U204',	
'U247',	
'U253')	


pct=c(30.1545,
30.5463,
30.8456,
30.3,
30.5748,
31.0232,
31.1,
30.3541,
30.4169,
31.0506,
33.0172,
33.1649,
33.4889,
33.0271,
33.1793,
33.8233,
33.1,
33.1893,
34.9525, 35)


library(xlsx)
dfct=read.xlsx('/data/Covid/swabseq/analysis/v46/NPCts.xlsx', 1)
names(dfct)[1]='virus_identity'

#dfct=data.frame(virus_identity=pid, Ct=pct)
dfsM=merge(dfct, dfs, by='virus_identity', all.y=T)
dfsM=dfsM %>%filter(Plate_ID=='Plate3')
dfsM=dfsM[grep('^U|^N' , dfsM$virus_identity),]

dfsM$NPResult=grepl('^U', dfsM$virus_identity)
dfsM$S.gene.Ct..KF.ABI.7500.[is.na(dfsM$S.gene.Ct..KF.ABI.7500.)]=40
    #names(dfsF)[38]='Ct'
#names(dfsF)[41]='NPResult'
dfsM= dfsM %>%filter(Stotal>2000)
dfsM$PCR_cycles=40

gt=ggplot(dfsM, aes(x=S.gene.Ct..KF.ABI.7500., y=S2_normalized_to_S2_spike))+
geom_quasirandom(size=2, alpha=.75)+
    scale_y_log10() + annotation_logticks(sides="l")+
    #geom_hline(yintercept=3e-3, color='red')+
    theme_bw()+geom_hline(yintercept=3e-3)+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    facet_grid(~NPResult, scales='free_x')+
    ggtitle('Purified NP, high Ct')
    

gt2=ggplot_gtable(ggplot_build(gt))
gt2$widths[5]=.3*gt2$widths[5]
library(ggpubr)
ggarrange(gt2)

out=dfsM[match(unique(as.character(dfsM$virus_identity)), dfsM$virus_identity),]
out[8,]= dfsM[15,]
out[13,]=dfsM[25,]
out %>% write.csv(paste0(outdir, 'NP_highCT.csv'))

gt=ggplot(out, aes(x=S.gene.Ct..KF.ABI.7500., y=S2_normalized_to_S2_spike))+
geom_quasirandom(size=2, alpha=.75)+
    scale_y_log10() + annotation_logticks(sides="l")+
    #geom_hline(yintercept=3e-3, color='red')+
    theme_bw()+geom_hline(yintercept=3e-3)+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    facet_grid(~NPResult, scales='free_x')+
    ggtitle('Purified NP, high Ct')
    

gt2=ggplot_gtable(ggplot_build(gt))
gt2$widths[5]=.3*gt2$widths[5]
library(ggpubr)
ggarrange(gt2)
ggsave(paste0(outdir,'purifiedNP_v_Ct.png'))
