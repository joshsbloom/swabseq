swabseq.dir="/mnt/e/swabseq/"
swabseq.dir="/data/Covid/swabseq/"
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v46/')
outdir=paste0(swabseq.dir, 'analysis/v46/')

dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T, Stotal_filt=2000, input=384)
df=dfL$df
dfs=dfL$dfs
#usable swabseq reads
#sum(df$Count)
#17477902

sum(dfs$S2)+sum(dfs$S2_spike)
sum(dfs$RPP30)

titl='v46 - NP Purified 40 cycles'

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



#Combine with higher CT samples 

out$virus_ident2=''
out$virus_ident2[out$NPResult==TRUE]='Positive'
out$virus_ident2[out$NPResult==FALSE]='Negative'





rundir=paste0(swabseq.dir, 'runs/v19/')
dfsR=mungeTables(paste0(rundir, 'countTable.RDS'))
dfsR$Plate_384=droplevels(dfsR$Plate_384)
dfsR$virus_copy[dfsR$virus_copy=='Negative Patient']='0'
dfsR$virus_copy=droplevels(dfsR$virus_copy) #Plate_384=droplevels(dfsR$Plate_384)
dfsR$virus_copy=factor(dfsR$virus_copy, levels(dfsR$virus_copy)[order(as.numeric(levels(dfsR$virus_copy)))])
#dfsR$virus_copy=as.character(dfsR$virus_copy)
dfsR$virus_ident2=dfsR$virus_identity
dfsR$virus_ident2[grepl('^N|^P', dfsR$virus_identity)]='Negative'
dfsR$virus_ident2[grepl('^U', dfsR$virus_identity)]='Positive'
dfsR=dfsR[dfsR$SARS_COV_2_Detected!='Inconclusive',]

highCT=dfsR %>% filter(virus_ident2=='Negative' | virus_ident2=='Positive') %>% filter(SARS_COV_2_Detected!='Inconclusive')


comb=data.frame(virus_ident2=c(highCT$virus_ident2, out$virus_ident2),
            S2_normalized_to_S2_spike=c(highCT$S2_normalized_to_S2_spike, out$S2_normalized_to_S2_spike))

fig2e=dfsR %>% filter(virus_ident2=='Negative' | virus_ident2=='Positive') %>% filter(SARS_COV_2_Detected!='Inconclusive') %>% 
ggplot(comb,aes(x=virus_ident2, y=S2_normalized_to_S2_spike))+ #, color=SARS_COV_2_Detected))+ #spike+S2, group=virus_copy))+ #log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks(sides='l') + ylab('(S2+1)/(S2 spike+1)')+ #xlab('patient SARs-CoV-2 status')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+ xlab("")+
    geom_hline(yintercept=.003, color='red') +

    theme_bw()+ggtitle('Nasal Swab, Purified RNA') 
ggsave('/data/Covid/swabseq/analysis/v46/purified_combined.png')

    #facet_wrap(~virus_ident2)+ 
    #scale_color_viridis(option = 'plasma')+
    +    
