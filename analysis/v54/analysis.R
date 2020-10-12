swabseq.dir="/mnt/e/swabseq/"
swabseq.dir="/data/Covid/swabseq/"
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v54/')
outdir=paste0(swabseq.dir, 'analysis/v54/')

dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T, Stotal_filt=500, input=384)
df=dfL$df
dfs=dfL$dfs
#usable swabseq reads
#sum(df$Count)
#17477902
sum(dfs$S2)
sum(dfs$S2_spike)
sum(dfs$RPP30)

#dfa=dfs %>% filter(Plate_ID=='Plate1') %>% filter(quadrant_96=='D') 
#p1.s2=sum(dfa$S2)
#p1.s2_spike=sum(dfa$S2_spike)
#p1.rpp30=sum(dfa$RPP30)
#
#
#dfa=dfs %>% filter(Plate_ID=='Plate2') %>% filter(quadrant_96=='B') 
#
#p2.s2=sum(dfa$S2)
#p2.s2_spike=sum(dfa$S2_spike)
#p2.rpp30=sum(dfa$RPP30)
#par(mfrow=c(2,1))
#barplot(c(p1.s2,p1.s2_spike, p1.rpp30)/sum(c(p1.s2,p1.s2_spike, p1.rpp30)))
#barplot(c(p2.s2,p2.s2_spike, p2.rpp30)/sum(c(p2.s2,p2.s2_spike, p2.rpp30)))


titl='v54 - Gamma gradient taqpath vs Luna /  MNS NS'

#plate visualization 
#plate visualization 
df %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID+Description) +
  scale_fill_viridis_c(option = 'plasma')+ggtitle(titl)
ggsave(paste0(outdir,'plateVis_all_indices.png'))


#LoD Luna
table((dfs %>% filter(quadrant_96=='A') %>% 
    filter(Plate_ID=='Plate2') %>% filter(virus_identity!='0') %>%
     filter(virus_identity!='TE'))$SARS_COV_2_Detected)

a=dfs %>% filter(quadrant_96=='A') %>% 
    filter(Plate_ID=='Plate2') %>% filter(virus_identity!='0') %>%
    # filter(virus_identity!='TE') %>%
    filter(Stotal>500) %>%
ggplot(aes(x=log10(as.numeric(as.character(virus_copy))
                  +1), y=S2_normalized_to_S2_spike))+
geom_quasirandom(size=2, alpha=.75)+
    scale_y_log10(limits=(c(1e-5,100))) + annotation_logticks(sides="l")+
    #geom_hline(yintercept=3e-3, color='red')+
    theme_bw()+geom_hline(yintercept=3e-3)+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    #facet_grid(~NPResult, scales='free_x')+
    theme(axis.text.x = element_text(angle = 90))+
    ggtitle('Gamma Saliva LoD _ Luna + UDG index plate 2')

b=dfs %>% filter(quadrant_96=='A') %>% 
    filter(Plate_ID=='Plate4') %>% filter(virus_identity!='0') %>%
    # filter(virus_identity!='TE') %>%
    filter(Stotal>500) %>%
ggplot(aes(x=log10(as.numeric(as.character(virus_copy))
                  +1), y=S2_normalized_to_S2_spike))+
geom_quasirandom(size=2, alpha=.75)+
    scale_y_log10(limits=(c(1e-5,100))) + annotation_logticks(sides="l")+
    #geom_hline(yintercept=3e-3, color='red')+
    theme_bw()+geom_hline(yintercept=3e-3)+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    #facet_grid(~NPResult, scales='free_x')+
    theme(axis.text.x = element_text(angle = 90))+
    ggtitle('Gamma Saliva LoD _ Taqpath index plate 4')

library(ggpubr)
ggarrange(a,b)
ggsave(paste0(outdir,'Gamma_saliva_LoD_NEB_vs_Taqpath.png'))



pSampAsheLuna=
    dfs %>% filter(quadrant_96=='B') %>% 
    filter(Plate_ID=='Plate2') %>% filter(virus_identity!='0') %>%
    # filter(virus_identity!='TE') %>%
    filter(Stotal>500) %>%
ggplot(aes(x=virus_identity, y=S2_normalized_to_S2_spike))+
geom_quasirandom(size=2, alpha=.75)+
    scale_y_log10(limits=(c(1e-5,100))) + annotation_logticks(sides="l")+
    #geom_hline(yintercept=3e-3, color='red')+
    theme_bw()+geom_hline(yintercept=3e-3)+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    #facet_grid(~NPResult, scales='free_x')+
    theme(axis.text.x = element_text(angle = 90))+
    ggtitle('Ashe + ED _  Luna + UDG index plate 2')


pSampAsheTP=
    dfs %>% filter(quadrant_96=='B') %>% 
    filter(Plate_ID=='Plate4') %>% filter(virus_identity!='0') %>%
    # filter(virus_identity!='TE') %>%
    filter(Stotal>500) %>%
ggplot(aes(x=virus_identity, y=S2_normalized_to_S2_spike))+
geom_quasirandom(size=2, alpha=.75)+
    scale_y_log10(limits=(c(1e-5,100))) + annotation_logticks(sides="l")+
    #geom_hline(yintercept=3e-3, color='red')+
    theme_bw()+geom_hline(yintercept=3e-3)+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    #facet_grid(~NPResult, scales='free_x')+
    theme(axis.text.x = element_text(angle = 90))+
    ggtitle('Ashe + ED _  Taqpath + UDG index plate 2')

ggarrange(pSampAsheLuna,pSampAsheTP, nrow=2)
ggsave(paste0(outdir,'saliva_taqpath_v_luna_psamp.png'))





library(gdata)
library(ggrepel)
pd=read.xls(paste0(swabseq.dir,'analysis/v37/', 'EUA_data.xlsx'), stringsAsFactors=F)
pd$Sample[2]='R-003'
pd$Sample[19]='R-006'
#dfsF =dfs %>%  filter(Plate_ID=='Plate15')
#dfsF=merge(dfsF, pd, by.x='virus_identity', by.y='Sample', all.x=T)
#for genapsys
dfsF=merge(dfs, pd, by.x='virus_identity', by.y='Sample', all.x=T)
names(dfsF)[38]='Ct'
names(dfsF)[41]='NPResult'

dfsFA=dfsF %>% filter(NPResult!='NA') %>% filter(quadrant_96 !='A') %>% filter(quadrant_96 !='B')
gt=ggplot(dfsFA, aes(x=Ct, y=S2_normalized_to_S2_spike, color=SARS_COV_2_Detected, label=virus_identity))+
geom_quasirandom(size=2)+
    scale_y_log10() + annotation_logticks(sides="l")+geom_label_repel()+
    geom_hline(yintercept=3e-3, color='red')+theme_bw()+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    facet_grid(~NPResult+quadrant_96, scales='free_x')+ggtitle('MNS NS taqpath')

gt2=ggplot_gtable(ggplot_build(gt))
gt2$widths[5]=.4*gt2$widths[5]
gt2$widths[7]=.4*gt2$widths[7]
grid.draw(gt2)
ggsave(paste0(outdir,'MNS_NS_taqpath.png'))







#ggsave(paste0(outdir,'Gamma_saliva_LoD.png'))


#one more inconclusive for taqpath
table((dfs %>% filter(quadrant_96=='C' | quadrant_96=='D') %>% 
    filter(Plate_ID=='Plate1') %>% filter(virus_identity!='IGNORE') %>%
     filter(virus_identity!='TE') )$SARS_COV_2_Detected)

b=dfs %>% filter(quadrant_96=='C' | quadrant_96=='D') %>% 
    filter(Plate_ID=='Plate1') %>% filter(virus_identity!='IGNORE') %>%
     filter(virus_identity!='TE') %>%

    filter(Stotal>500) %>%
ggplot(aes(x=log10(as.numeric(as.character(virus_copy))
                  +1), y=S2_normalized_to_S2_spike))+
geom_quasirandom(size=2, alpha=.75)+
    scale_y_log10(limits=(c(1e-6,100))) + annotation_logticks(sides="l")+
    #geom_hline(yintercept=3e-3, color='red')+
    theme_bw()+geom_hline(yintercept=3e-3)+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    #facet_grid(~NPResult, scales='free_x')+
    theme(axis.text.x = element_text(angle = 90))+
    ggtitle('Gamma Saliva LoD _ Taqpath _ index plate 1')
#ggsave(paste0(outdir,'Gamma_saliva_LoD.png'))



#lusis experiment
table_out=dfs %>% filter(quadrant_96=='A' | quadrant_96=='B') %>% 
    filter(Plate_ID=='Plate1')%>% filter(virus_identity!='IGNORE') # %>%
#    filter(virus_identity=!='TE')
table_out$SARS_COV_2_Detected=table_out$S2_normalized_to_S2_spike>.003

table_out %>% write.csv(paste0(outdir,'mouseResults.csv'))


#empirical estimate of s2 spike
dfa=dfs %>% filter(quadrant_96=='A' | quadrant_96=='B') %>% 
    filter(Plate_ID=='Plate3') %>% filter(virus_identity!='IGNORE') %>%
     filter(virus_identity!='TE') %>%
    filter(Stotal>500)
plot(log10(as.numeric(as.character(dfa$virus_copy))/285+1),
           log10(dfa$S2_normalized_to_S2_spike))
#abline(lm(log10(dfa$S2_normalized_to_S2_spike)~log10(as.numeric(as.character(dfa$virus_copy))/285+1)))
dfa2=dfa[dfa$virus_copy!='0',]

abline(lm(log10(dfa2$S2_normalized_to_S2_spike)~log10(as.numeric(as.character(dfa2$virus_copy))/285+1)))
