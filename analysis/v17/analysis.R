library(ShortRead)
library(stringdist)
library(tidyverse)
library(ggbeeswarm)
library(viridis)

#for example ... 
rundir='/data/Covid/swabseq/runs/v17/'
outdir='/data/Covid/swabseq/analysis/v17/'

countTables=readRDS(paste0(rundir, 'countTable.RDS')) 
df=do.call('rbind', countTables)
df$virus_copy=as.factor(df$virus_copy) 
df$Col=as.factor(gsub('^.', '', df$Sample_Well))
df$Row=factor(gsub('..$', '', df$Sample_Well), levels=rev(toupper(letters[1:8])))
df$Sample=paste0(df$Plate_ID, '-' ,df$Sample_Well)
df$Plate_ID=as.factor(df$Plate_ID)
df$Plate_ID=factor(df$Plate_ID, levels(df$Plate_ID)[order(as.numeric(gsub('Plate', '', levels(df$Plate_ID))))])  
df$Plate_384=as.factor(df$Plate_384)
df$amplicon=factor(df$amplicon, level=c('S2', 'S2_spike', 'RPP30', 'RPP30_spike'))

#plate visualization 
df %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID) +
  scale_fill_viridis_c(option = 'plasma')
ggsave(paste0(outdir,'plateVis.png'))


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
  filter(Plate_ID!='Plate8') %>%
  mutate(SARS_COV_2_Detected=S2_normalized_to_S2_spike>.003)
dfs$SARS_COV_2_Detected[!dfs$RPP30_Detected]='Inconclusive'
dfs$SARS_COV_2_Detected[dfs$Stotal<2000]='Inconclusive'
dfs$Plate_384=droplevels(dfs$Plate_384)
dfs$virus_copy[dfs$virus_copy=='Negative Patient']='0.0'
dfs$virus_copy=as.numeric(as.character(dfs$virus_copy))*1000
dfs$virus_copy=droplevels(as.factor(dfs$virus_copy))
dfs$virus_copy=factor(dfs$virus_copy, levels(dfs$virus_copy)[order(as.numeric(levels(dfs$virus_copy)))])
dfs$virus_ident2=dfs$virus_identity
dfs$virus_ident2[grepl('^N|^P', dfs$virus_identity)]='NegPatient'
dfs$virus_copy=(as.numeric(as.character(dfs$virus_copy)))
dfs$zeros=factor(dfs$virus_copy==0, levels=c(TRUE, FALSE))

dfs$virus_copy[dfs$virus_copy==0]=1
dfs %>% filter(virus_ident2!='D1' & virus_ident2!='D2') %>% filter(Stotal>2000) %>% 
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike))+
    geom_quasirandom(alpha=.75)+
    facet_grid(~zeros, scales='free_x', space='free')+
    scale_y_log10() + annotation_logticks() + ylab('(S2 + 1)/(S2 spike + 1)')+xlab('log10(copies/mL )')+
     scale_x_log10()+
    #scale_x_continuous(trans=scales::pseudo_log_trans(base=10))+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+ggtitle('LoD')

dfs %>% filter(virus_ident2!='D1' & virus_ident2!='D2') %>% filter(Stotal>2000) %>% 
ggplot(aes(x=virus_copy, y=S2+1))+
    geom_quasirandom(alpha=.75)+
    facet_grid(~zeros, scales='free_x', space='free')+
    scale_y_log10(limits=c(1,1e5)) + annotation_logticks() + ylab('(S2 + 1)')+xlab('log10(copies/mL )')+
     scale_x_log10()+
     #scale_x_continuous(trans=scales::pseudo_log_trans(base=10))+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+ggtitle('LoD')

dfs %>% filter(virus_ident2!='D1' & virus_ident2!='D2') %>% filter(Stotal>2000) %>% 
ggplot(aes(x=virus_copy, y=S2_spike+1))+
    geom_quasirandom(alpha=.75)+
    facet_grid(~zeros, scales='free_x', space='free')+
    scale_y_log10(limits=c(1,1e5)) + annotation_logticks() + ylab('(S2 spike + 1)')+xlab('log10(copies/mL )')+
     scale_x_log10()+
     #scale_x_continuous(trans=scales::pseudo_log_trans(base=10))+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+ggtitle('LoD')






#dilution subset
dfs= df %>%filter(amplicon=='S2') %>%  
  #count(Sample_Well, wt=Count, name='S2_total_across_all_wells') %>%
  right_join(df)
dfs= dfs %>% count(Sample, wt=Count, name='well_total') %>%
  right_join(dfs) %>% 
  filter(Plate_ID!='Plate8') %>%
  filter(well_total>1000) %>%
  #filter(Plate_384=='3') %>%
  select(-mergedIndex, -Sample_ID, -index, -index2 ) %>% 
  spread(amplicon, Count) %>% 
  mutate(S2_normalized_to_S2_spike=(S2+1)/(S2_spike+1))
dfs$Plate_384=droplevels(dfs$Plate_384)
dfs$virus_copy[dfs$virus_copy=='Negative Patient']='0.0'
dfs$virus_copy=as.numeric(as.character(dfs$virus_copy))*1000
dfs$virus_copy=droplevels(as.factor(dfs$virus_copy))
dfs$virus_copy=factor(dfs$virus_copy, levels(dfs$virus_copy)[order(as.numeric(levels(dfs$virus_copy)))])
dfs$virus_ident2=dfs$virus_identity
dfs$virus_ident2[grepl('^N|^P', dfs$virus_identity)]='NegPatient'
dfs$virus_copy=as.numeric(as.character(dfs$virus_copy))




ggplot(dfs, aes(x=virus_copy+1, y=S2_normalized_to_S2_spike, group=Plate_ID, color=Plate_ID))+geom_quasirandom(alpha=.75)+
    facet_wrap(~virus_ident2+lysate) +scale_x_log10()+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2_spike+1)')+xlab('log10(copies/mL +1)')+

    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+ggtitle('LoD')














ggplot(dfs, aes(x=virus_copy, y=S2+1, group=virus_copy)) + #, color=log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75)+
    facet_wrap(~virus_ident2+lysate+Plate_ID)+  scale_y_log10() + annotation_logticks() + ylab('(S2+1)')+xlab('copies/mL')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
    scale_color_viridis(option = 'plasma')+
    theme_bw()+ggtitle('LoD')
ggsave(paste0(rundir,'S2.png'))

ggplot(dfs, aes(x=virus_copy, y=S2_spike+1, group=virus_copy, color=Plate_ID))+geom_quasirandom(alpha=.75)+
    facet_wrap(~virus_ident2+lysate+Plate_ID) +
     scale_y_log10() + annotation_logticks() + ylab('(S2_spike+1)')+xlab('copies/mL')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+ggtitle('LoD')
ggsave(paste0(rundir,'S2_spike.png'))


ggplot(dfs, aes(x=virus_copy, y=RPP30+1, group=virus_copy, color=Plate_ID))+geom_quasirandom(alpha=.75)+
    facet_wrap(~virus_ident2+lysate+Plate_ID) +
     scale_y_log10() + annotation_logticks() + ylab('RPP30')+xlab('copies/mL')+

    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+ggtitle('LoD')
ggsave(paste0(rundir,'RPP30.png'))



ggplot(dfs, aes(x=virus_copy+1, y=S2_normalized_to_S2_spike, group=Plate_ID, color=Plate_ID))+geom_quasirandom(alpha=.75)+
    facet_wrap(~virus_ident2+lysate) +scale_x_log10()+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2_spike+1)')+xlab('log10(copies/mL +1)')+

    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+ggtitle('LoD')
ggsave(paste0(rundir,'S2_to_S2_spike.png'))
