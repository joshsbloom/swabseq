library(ShortRead)
library(stringdist)
library(tidyverse)
library(ggbeeswarm)
library(viridis)
library(ggpubr)

swabseq.dir='/data/Covid/swabseq/'
outdir=paste0(swabseq.dir, 'analysis/manuscript/')

rundir=paste0(swabseq.dir, 'runs/v11/')

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
  #filter(Plate_ID!='Plate8') %>%
  mutate(SARS_COV_2_Detected=S2_normalized_to_S2_spike>.003)
dfs$SARS_COV_2_Detected[!dfs$RPP30_Detected]='Inconclusive'
dfs$SARS_COV_2_Detected[dfs$Stotal<2000]='Inconclusive' 
levels(dfs$virus_copy)=as.character(round(as.numeric(levels(dfs$virus_copy))))

fig1e1=
dfs%>%filter(lysate=='MTS_TE-1_to_2') %>% 
ggplot(aes(x=virus_copy, y=S2+1))+geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)')+xlab('GCE per reaction')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()
fig1e2=
dfs%>%filter(lysate=='MTS_TE-1_to_2') %>% 
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike))+geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 Spike + 1)')+xlab('GCE per reaction')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()

ggarrange(fig1e1,fig1e2)
ggsave('/data/Covid/swabseq/analysis/manuscript/fig1e.png')
ggsave('/data/Covid/swabseq/analysis/manuscript/fig1e.svg')



