library(ShortRead)
library(stringdist)
library(tidyverse)
library(ggbeeswarm)
library(viridis)
library(ggpubr)

swabseq.dir='/data/Covid/swabseq/'
outdir=paste0(swabseq.dir, 'analysis/manuscript/')

mungeTables=function(tables.RDS){
    countTables=readRDS(tables.RDS) #paste0(rundir, 'countTable.RDS')) 
    df=do.call('rbind', countTables)
    df$virus_copy=as.factor(df$virus_copy) 
    df$Col=as.factor(gsub('^.', '', df$Sample_Well))
    df$Row=factor(gsub('..$', '', df$Sample_Well), levels=rev(toupper(letters[1:8])))
    df$Sample=paste0(df$Plate_ID, '-' ,df$Sample_Well)
    df$Plate_ID=as.factor(df$Plate_ID)
    df$Plate_ID=factor(df$Plate_ID, levels(df$Plate_ID)[order(as.numeric(gsub('Plate', '', levels(df$Plate_ID))))])  
    if(!is.null(df$Plate_384)) {    df$Plate_384=as.factor(df$Plate_384) }
    if(length(table(df$amplicon))==3){
        df$amplicon=factor(df$amplicon, level=c('S2', 'S2_spike', 'RPP30'))
    } else {    df$amplicon=factor(df$amplicon, level=c('S2', 'S2_spike', 'RPP30', 'RPP30_spike')) }

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
    return(dfs)
}





#Figure 1 
rundir=paste0(swabseq.dir, 'runs/v11/')
dfs=mungeTables(paste0(rundir, 'countTable.RDS'))
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


# note can reuse fig1e2 as figure 2c
#LoD Purified 
rundir=paste0(swabseq.dir, 'runs/v19/')
dfsR=mungeTables(paste0(rundir, 'countTable.RDS'))
dfsR$Plate_384=droplevels(dfsR$Plate_384)
dfsR$virus_copy[dfsR$virus_copy=='Negative Patient']='0'
dfsR$virus_copy=as.factor(round(as.numeric(as.character(dfsR$virus_copy))*28/1000))
dfsR$virus_copy=droplevels(dfsR$virus_copy) #Plate_384=droplevels(dfsR$Plate_384)
dfsR$virus_copy=factor(dfsR$virus_copy, levels(dfsR$virus_copy)[order(as.numeric(levels(dfsR$virus_copy)))])
#dfsR$virus_copy=as.character(dfsR$virus_copy)
dfsR$virus_ident2=dfsR$virus_identity
dfsR$virus_ident2[grepl('^N|^P', dfsR$virus_identity)]='NegPatient'
dfsR$virus_ident2[grepl('^U', dfsR$virus_identity)]='PosPatient'

fig2a=
dfsR %>% filter(virus_ident2!='PosPatient') %>% filter(Stotal>2000) %>%
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+ #log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2_spike+1)')+xlab('GCE per reaction')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+ ggtitle('LoD - Purified RNA')+
    theme_bw()

#LoD Saliva 
rundir=paste0(swabseq.dir, 'runs/v13/')
dfs=mungeTables(paste0(rundir, 'countTable.RDS'))

dfsc = dfs %>% filter(lysate=='NS_Saliva_1_to_2') %>% filter(MasterMixVol=='20uL') %>% filter(Stotal>2000) 
#really, R?
dfsc$virus_copy=factor(as.character(round(as.numeric(as.character(dfsc$virus_copy)))), levels=c(sort(unique(as.numeric(as.character(dfsc$virus_copy))))))

fig2b=ggplot(dfsc, aes(x=virus_copy, y=S2_normalized_to_S2_spike))+geom_quasirandom(alpha=.75)+
     scale_y_log10() + annotation_logticks(sides='l') + ylab('(S2+1)/(S2_spike+1)')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+theme_bw()+
    xlab('GCE per reaction')+
    ggtitle('LoD - Saliva in Normal Saline 1:2')
fig2c=fig1e2+ggtitle('LoD - Mid Turbinate 1:2 TE')
ggarrange(fig2a,fig2b,fig2c, labels=c('A', 'B', 'C'),ncol=3)
ggsave('/data/Covid/swabseq/analysis/manuscript/fig2abc.png')
ggsave('/data/Covid/swabseq/analysis/manuscript/fig2abc.svg')



#quantile(dfs$S2_normalized_to_S2_spike[dfs$virus_ident2=='NegPatient'],.995)

#dfsRe=dfsR
#dfsRe=dfsRe[-(which(grepl('01|02|03|04|05|06', dfsRe$Col))), ]
#| dfsRe$virus_ident2=='TE'),]
#dfsRe$virus_copy[dfsRe$virus_ident2=='PosPatient']='Positive Patient'
#dfsRe$virus_copy=factor(dfsRe$virus_copy, unique(dfsRe$virus_copy)[c(2,8,7,6,5,4,3,1)]) #levels(dfsRe$virus_copy)[order(as.numeric(levels(dfsR$virus_copy)))])
#dfsRe %>% write_csv(paste0(rundir, 'PreliminaryLoD.csv'))


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





