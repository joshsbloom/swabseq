library(gdata)
library(data.table)
ed.table=read.xls('/data/Covid/swabseq/ED_Samples_Deidentified.xlsx', stringsAsFactors=F)

#v26
swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v26/')
outdir=paste0(swabseq.dir, 'analysis/v26/')
dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfs26=dfL$dfs


swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v27/')
outdir=paste0(swabseq.dir, 'analysis/v27/')

dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfs27=dfL$dfs

swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v28/')
outdir=paste0(swabseq.dir, 'analysis/v28/')

dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfs28=dfL$dfs

dfs=rbindlist(list(v26=dfs26, v27=dfs27, v28=dfs28), idcol='experiment')

edo=dfs[as.character(dfs$virus_identity) %in% ed.table$Matrix.Tube.Barcode | as.character(dfs$virus_identity) %in% ed.table$Deidenfied.ID,]
mu=match(edo$virus_identity, ed.table$Matrix.Tube.Barcode)
edo$virus_identity[!is.na(mu)]=na.omit(ed.table$Deidenfied.ID[mu])

edo=merge(edo, ed.table, all.x=T, by.x='virus_identity', by.y='Deidenfied.ID')

ct.char=as.character(edo$Ct)
ct.char[!is.na(ct.char)]=paste0("(",ct.char[!is.na(ct.char)])
ct.char[!is.na(ct.char)]=paste0(ct.char[!is.na(ct.char)] ,")")

ct.char[is.na(ct.char)]=''
ct.char2=paste(edo$virus_identity, ct.char)
edo$ct.char=ct.char2

edo %>%  filter(Stotal>1000) %>%  
    ggplot(aes(x=ct.char, y=S2_normalized_to_S2_spike,color=NPResult))+
    geom_point(alpha=.75, size=2)+
     #geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~experiment, scales='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('ED')#+
ggsave(paste0(outdir,'ED_Data_by_experiment.png'))

edo %>%  filter(Stotal>2000) %>%  filter(treatment!="TE+RNA Sec + QP") %>% filter(!(experiment=='v26' & Plate_ID!='Plate13')) %>%
    ggplot(aes(x=ct.char, y=S2_normalized_to_S2_spike, color=NPResult))+
    geom_point(alpha=.75, size=2)+
     #geom_quasirandom(alpha=.75, size=2)+
    scale_y_log10() + annotation_logticks(sides='l') + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('ED')#+
ggsave(paste0(outdir,'ED_Data_combined_v26_plate13_v27_v28.png'))


#for UCLA present
edof=edo %>%  filter(Stotal>2000) %>%  filter(treatment!="TE+RNA Sec + QP") %>% 
    filter(!(experiment=='v26' & Plate_ID!='Plate13')) %>% filter(!grepl('S-031', ct.char)) %>%
    filter(!grepl('S-005', ct.char))
    #filter(!grepl('S-022', ct.char))
#edof$NPResult=factor(edof$NPResult)   
#edodf$ct.char=factor(edodf$ct.char)
edof$NP_result=edof$NPResult=='DETECTED'

edof %>% 
    ggplot(aes(x=NP_result, y=S2_normalized_to_S2_spike, color=grepl('33', ct.char)))+
    #facet_wrap(~NP_result, scales='free_x')+
    scale_color_manual(values=c('black', 'orange'), name='Ct>32')+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    scale_y_log10() + annotation_logticks(sides='l') + xlab('SARS-CoV-2 detected by NP swab') +
    ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
     geom_hline(yintercept=5e-4, color='red')+
     ggtitle('ED Saliva Study')#+
ggsave(paste0(outdir,'ED_Data_plot_for_ucla_presentation.png'))

 
edo %>%  filter(Stotal>1000) %>%  
    ggplot(aes(x=virus_identity, y=S2_normalized_to_S2_spike, fill=))+
    geom_point(alpha=.75, size=2)+
     #geom_quasirandom(alpha=.75, size=2)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('ED_only')#+


dfs27 %>%  filter(Plate_ID=='Plate6' | Plate_ID=='Plate11') %>%  
    ggplot(aes(x=virus_identity, y=Stotal+1))+
    geom_point(alpha=.75, size=2)+

    facet_wrap(~lysate+treatment, scales='free_x')+
    scale_y_log10() + ylab('(S2 + S2 spike)')+xlab('identity')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


 dfs28  %>%   filter(Stotal>1000) %>% 
    ggplot(aes(x=virus_identity, y=S2_normalized_to_S2_spike))+
    geom_point(alpha=.75, size=2)+
     #geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~lysate+treatment, scales='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
