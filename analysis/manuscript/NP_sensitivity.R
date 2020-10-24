#v19 prelim LoD
#500, 1000,2000,and 8000
rundir=paste0(swabseq.dir, 'runs/v19/')
dfsR=mungeTables(paste0(rundir, 'countTable.RDS'))
dfsR$Plate_384=droplevels(dfsR$Plate_384)
dfsR$virus_copy[dfsR$virus_copy=='Negative Patient']='0'
dfsR$virus_copy=as.factor(round(as.numeric(as.character(dfsR$virus_copy)))) #*28/1000))
dfsR$virus_copy=droplevels(dfsR$virus_copy) #Plate_384=droplevels(dfsR$Plate_384)
dfsR$virus_copy=factor(dfsR$virus_copy, levels(dfsR$virus_copy)[order(as.numeric(levels(dfsR$virus_copy)))])
#dfsR$virus_copy=as.character(dfsR$virus_copy)
dfsR$virus_ident2=dfsR$virus_identity
dfsR$virus_ident2[grepl('^N|^P', dfsR$virus_identity)]='NegPatient'
dfsR$virus_ident2[grepl('^U', dfsR$virus_identity)]='PosPatient'
dfsR=dfsR %>% filter(virus_ident2!='PosPatient') %>% filter(Stotal>2000)

ggplot(dfsR, aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+ #log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike+1)')+xlab('GCE per mL')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+ ggtitle('Nasal Swab, Purified RNA')+
    theme_bw()


A=dfsR %>% filter(dfsR$virus_copy=='500' | dfsR$virus_copy=='2000')


#v20 confirmatory LoD, just 500

swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v20/')
outdir=paste0(swabseq.dir, 'analysis/EUA/')
dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfsC=dfL$dfs

#filters and reformat for plots 
swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v20/')
outdir=paste0(swabseq.dir, 'analysis/EUA/')
dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfsC=dfL$dfs
dfsC= dfsC %>% 
  filter(Plate_ID=='Plate8' ) %>%filter(Stotal>2000)
dfsC$virus_copy=droplevels(dfsC$virus_copy) #Plate_384=droplevels(dfsC$Plate_384)
dfsC$virus_copy=factor(dfsC$virus_copy, levels(dfsC$virus_copy)[order(as.numeric(levels(dfsC$virus_copy)))])
#dfsC=dfsC[dfsC$virus_copy!='125',]

ggplot(dfsC, aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+ #log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2_spike+1)')+xlab('copies/mL')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
    theme_bw()+ggtitle('LoD Confirmation')+
    geom_hline(yintercept=.003, color='red')
B=dfsC %>% filter(dfsC$virus_copy=='500' | dfsC$virus_copy=='2000')


#v21 500 and 2000
swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v21/')
outdir=paste0(swabseq.dir, 'analysis/EUA/')
dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T)
df=dfL$df
dfsC=dfL$dfs
dfsC= dfsC %>% 
    filter(Plate_ID=='Plate8' | Plate_ID=='Plate11') %>%filter(Stotal>2000)
dfsC$virus_copy=droplevels(dfsC$virus_copy) #Plate_384=droplevels(dfsC$Plate_384)
dfsC$virus_copy=factor(dfsC$virus_copy, levels(dfsC$virus_copy)[order(as.numeric(levels(dfsC$virus_copy)))])
#dfsC=dfsC[dfsC$virus_copy!='125',]
ggplot(dfsC, aes(x=virus_copy, y=S2_normalized_to_S2_spike, group=virus_copy))+ #log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2_spike+1)')+xlab('copies/mL')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
    theme_bw()+ggtitle('LoD Confirmation')+
    geom_hline(yintercept=.003, color='red')
CC=dfsC %>% filter(dfsC$virus_copy=='500' | dfsC$virus_copy=='2000')

three.expt=data.frame(
           experiment=c(rep('19', nrow(A)), rep('20', nrow(B)), rep('21', nrow(CC))),
           virus_copy=as.character(c(as.character(A$virus_copy), as.character(B$virus_copy), as.character(CC$virus_copy))),
           ratio=c(A$S2_normalized_to_S2_spike,B$S2_normalized_to_S2_spike,CC$S2_normalized_to_S2_spike))


three.expt %>% write.csv("/data/Covid/swabseq/analysis/manuscript/NP_variance.txt")

ggplot(three.expt, aes(x=virus_copy, y=ratio))+ #log10(S2_total_across_all_wells)))+
    geom_quasirandom(alpha=.75)+facet_grid(~experiment)+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2_spike+1)')+xlab('copies/mL')+
    theme(axis.text.x = element_text(angle = 90, vjust=0.3))+
    theme_bw()+ggtitle('LoD Confirmation')+
    geom_hline(yintercept=.003, color='red')

