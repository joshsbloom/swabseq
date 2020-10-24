swabseq.dir="/mnt/e/swabseq/"
swabseq.dir="/data/Covid/swabseq/"
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v56/')
outdir=paste0(swabseq.dir, 'analysis/v56/')

dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T, Stotal_filt=500, input=384)
df=dfL$df
dfs=dfL$dfs
#usable swabseq reads
#sum(df$Count)
#17477902
sum(dfs$S2)
sum(dfs$S2_spike)
sum(dfs$RPP30) 

dfs %>% filter(virus_identity=='MNS NS') %>% 
    # filter(virus_identity!='TE') %>%
    #filter(Stotal>500) %>%
ggplot(aes(x=virus_copy, y=S2_normalized_to_S2_spike,color=Stotal>500))+
geom_quasirandom(size=2, alpha=.75)+
    scale_y_log10() + annotation_logticks(sides="l")+
    #geom_hline(yintercept=3e-3, color='red')+
    theme_bw()+geom_hline(yintercept=3e-3)+
    ylab('(S2 + 1)/(S2_spike + 1)')+
    #facet_grid(~NPResult, scales='free_x')+
    theme(axis.text.x = element_text(angle = 90))+
    ggtitle('MNS NS 1:4 - gamma')

