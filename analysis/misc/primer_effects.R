library(gdata)
library(data.table)
#ed.table=read.xls('/data/Covid/swabseq/EDSamples_073120.xlsx', stringsAsFactors=F) #Deidentified.xlsx', stringsAsFactors=F)
#ed.table=read.delim('/data/Covid/swabseq/EDSamples_073120.csv', stringsAsFactors=F) #Deidentified.xlsx', stringsAsFactors=F)
#ed.table=read.delim('/data/Covid/swabseq/EDSamples_073120.csv', stringsAsFactors=F, header=T, sep=',')
#v26

swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))




runs=paste0('v', c(22:36))
dfsL=list()
for(r in runs) {
    rundir=paste0(swabseq.dir, paste0('runs/',r, '/'))#v26/')
    # outdir=paste0(swabseq.dir, 'analysis/v26/')
    dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T, Stotal_filt=1000)
    df=dfL$df
    dfsL[[r]]=dfL$dfs
}
dfs=rbindlist(dfsL, idcol='experiment', fill=T)

da=dfs %>% filter(Description!="")  %>% filter(Description!=" ") %>% filter(Description!="0")

library(lme4)
library(lmerTest)
sm=lmer(log10(Stotal+1)~(1|Plate_ID)+(1|experiment)+(1|virus_copy)+S_index -1, data=da)
sm2=lmer(log10(Stotal+1)~(1|Plate_ID)+(1|experiment)+(1|virus_copy)+S_index2-1, data=da)
sm3=lmer(log10(Stotal+1)~(1|Plate_ID)+(1|experiment)+Sample-1, data=da)
sm4=lmer(log10(Stotal+1)~(1|Plate_ID)+(1|experiment)+Sample_Well-1, data=da)

i7=summary(sm,ddf='lme4')$coefficients
i5=summary(sm2,ddf='lme4')$coefficients
library(gdata)
par(mfrow=c(2,1))
plotCI(1:384, i7[,1],2*i7[,2], ylim=c(0,5.5),sfrac=0, gap=.1,
       main='i7 effect on log10(S2 + S2 spike)',
       ylab='log10(S2+S2spike)')
plotCI(1:384, i5[,1],2*i5[,2], ylim=c(0,5.5),sfrac=0, gap=.1,
       ylab='log10(S2+S2spike)',
       main='i5 effect on log10(S2 + S2 spike)')

iS=summary(sm3,ddf='lme4')$coefficients
plotCI(1:1536, iS[,1],2*iS[,2], ylim=c(0,7),sfrac=0, 
       ylab='log10(S2+S2spike)',
       main='i5 effect on log10(S2 + S2 spike)')

iSw=summary(sm4,ddf='lme4')$coefficients
plotCI(1:96, iSw[,1],2*iSw[,2], ylim=c(0,7),sfrac=0, 
       ylab='log10(S2+S2spike)',
       main='i5 effect on log10(S2 + S2 spike)')


x=summary(sm, ddf='Kenward-Roger')

da%>%ggplot(aes(x=Sample, y=Stotal))+
    #geom_point(alpha=.75, size=2)+
     #geom_quasirandom(alpha=.75, size=2)+
    #facet_wrap(~NPResult.x, scales='free_x')+
    stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange") +
  stat_summary(fun.y = mean,
               geom = "line")+
    scale_y_log10() + annotation_logticks(sides='l') +
    ylab('(Stotal)')+
    #ylim(c(0,10000))+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('S2+ S2 spike v22-v36')



edo=dfs[as.character(dfs$virus_identity) %in% ed.table$Matrix.Tube.Barcode | as.character(dfs$virus_identity) %in% ed.table$Deidenfied.ID,]
mu=match(edo$virus_identity, ed.table$Matrix.Tube.Barcode)
edo$virus_identity[!is.na(mu)]=na.omit(ed.table$Deidenfied.ID[mu])

edo=merge(edo, ed.table, all.x=T, by.x='virus_identity', by.y='Deidenfied.ID')

edo %>%  #filter(Stotal>1000) %>%  
    ggplot(aes(x=ct.char, y=S2_normalized_to_S2_spike,color=experiment, fill=Stotal>1000))+
    geom_point(alpha=.75, size=2)+
     #geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~NPResult.x, scales='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('ED')#+



all.ed=edo[,c(1,2,3,4,5,6,8,25:35,44,57)]

all.ed %>% write.csv('/data/Covid/swabseq/analysis/ED/allED.csv')

edo$Ct=as.numeric(edo$Ct.y)
edo$Ct[is.na(edo$Ct)]=40

edd=edo %>% filter(experiment=='v34') 
edd$Ct=as.numeric(edd$Ct.y)
edd$Ct[is.na(edd$Ct)]=40
edd$Ct[edd$Ct==0]=40
edd$Stot1000=edd$Stotal>1000
edd$NPResult=relevel(factor(edd$NPResult), 'ND')
edd%>%filter(Ct>0) %>% filter(Description!='729 ED' ) %>% filter(!grepl("\\(Should", X.5)) %>%
ggplot(aes(x=Ct, y=S2_normalized_to_S2_spike,label=virus_identity,color=Stot1000))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    geom_label_repel()+
    facet_grid(~NPResult, scales='free_x')+
    scale_x_reverse()+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
    #geom_smooth(method='lm')+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('ED Ratio vs Ct')





edd%>%filter(Ct>0) %>% filter(Description!='729 ED' ) %>% filter(!grepl("\\(Should", X.5)) %>% #filter(Stotal>1000) %>%
ggplot(aes(x=Ct, y=S2_normalized_to_S2_spike,label=virus_identity, color = Stotal>1000))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    #geom_label_repel()+
    facet_grid(~NPResult, scales='free_x')+
    scale_x_reverse()+
    scale_y_log10() + annotation_logticks(sides='l') + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
    #geom_smooth(method='lm')+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('ED Saliva vs NP study')



x=edd%>%filter(Ct>0) %>% filter(Description!='729 ED' ) %>% filter(!grepl("\\(Should", X.5)) %>% filter(Stotal>1000)
cor.test(x$Ct[x$NPResult=='DETECTED'], log10(x$S2_normalized_to_S2_spike[x$NPResult=='DETECTED']) )




edo%>%filter(Ct>0) %>%
ggplot(aes(x=Ct, y=S2_normalized_to_S2_spike,color=virus_identity))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~NPResult.x, scales='free_x')+
    scale_x_reverse()+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('ED Ratio vs Ct')

edo%>%
ggplot(aes(x=Ct, y=S2_normalized_to_S2_spike,color=experiment))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~Stot1000, scales='free_x')+
    scale_x_reverse()+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('ED Ratio vs Ct')




edo %>% filter(experiment=='v34') %>% filter(Stotal>1000) %>% filter(Ct>0) %>% 
    ggplot(aes(x=Ct, y=S2_normalized_to_S2_spike))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~Stot1000, scales='free_x')+
    scale_x_reverse()+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
    geom_smooth(method='lm')+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ggtitle('ED Ratio vs Ct , v34')


all.ed$SARS_COV_2_prob=predict(classifier, newdata=all.ed, type='response')
all.ed %>% write.csv(paste0(outdir, 'ED_rerun_logReg.csv'))











edd%>%filter(Ct>0) %>%filter(Description!='729 ED' ) %>% 
    #filter(Col!='06'|Col!='07'|Col!='08') %>%
    #filter(Stotal>1000) %>%
    ggplot(aes(x=virus_identity, y=S2_normalized_to_S2_spike,label= color=Stotal>1000))+
    #geom_point(alpha=.75, size=2)+
    geom_quasirandom(alpha=.75, size=2)+
    facet_wrap(~Description, scales='free_x')+
    scale_y_log10() + annotation_logticks() + ylab('(S2+1)/(S2 spike + 1)')+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



 cor.test(edd$Ct, log10(edd$S2_normalized_to_S2_spike), method='pearson')
c

tt=edo[edo$Ct>0,]



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

swabseq.dir='/data/Covid/swabseq/'
source(paste0(swabseq.dir, 'code/helper_functions.R'))
rundir=paste0(swabseq.dir, 'runs/v29/')
outdir=paste0(swabseq.dir, 'analysis/v29/')

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
    ggplot(aes(x=ct.char, y=S2_normalized_to_S2_spike, color=grepl('33', ct.char)))+
    facet_wrap(~NP_result, scales='free_x')+
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
