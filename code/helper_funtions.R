library(tidyverse)
library(ggbeeswarm)
library(viridis)

mungeTables=function(tables.RDS,lw=F){
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
    if(lw) { return(list(df=df,dfs=dfs))} 
    return(dfs)
    
}

