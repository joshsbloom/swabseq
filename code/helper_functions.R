library(tidyverse)
library(ggbeeswarm)
library(viridis)

mungeTables=function(tables.RDS,lw=F,Stotal_filt=2000,input=96){
    countTables=readRDS(tables.RDS) #paste0(rundir, 'countTable.RDS')) 
    df=do.call('rbind', countTables)
    df$virus_copy=as.factor(df$virus_copy) 
    df$Col=as.factor(gsub('^.', '', df$Sample_Well))
    if(input==96) {
    df$Row=factor(gsub('..$', '', df$Sample_Well), levels=rev(toupper(letters[1:8])))
    }
    if(input==384) {
    df$Row=factor(gsub('..$', '', df$Sample_Well), levels=rev(toupper(letters[1:16])))

    }
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
      right_join(dfs)
    #modify code to track indices
    s2.indices=dfs%>%filter(amplicon=='S2')%>%select(Sample, index,index2)
    names(s2.indices)[c(2,3)]=paste0('S_', names(s2.indices)[c(2,3)])
    rpp.indices=dfs%>%filter(amplicon=='RPP30')%>%select(Sample, index,index2)
    names(rpp.indices)[c(2,3)]=paste0('R_', names(rpp.indices)[c(2,3)])
    df.i=right_join(s2.indices, rpp.indices)
    
    dfs=dfs %>% 
      select(-mergedIndex, -Sample_ID, -index, -index2 ) %>% right_join(df.i) %>%
      spread(amplicon, Count) %>% 
      mutate(S2_normalized_to_S2_spike=(S2+1)/(S2_spike+1))%>%
      mutate(RPP30_Detected=RPP30>10) %>%  
      #filter(Plate_ID!='Plate8') %>%
      mutate(SARS_COV_2_Detected=S2_normalized_to_S2_spike>.003)
   
    dfs$SARS_COV_2_Detected[!dfs$RPP30_Detected]='Inconclusive'
    dfs$SARS_COV_2_Detected[dfs$Stotal<Stotal_filt]='Inconclusive' 
    if(lw) { return(list(df=df,dfs=dfs))} 
    return(dfs)
    
}


make_hamming1_sequences=function(x) {
    eseq=s2c(x)
    eseqs=c(c2s(eseq))
    for(i in 1:length(eseq)){
        eseq2=eseq
        eseq2[i]='A'
        eseqs=c(eseqs,c2s(eseq2))
        eseq2[i]='C'
        eseqs=c(eseqs,c2s(eseq2))
        eseq2[i]='T'
        eseqs=c(eseqs,c2s(eseq2))
        eseq2[i]='G'
        eseqs=c(eseqs,c2s(eseq2))
        eseq2[i]='N'
        eseqs=c(eseqs,c2s(eseq2))

    }
    eseqs=unique(eseqs)
    return(eseqs)
}


# error correct the indices and count amplicons using amatch
errorCorrectIdxAndCountAmpliconsOld=function(rid, count.table, ind1,ind2,e=1, threads=threads){
        # get set of unique expected index1 and index2 sequences
        index1=unique(count.table$index)
        index2=unique(count.table$index2)
        # for subset of reads where matching amplicon of interest (rid)
        # match observed index sequences to expected sequences allowing for e hamming distance
        i1m=amatch(ind1[rid],index1, method='hamming', maxDist=e, matchNA=F, nthread=threads)
        i2m=amatch(ind2[rid],index2, method='hamming', maxDist=e, matchNA=F, nthread=threads)
        # combine the error corrected indices together per read
        idm=paste0(index1[i1m], index2[i2m])
        #match error corrected indices to lookup table and count
        tS2=table(match(idm, count.table$mergedIndex))
        #get index in lookup table for indices with at least one observed count
        tbix=match(as.numeric(names(tS2)), 1:nrow(count.table))
        #increment these samples by count per sample
        count.table$Count[tbix]=as.vector(tS2)+count.table$Count[tbix]
        return(count.table)
}

#error correct the indices and count amplicons, equivalent base R, faster
errorCorrectIdxAndCountAmplicons=function(rid, count.table, ind1,ind2){
    # get set of unique expected index1 and index2 sequences
    index1=unique(count.table$index)
    index2=unique(count.table$index2)

    i1h=lapply(index1, make_hamming1_sequences)
    names(i1h)=index1
    ih1=Biobase::reverseSplit(i1h)
    ih1.elements=names(ih1)
    ih1.indices=as.vector(unlist(ih1))
    i1m=ih1.indices[match(ind1[rid],ih1.elements)]

    i2h=lapply(index2, make_hamming1_sequences)
    names(i2h)=index2
    ih2=Biobase::reverseSplit(i2h)
    ih2.elements=names(ih2)
    ih2.indices=as.vector(unlist(ih2))
    i2m=ih2.indices[match(ind2[rid],ih2.elements)]
    
    idm=paste0(i1m,i2m)
    tS2=table(match(idm, count.table$mergedIndex))
    tbix=match(as.numeric(names(tS2)), 1:nrow(count.table))
    count.table$Count[tbix]=as.vector(tS2)+count.table$Count[tbix]
    return(count.table)
}



#add Row384 and Col384 information for 384-well plates used
add384Mapping=function(df){
    #move this to helper functions
    col384=sprintf('%02d', 1:24)
    row384=toupper(letters[1:16])
    col384L=list(
    'A'=col384[seq(1,24,2)],
    'B'=col384[seq(2,24,2)],
    'C'=col384[seq(1,24,2)],
    'D'=col384[seq(2,24,2)])
    col384L=lapply(col384L, function(x) { names(x)=sprintf('%02d', 1:12); return(x); })

    row384L=list(
    'A'=row384[seq(1,16,2)],
    'B'=row384[seq(1,16,2)],
    'C'=row384[seq(2,16,2)],
    'D'=row384[seq(2,16,2)])
    row384L=lapply(row384L, function(x) { names(x)=toupper(letters[1:8]); return(x); })

    df$Row384=''
    df$Col384=''
    df$Row384[df$Plate_384_Quadrant=='A']=as.character(row384L[['A']][as.character(df$Row[df$Plate_384_Quadrant=='A'])])
    df$Row384[df$Plate_384_Quadrant=='B']=as.character(row384L[['B']][as.character(df$Row[df$Plate_384_Quadrant=='B'])])
    df$Row384[df$Plate_384_Quadrant=='C']=as.character(row384L[['C']][as.character(df$Row[df$Plate_384_Quadrant=='C'])])
    df$Row384[df$Plate_384_Quadrant=='D']=as.character(row384L[['D']][as.character(df$Row[df$Plate_384_Quadrant=='D'])])
    df$Col384[df$Plate_384_Quadrant=='A']=as.character(col384L[['A']][as.character(df$Col[df$Plate_384_Quadrant=='A'])])
    df$Col384[df$Plate_384_Quadrant=='B']=as.character(col384L[['B']][as.character(df$Col[df$Plate_384_Quadrant=='B'])])
    df$Col384[df$Plate_384_Quadrant=='C']=as.character(col384L[['C']][as.character(df$Col[df$Plate_384_Quadrant=='C'])])
    df$Col384[df$Plate_384_Quadrant=='D']=as.character(col384L[['D']][as.character(df$Col[df$Plate_384_Quadrant=='D'])])

    filled.plates=df %>% filter(Description!='' & Description!=' ') %>% filter(Plate_384!='')
    #filled.plates$Row384=droplevels(factor(filled.plates$Row384))
    #filled.plates$Col384=droplevels(factor(filled.plates$Col384))

    filled.plates$Row384=factor(filled.plates$Row384, levels=c(rev(toupper(letters[1:16]))))
    return(filled.plates)
}
