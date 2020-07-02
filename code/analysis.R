library(ShortRead)
library(stringdist)
library(tidyverse)
library(ggbeeswarm)
library(viridis)

#for example ... 
rundir='/data/Covid/swabseq/runs/SwabSeqV17/'
basespaceID="195853687"

#requires BCLs (or basespaceID), python3.7 in path, bcl2fastq in path and SampleSheet.csv (or xls)

# Generate sample sheet from XLS file ----------------------------------------------------------------
# using octant python script [skip this section if you already have SampleSheet.csv in rundir/
# see https://github.com/octantbio/SwabSeq
sampleXLS=paste0(rundir,'SwabSeq.xlsx')
sampleCSV=paste0(rundir,'SampleSheet.csv')
setwd(rundir)
system(paste("python3.7 ../code/platemap2samp.py", sampleXLS,  '--o', sampleCSV))
#-----------------------------------------------------------------------------------------------------

#Pull BCLs from basespace [skip this section if you already placed bcls in rundir/bcls/] ------------
#if running miseq then paste run id here
#if miseq run then grab from basespace, otherwise place bcls here and skip lines 12-23 
system(paste("bs download run --id", basespaceID, "-o bcls/"))
# run bcl2fastq to generate fastq.gz files (no demux is happening here)
setwd(paste0(rundir,'bcls/'))
#note this is using 64 threads and running on a workstation, reduce threads if necessary
system("bcl2fastq --runfolder-dir . --output-dir out/ --create-fastq-for-index-reads  --ignore-missing-bcl --use-bases-mask=Y26,I10,I10 --processing-threads 64 --no-lane-splitting --sample-sheet /dev/null")
#-----------------------------------------------------------------------------------------------------

# read in n reads at a time (reduce if RAM limited)
nbuffer=3e7

# error correct the indices and count amplicons
errorCorrectIdxAndCountAmplicons=function(rid, count.table, ind1,ind2,e=1){
        # get set of unique expected index1 and index2 sequences
        index1=unique(count.table$index)
        index2=unique(count.table$index2)
        # for subset of reads where matching amplicon of interest (rid)
        # match observed index sequences to expected sequences allowing for e hamming distance
        i1m=amatch(ind1[rid],index1, method='hamming', maxDist=e, matchNA=F)
        i2m=amatch(ind2[rid],index2, method='hamming', maxDist=e, matchNA=F)
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

#expected amplicons, note will update for RPP30 spike-in 
#RPP
#CGCAGA gccttcaggtcagaacccgc
#RPP spike
#GCGTCA gccttcaggtcagaacccgc
amplicons=list(
    S2=      'TATCTTCAACCTAGGACTTTTCTATT',
    S2_spike='ATAGAACAACCTAGGACTTTTCTATT',
    RPP30='CGCAGAGCCTTCAGGTCAGAACCCGC'
)

# Munginging sample sheet-------------------------------------------------------------------
ss=read.delim(paste0(rundir,'SampleSheet.csv'), stringsAsFactors=F, skip=14, sep=',')
ss$mergedIndex=paste0(ss$index, ss$index2)

# subset of indices for S2/S2 spike
ssS=ss[grep('-1$', ss$Sample_ID),]
#subset of indices for RPP30
ssR=ss[grep('-2$', ss$Sample_ID),]

#initalize output count tables ------------------------------------------------------------
S2.table=ssS; S2_spike.table=ssS; RPP30.table=ssR
S2.table$Count=0; S2_spike.table$Count=0; RPP30.table$Count=0
S2.table$amplicon='S2'; S2_spike.table$amplicon='S2_spike'; RPP30.table$amplicon='RPP30'
#------------------------------------------------------------------------------------------

fastq_dir  <- paste0(rundir, 'bcls/out/')
in.fileI1  <- paste0(fastq_dir, 'Undetermined_S0_I1_001.fastq.gz')
in.fileI2  <- paste0(fastq_dir, 'Undetermined_S0_I2_001.fastq.gz')
in.fileR1  <- paste0(fastq_dir, 'Undetermined_S0_R1_001.fastq.gz')
        
i1 <- FastqStreamer(in.fileI1, nbuffer, readerBlockSize = 1e9, verbose = T)
i2 <- FastqStreamer(in.fileI2, nbuffer, readerBlockSize = 1e9, verbose = T)
r1 <- FastqStreamer(in.fileR1, nbuffer, readerBlockSize = 1e9, verbose = T)

repeat{
    rfq1 <- yield(i1) 
    if(length(rfq1) == 0 ) { break }
    rfq2 <- yield(i2) 
    rfq3 <- yield(r1) 
    ind1 <- sread(rfq1)
    ind2 <- sread(rfq2)
    rd1  <- sread(rfq3)
                       
    #match amplicon
    amp.match=lapply(amplicons, function(x) {amatch(rd1, x, method='hamming', maxDist=1, matchNA=F)})

    #summary
    print(sapply(amp.match, function(x) sum(!is.na(x))))

    #convert to indices
    am.mat=do.call('cbind', amp.match)
    per.amplicon.row.index=apply(!is.na(am.mat), 2, function(x) which(x==TRUE))
    
    #for each amplicon of interest count up reads where indices match expected samples
    S2.table       = errorCorrectIdxAndCountAmplicons(per.amplicon.row.index$S2, S2.table, ind1,ind2, 1)
    S2_spike.table = errorCorrectIdxAndCountAmplicons(per.amplicon.row.index$S2_spike, S2_spike.table, ind1,ind2,1)
    RPP30.table    = errorCorrectIdxAndCountAmplicons(per.amplicon.row.index$RPP30, RPP30.table, ind1,ind2,1)
}
close(i1); close(i2); close(r1);


# re-formatting ... output count table 
df=rbind(S2.table,S2_spike.table,RPP30.table)
df$virus_copy=as.factor(df$virus_copy) 
df$Col=as.factor(gsub('^.', '', df$Sample_Well))
df$Row=factor(gsub('..$', '', df$Sample_Well), levels=rev(toupper(letters[1:8])))
df$Sample=paste0(df$Plate_ID, '-' ,df$Sample_Well)
df$Plate_ID=as.factor(df$Plate_ID)
df$Plate_ID=factor(df$Plate_ID, levels(df$Plate_ID)[order(as.numeric(gsub('Plate', '', levels(df$Plate_ID))))])  
df$Plate_384=as.factor(df$Plate_384)
df$amplicon=factor(df$amplicon, level=c('S2', 'S2_spike', 'RPP30'))
df %>%
  write_csv(paste0(rundir, 'annotated_df.csv')) 


#plate visualization 
df %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
  geom_raster() +
  coord_equal() +
  facet_grid(amplicon~Plate_384+Plate_ID) +
  scale_fill_viridis_c(option = 'plasma')
ggsave(paste0(rundir,'plateVis.png'))


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
  mutate(S2_normalized_to_S2_spike=(S2)/(S2_spike))%>%
  mutate(RPP30_Detected=RPP30>10) %>%  
  mutate(SARS_COV_2_Detected=S2_normalized_to_S2_spike>.003)
dfs$SARS_COV_2_Detected[!RPP30_Detected]='Inconclusive'
dfs$SARS_COV_2_Detected[dfs$Stotal<2000]='Inconclusive'
dfs %>%
  write_csv(paste0(rundir, 'Calls_per_sample.csv')) 
