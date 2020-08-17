#usage: countAmplicons.R [--] [--help] [--opts OPTS] [--rundir RUNDIR]
#       [--basespaceID BASESPACEID] [--countsOnly COUNTSONLY]
#       [--extendedAmplicons EXTENDEDAMPLICONS]
#
#utility to count amplicons for SwabSeq
#
#flags:
#  -h, --help               show this help message and exit
#
#optional arguments:
#  -x, --opts               RDS file containing argument values
#  -r, --rundir             path containing SampleSheet [default: .]
#  -b, --basespaceID        BaseSpace Run ID [default: 1000]
#  -c, --countsOnly         only output table of counts [default: TRUE]
#  -e, --extendedAmplicons  additional swabseq amplicons [default:
#                           FALSE]


#parse command line arguments and exit fast if run with -h
library(argparser)
p = arg_parser("utility to count amplicons for SwabSeq")
p=add_argument(p,"--rundir",  default=".", help="path containing SampleSheet")
p=add_argument(p,"--basespaceID",  default=1000, help="BaseSpace Run ID")
#p=add_argument(p,"--countsOnly",  default=TRUE, help="only output table of counts")
p=add_argument(p,"--extendedAmplicons", default=F, help="additional swabseq amplicons")
#p=add_argument(p,"--lbuffer", default=30000001, help="how many reads to load into ram at a time")
#p=add_argument(p,"--threads", default=64, help="number of threads for bcl2fastq")
args=parse_args(p) 

#load required packages
library(ShortRead)
library(stringdist)
library(tidyverse)
library(ggbeeswarm)
library(viridis)
#for example ... 
#rundir='/data/Covid/swabseq/runs/SwabSeqV17/'
#basespaceID="195853687"
rundir=args$rundir
basespaceID=args$basespaceID
outputCountsOnly=args$countsOnly
extendedAmplicons=args$extendedAmplicons
#requires BCLs (or basespaceID), python3.7 in path, bcl2fastq in path and SampleSheet.csv (or xls)

# Generate sample sheet from XLS file ----------------------------------------------------------------
# using octant python script [skip this section if you already have SampleSheet.csv in rundir/
# see https://github.com/octantbio/SwabSeq
sampleXLS=paste0(rundir,'SwabSeq.xlsx')
sampleCSV=paste0(rundir,'SampleSheet.csv')
setwd(rundir)
# if csv sample sheet doesn't exist, make it 
if(!file.exists(sampleCSV)) {
    system(paste("python3.7 ../../code/platemap2samp.py", sampleXLS,  '--o', sampleCSV))
}
#-----------------------------------------------------------------------------------------------------

# if fastqs don't exist grab them from basespace
fastqR1  <- paste0(rundir, 'bcls/out/Undetermined_S0_R1_001.fastq.gz')
if(!file.exists(fastqR1)) { 
    #Pull BCLs from basespace [skip this section if you already placed bcls in rundir/bcls/] ------------
    #if running miseq then paste run id here
    #if miseq run then grab from basespace, otherwise place bcls here and skip lines 12-23 
    system(paste("bs download run --id", basespaceID, "-o bcls/"))
    # run bcl2fastq to generate fastq.gz files (no demux is happening here)
    setwd(paste0(rundir,'bcls/'))
    #note this is using 64 threads and running on a workstation, reduce threads if necessary
    system(paste("bcl2fastq --runfolder-dir . --output-dir out/ --create-fastq-for-index-reads  --ignore-missing-bcl --use-bases-mask=Y26,I10,I10 --processing-threads 8 --no-lane-splitting --sample-sheet /dev/null"))
    #for sharing, make tar of bcls directory
    system(paste("tar -cvf", paste0(rundir,'bcls.tar'), "../bcls/"))
    #-----------------------------------------------------------------------------------------------------
}

# read in n reads at a time (reduce if RAM limited)
nbuffer=3e7
#nbuffer=as.numeric(args$lbuffer)
#print(nbuffer)


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
if(extendedAmplicons) {
    amplicons=list(
        S2=      'TATCTTCAACCTAGGACTTTTCTATT',
        S2_spike='ATAGAACAACCTAGGACTTTTCTATT',
        RPP30='CGCAGAGCCTTCAGGTCAGAACCCGC',
        RPP30_spike='GCGTCAGCCTTCAGGTCAGAACCCGC'
    )
}

# Munginging sample sheet-------------------------------------------------------------------
ss=read.delim(paste0(rundir,'SampleSheet.csv'), stringsAsFactors=F, skip=14, sep=',')
ss$mergedIndex=paste0(ss$index, ss$index2)

# this code would be obviated if indices designate wells, for most analyses here there are different indices for s2/s2spike and rpp30
# subset of indices for S2/S2 spike
ssS=ss[grep('-1$', ss$Sample_ID),]
#subset of indices for RPP30
ssR=ss[grep('-2$', ss$Sample_ID),]

#initalize output count tables ------------------------------------------------------------
S2.table=ssS; S2_spike.table=ssS; RPP30.table=ssR; RPP30_spike.table=ssR;
S2.table$Count=0; S2_spike.table$Count=0; RPP30.table$Count=0; RPP30_spike.table$Count=0
S2.table$amplicon='S2'; S2_spike.table$amplicon='S2_spike'; RPP30.table$amplicon='RPP30'; RPP30_spike.table$amplicon='RPP30_spike'
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
    if(extendedAmplicons){
        RPP30_spike.table    = errorCorrectIdxAndCountAmplicons(per.amplicon.row.index$RPP30_spike, RPP30_spike.table, ind1,ind2,1)
    }

}
close(i1); close(i2); close(r1);
results=list(S2.table=S2.table, S2_spike.table=S2_spike.table, RPP30.table=RPP30.table)
if(extendedAmplicons){ results=list(S2.table=S2.table, S2_spike.table=S2_spike.table, RPP30.table=RPP30.table,RPP30_spike.table=RPP30_spike.table) }

do.call('rbind', results) %>% write_csv(paste0(rundir, 'countTable.csv')) 
saveRDS(results, file=paste0(rundir, 'countTable.RDS'),version=2)

# this is moved to mungeTables() in helper_functions.R
#if(outputCountsOnly==FALSE) {
#    results=readRDS(paste0(rundir, 'countTable.RDS'))
#    attach(results)
#    # re-formatting ... output count table 
#    df=rbind(S2.table,S2_spike.table,RPP30.table)
#    df$virus_copy=as.factor(df$virus_copy) 
#    df$Col=as.factor(gsub('^.', '', df$Sample_Well))
#    df$Row=factor(gsub('..$', '', df$Sample_Well), levels=rev(toupper(letters[1:8])))
#    df$Sample=paste0(df$Plate_ID, '-' ,df$Sample_Well)
#    df$Plate_ID=as.factor(df$Plate_ID)
#    df$Plate_ID=factor(df$Plate_ID, levels(df$Plate_ID)[order(as.numeric(gsub('Plate', '', levels(df$Plate_ID))))])  
#    df$Plate_384=as.factor(df$Plate_384)
#    #to do update this for extended amplicon set
#    df$amplicon=factor(df$amplicon, level=c('S2', 'S2_spike', 'RPP30'))
#    df %>%
#      write_csv(paste0(rundir, 'annotated_df.csv')) 
#
#
#    #plate visualization 
#    df %>%
#      ggplot(aes(x=Col, y=Row, fill=log10(Count))) + 
#      geom_raster() +
#      coord_equal() +
#      facet_grid(amplicon~Plate_384+Plate_ID) +
#      scale_fill_viridis_c(option = 'plasma')
#    ggsave(paste0(rundir,'plateVis.png'))
#
#
#
#    #assay results
#    dfs= df %>%filter(amplicon=='S2') %>%  
#      count(Sample_Well, wt=Count, name='S2_total_across_all_wells') %>%
#      right_join(df)
#    dfs= df %>%filter(amplicon=='S2'|amplicon=='S2_spike') %>%  
#      count(Sample, wt=Count, name='Stotal') %>%
#      right_join(dfs)
#    dfs= dfs %>% count(Sample, wt=Count, name='well_total') %>%
#      right_join(dfs) %>% 
#      select(-mergedIndex, -Sample_ID, -index, -index2 ) %>% 
#      spread(amplicon, Count) %>% 
#      mutate(S2_normalized_to_S2_spike=(S2+1)/(S2_spike+1))%>%
#      mutate(RPP30_Detected=RPP30>10) %>%  
#      mutate(SARS_COV_2_Detected=S2_normalized_to_S2_spike>.003)
#    dfs$SARS_COV_2_Detected[!dfs$RPP30_Detected]='Inconclusive'
#    dfs$SARS_COV_2_Detected[dfs$Stotal<2000]='Inconclusive'
#    dfs %>%
#      write_csv(paste0(rundir, 'Calls_per_sample.csv')) 
#}
