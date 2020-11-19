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
p=add_argument(p,"--lbuffer", default=30000001, help="how many reads to load into ram at a time")
p=add_argument(p,"--threads", default=8, help="number of threads for bcl2fastq")
p=add_argument(p,"--generateSampleSheet", default=F, help="use barcode csv file to generate sample sheet")

#load required packages
library(ShortRead)
library(stringdist)
library(tidyverse)
library(ggbeeswarm)
library(viridis)
library(reader)

args=parse_args(p) 
rundir=args$rundir
print(rundir)
basespaceID=args$basespaceID
outputCountsOnly=args$countsOnly
extendedAmplicons=args$extendedAmplicons
threads=args$threads
nbuffer=as.integer(as.numeric(args$lbuffer))
gSampleSheet=args$generateSampleSheet
#nbuffer=as.numeric(args$lbuffer) 13e7
print(nbuffer)

#extract bcl directory name ----------------------------------------------------------------------------
ff=find.file('RTAComplete.txt', dir=paste(rundir, list.files(rundir), '/', sep=''))
if(ff=="") {
    #Pull BCLs from basespace [skip this section if you already placed bcls in rundir/bcls/] ------------
    #if running miseq then paste run id here
    #if miseq run then grab from basespace, otherwise place bcls here and skip lines 12-23 
    if(basespaceID!=1000) {   system(paste("bs download run --id", basespaceID, "-o bcls/")) }
}

ff=find.file('RTAComplete.txt', dir=paste(rundir, list.files(rundir), '/', sep=''))
if(ff==""){
    print("bcl directory not found")
    quit(save="no")
}
ff=strsplit(ff,'/')[[1]]
bcl.dirname=ff[length(ff)-1]
bcl.dir=paste0(rundir,bcl.dirname,'/')
#---------------------------------------------------------------------------------------------------------

setwd(rundir)
source('../../code/helper_functions.R')
source('../../code/create_samplesheet.R')

if(gSampleSheet) { 
    #flowcell=unzipDirAndFixKey(rundir)
    #print(flowcell)
    makeSS(rundir,bcl.dir) 
}

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
#move to bcl directory
setwd(bcl.dir)

# if fastqs don't exist run bcl2fastq to generate them 
fastqR1  <- paste0(bcl.dir, 'out/Undetermined_S0_R1_001.fastq.gz')
if(!file.exists(fastqR1)) { 
      
    # run bcl2fastq to generate fastq.gz files (no demux is happening here)
    #note this is using 64 threads and running on a workstation, reduce threads if necessary
    system(paste("bcl2fastq --runfolder-dir . --output-dir out/ --create-fastq-for-index-reads  --ignore-missing-bcl --use-bases-mask=Y26,I10,I10 --processing-threads", threads, "--no-lane-splitting --sample-sheet /dev/null"))
    #for sharing, make tar of bcls directory
    system(paste("tar -cvf", paste0(rundir,'bcls.tar'),  bcl.dir))
    #-----------------------------------------------------------------------------------------------------
}

# read in n reads at a time (reduce if RAM limited)


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
if(sum(grepl('-1$', ss$Sample_ID))==0){
    ssS=ss
    ssR=ss
} else {
    ssS=ss[grep('-1$', ss$Sample_ID),]
    #subset of indices for RPP30
    ssR=ss[grep('-2$', ss$Sample_ID),]
}

#initalize output count tables ------------------------------------------------------------
count.tables=list()
for(a in names(amplicons)){
    if(grepl('^S',a)){    count.tables[[a]]=ssS   } 
    if(grepl('^R',a)){    count.tables[[a]]=ssR   } 
        count.tables[[a]]$Count=0
        count.tables[[a]]$amplicon=a
}
#--------------------------------------------------------------------------------------------

if(!file.exists(paste0(rundir, 'countTable.RDS'))) { 

    fastq_dir  <- paste0(bcl.dir, 'out/')
    in.fileI1  <- paste0(fastq_dir, 'Undetermined_S0_I1_001.fastq.gz')
    in.fileI2  <- paste0(fastq_dir, 'Undetermined_S0_I2_001.fastq.gz')
    in.fileR1  <- paste0(fastq_dir, 'Undetermined_S0_R1_001.fastq.gz')
            
    i1 <- FastqStreamer(in.fileI1, nbuffer, readerBlockSize = 1e9, verbose = T)
    i2 <- FastqStreamer(in.fileI2, nbuffer, readerBlockSize = 1e9, verbose = T)
    r1 <- FastqStreamer(in.fileR1, nbuffer, readerBlockSize = 1e9, verbose = T)

    amp.match.summary.table=rep(0, length(amplicons)+1)
    names(amp.match.summary.table)=c(names(amplicons),'no_align')

    repeat{
        rfq1 <- yield(i1) 
        if(length(rfq1) == 0 ) { break }
        rfq2 <- yield(i2) 
        rfq3 <- yield(r1) 
        ind1 <- sread(rfq1)
        ind2 <- sread(rfq2)
        rd1  <- sread(rfq3)
        
        #match amplicon
        amph1=lapply(amplicons, make_hamming1_sequences)
        amph1=Biobase::reverseSplit(amph1)
        amph1.elements=names(amph1)
        amph1.indices=as.vector(unlist(amph1))
        amp.match=amph1.indices[match(rd1, amph1.elements)]
        no_align=sum(is.na(amp.match))

        #summarize amplicon matches
        amp.match.summary=table(amp.match)
        amp.match.summary=amp.match.summary[match(names(amplicons),names(amp.match.summary))]
        amp.match.summary=c(amp.match.summary, no_align)
        names(amp.match.summary) <- c(names(amp.match.summary[-length(amp.match.summary)]),"no_align")
        amp.match.summary.table=amp.match.summary.table+amp.match.summary

        #convert to indices
        per.amplicon.row.index=lapply(names(amplicons), function(x) which(amp.match==x))
        names(per.amplicon.row.index)=names(amplicons)

        #for each amplicon of interest count up reads where indices match expected samples
        for(a in names(count.tables)){
          count.tables[[a]]= errorCorrectIdxAndCountAmplicons(per.amplicon.row.index[[a]], count.tables[[a]], ind1,ind2)
        }

    }
    close(i1); close(i2); close(r1);
    #results=list(S2.table=S2.table, S2_spike.table=S2_spike.table, RPP30.table=RPP30.table)
    #if(extendedAmplicons){ results=list(S2.table=S2.table, S2_spike.table=S2_spike.table, RPP30.table=RPP30.table,RPP30_spike.table=RPP30_spike.table) }
    names(count.tables)=paste0(names(count.tables), '.table')

    results=count.tables

    do.call('rbind', results) %>% write_csv(paste0(rundir, 'countTable.csv')) 
    saveRDS(results, file=paste0(rundir, 'countTable.RDS'),version=2)
    saveRDS(amp.match.summary.table, file=paste0(rundir, 'ampCounts.RDS'),version=2)
}

#BiocManager::install("Rqc")
#BiocManager::install("savR")
amp.match.summary.table=readRDS(paste0(rundir, 'ampCounts.RDS'))
results=readRDS(paste0(rundir, 'countTable.RDS'))

library(savR)
sav=savR(bcl.dir)
tMet=tileMetrics(sav)
phiX=mean(tMet$value[tMet$code=='300'])
clusterPF=mean(tMet$value[tMet$code=='103']/tMet$value[tMet$code=='102'], na.rm=T)
clusterDensity=mean(tMet$value[tMet$code=='100']/1000)
clusterDensity_perLane=sapply(split(tMet, tMet$lane), function(x) mean(x$value[x$code=='100']/1000))    
seq.metrics=data.frame("totalReads"=format(sum(amp.match.summary.table),  big.mark=','),
                       "totalReadsPassedQC"=format(sum(amp.match.summary.table[!(names(amp.match.summary.table) %in% 'no_align')]), big.mark=','),
                       "phiX"=paste(round(phiX,2), "%"), "clusterPF"=paste(round(clusterPF*100,1), "%"),
                       "clusterDensity"=paste(round(clusterDensity,1), 'K/mm^2'), 
                       "clusterDensity_perLane"=paste(sapply(clusterDensity_perLane, round,1),collapse=' '))

#library(Rqc)
#qcRes = rqc(path = paste0(bcl.dir, 'out/'), pattern = ".fastq.gz", openBrowser=FALSE, n=1e5, workers=1)
#read_quality <- rqcCycleQualityBoxPlot(qcRes) + ylim(0,NA)
#seq_cont_per_cycle <- rqcCycleBaseCallsLinePlot(qcRes)
#read_freq_plot <- rqcReadFrequencyPlot(qcRes)
#base_calls_plot <- rqcCycleBaseCallsLinePlot(qcRes)


empty_well_set=c('', 'TBET', 'No Tube', 'NO TUBE', 'Empty', 'EMPTY', ' ', 'NA')

setwd(rundir)
dfL=mungeTables(paste0(rundir, 'countTable.RDS'),lw=T, Stotal_filt=500, input=384)
dwide=data.frame(dfL$dfs)
rsample=!(dwide$virus_identity%in%empty_well_set | is.na(dwide$virus_identity))

#dwide %>% filter(grepl('^\\d', virus_identity)) %>% write.csv(paste0(rundir, 'report_patient_samples.csv'))
#dwide %>% filter(grepl('MNS NS', virus_identity)) %>% write.csv(paste0(rundir, 'report_MNS_NS_LOD.csv'))

loo.key=read.csv('../../reference/MatrixWaterTubes.csv', header=F, stringsAsFactors=F, colClasses=c('character', 'character', 'character'))
names(loo.key)=c('type','Pos96', 'ID')

#!='' & dwide$virus_identity!='TBET') #& dwide$virus_identity!='TBET') & dwide$virus_copy!='0' 
results.summary=data.frame(
'TotalSamplesPocessed'=sum(rsample),                         
'Inconclusives'=sum(dwide$SARS_COV_2_Detected[rsample]=='Inconclusive'),
'NoVirusDetected'=sum(dwide$SARS_COV_2_Detected[rsample]=='FALSE'),
'VirusDetected'=sum(dwide$SARS_COV_2_Detected[rsample]=='TRUE'))

params <- list(
               loo.key=loo.key,
        experiment = strsplit(rundir,"/") %>% unlist() %>% tail(1),
        bcl.dir = bcl.dirname,                   
        amp.match.summary = amp.match.summary.table,
        seq.metrics=seq.metrics,
        results = results,
        # qcRes = qcRes,
        #read_quality = read_quality,
        #seq_cont_per_cycle = seq_cont_per_cycle,
        #read_freq_plot = read_freq_plot,
        #base_calls_plot = base_calls_plot,
        results.summary = results.summary,
        dlong=add96Mapping(dfL$df),
        dwide=add96Mapping(dwide),
        dsummary.reduced=add96Mapping(dwide) %>% filter(rsample) %>%  select(virus_identity, Plate_96_BC,Pos96, quadrant_96, Plate_384,S2,S2_spike,RPP30,SARS_COV_2_Detected)
 )

dir.create(paste0(rundir, '/results/'))
#output summary table as csv
rsample2= !(dwide$virus_identity%in%empty_well_set | is.na(dwide$virus_identity) | dwide$virus_identity %in% params$loo.key$ID)

params$dwide %>% filter(rsample2)  %>% write.csv(paste0(rundir,'/results/', params$experiment,'_report.csv'))


rmarkdown::render(
        input = "../../code/qc_report.Rmd",
        output_file = paste0(params$experiment,".html"),
        output_dir = paste0(rundir, 'results/'),
        params = params,
        envir = new.env(parent = globalenv())
    )







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
