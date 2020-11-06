library(plater)
library(XML)
library(seqinr)
library(data.table)
revcomp=function (x) {    toupper(c2s(rev(comp(s2c(x)))))}

makeSS=function(rundir, bcl.dir){
    #rundir='/data/Covid/swabseq/runs/vSS/'
    #bcl.dir='/data/Covid/swabseq/runs/vSS/201016_MN01365_0006_A000H37HTT/'
    #bcl.dirname='201016_MN01365_0006_A000H37HTT'
    xmlinfo=xmlToList(xmlParse(paste0(bcl.dir, 'RunParameters.xml')))
    chemistry=xmlinfo$Chemistry
    #MiniSeq High / MiniSeq Rapid High / NextSeq Mid / NextSeq High
    instrument=xmlinfo$ComputerName
    #MINISEQ / MINISEQA . /NB552456
    if(instrument=='MINISEQ')
        instrument_type='MiniSeq'
    if(instrument=='MINISEQA')
        instrument_type='MiniSeq'
    if(instrument=='NB552456')
        instrument_type='NextSeq'

    #don't reverse comp i5 for miniseq rapid or miseq / (note verify that chemistry toggle exists for xml for miseq, untested)
    i5RC.toggle=TRUE
    if(chemistry=="MiniSeq Rapid High" | chemistry=="MiSeq") {i5RC.toggle=F} 


    # look for zip file from box and unzip
    zipfile=list.files(rundir, pattern='.zip', full.names=F)
    if(identical(zipfile, character(0))){
       print(".zip directory not found")
        quit(save="no")
    }
    setwd(rundir)
    system(paste('unzip', zipfile))

    # locate key file 
    zipdir=paste0(rundir, gsub('.zip', '/', zipfile))
    setwd(zipdir)
    key.file=system('grep -iRl ",,1,1,2,2,3,3,4,4,5,5"', intern=T)
    if(identical(key.file, character(0))){
        print(".csv key file not found")
        quit(save="no")
    }

    #remove first column and row
    new.key.file=paste0(zipdir, 'keyfile.csv')
    system(paste("sed 1d", key.file, "| cut -f1 -d ',' --complement - >", new.key.file))

    # setup indices -----------------------------------------------------
    setwd(rundir)
    i7s=read_plates('../../reference/s2_r.csv', well_ids_column="Sample_Well")
    i7s=gather(i7s, Plate_ID, index, 3:ncol(i7s))
    i7s$index=as.vector(sapply(i7s$index, revcomp))
    i7s$Plate_ID=paste0('Plate', i7s$Plate_ID)
    i7s$Sample_ID=paste0(i7s$Plate_ID,'-', i7s$Sample_Well)

    i5s=read_plates('../../reference/s2_f.csv', well_ids_column="Sample_Well")
    i5s=gather(i5s, Plate_ID, index2, 3:ncol(i5s))
    if(i5RC.toggle) { i5s$index2=as.vector(sapply(i5s$index2, revcomp)) }
    i5s$Plate_ID=paste0('Plate', i5s$Plate_ID)
    i5s$Sample_ID=paste0(i5s$Plate_ID,'-', i5s$Sample_Well)

    i5s= i5s %>% select(Sample_ID, index2)
    i7s$bc_set='N1_S2_RPP30'

    index.key=right_join(i7s,i5s,by='Sample_ID') %>% select(-Plate) %>% select(Plate_ID,bc_set,Sample_Well,Sample_ID,index,index2)

    qA=apply(expand.grid(toupper(letters[seq(1,16,2)]), sprintf("%02d",seq(1,24,2))),1, paste0, collapse='')
    qB=apply(expand.grid(toupper(letters[seq(1,16,2)]), sprintf("%02d",seq(2,24,2))),1, paste0, collapse='')
    qC=apply(expand.grid(toupper(letters[seq(2,16,2)]), sprintf("%02d",seq(1,24,2))),1, paste0, collapse='')
    qD=apply(expand.grid(toupper(letters[seq(2,16,2)]), sprintf("%02d",seq(2,24,2))),1, paste0, collapse='')
    index.key$quadrant_96=''
    index.key$quadrant_96[index.key$Sample_Well %in% qA]='A'
    index.key$quadrant_96[index.key$Sample_Well %in% qB]='B'
    index.key$quadrant_96[index.key$Sample_Well %in% qC]='C'
    index.key$quadrant_96[index.key$Sample_Well %in% qD]='D'
    #-----------------------------------------------------------------

    # setup experiment key ---------------------------------------------------------
    #plate_ID will always refer to index plate location
    #Plate_384 will be a sample specific name
    experiment.key=read_plates(new.key.file, well_ids_column='Sample_Well')
    experiment.key=gather(experiment.key, Plate_384_long, virus_identity, 3:ncol(experiment.key))
    experiment.key$virus_identity[experiment.key$virus_identity=='#REF!']=''
    experiment.key$virus_identity=gsub('\"', '', experiment.key$virus_identity) 
    experiment.key$Plate_384_long=gsub("::", ": : ", experiment.key$Plate_384_long)

    #plate #, barcode, primer set"
    #PLATE 1:12346:1:TS01399002:TS01399030::
    ekd=tstrsplit(experiment.key$Plate_384_long, ':', names=c('Plate_384','Plate_384_BC', 'Plate_ID', 'Q1', 'Q2','Q3','Q4'))
    ekd$Plate_ID=paste0('Plate', ekd$Plate_ID)
    ekd$Q1[ekd$Q1==" "]=""
    ekd$Q2[ekd$Q2==" "]=""
    ekd$Q3[ekd$Q3==" "]=""
    ekd$Q4[ekd$Q4==" "]=""
    experiment.key$Plate_ID=ekd$Plate_ID
    experiment.key$Sample_ID=paste0(experiment.key$Plate_ID,'-',experiment.key$Sample_Well)
    experiment.key$Plate_384=ekd$Plate_384
    experiment.key$Plate_384_BC=ekd$Plate_384_BC
    experiment.key$Plate_96_BC=''
    experiment.key$Q1=ekd$Q1
    experiment.key$Q2=ekd$Q2
    experiment.key$Q3=ekd$Q3
    experiment.key$Q4=ekd$Q4

    experiment.key= experiment.key %>%select(-Plate, -Sample_Well, -Plate_ID) %>% right_join(index.key, by='Sample_ID')
    head(data.frame(experiment.key))

    for(p in unique(experiment.key$Plate_384)){
        experiment.key$Plate_96_BC[((experiment.key$Plate_384 %in% p) & (experiment.key$quadrant_96 %in% 'A'))]=
        experiment.key$Q1[((experiment.key$Plate_384 %in% p) & (experiment.key$quadrant_96 %in% 'A'))]
        
        experiment.key$Plate_96_BC[((experiment.key$Plate_384 %in% p) & (experiment.key$quadrant_96 %in% 'B'))]=
        experiment.key$Q2[((experiment.key$Plate_384 %in% p) & (experiment.key$quadrant_96 %in% 'B'))]

        experiment.key$Plate_96_BC[((experiment.key$Plate_384 %in% p) & (experiment.key$quadrant_96 %in% 'C'))]=
        experiment.key$Q3[((experiment.key$Plate_384 %in% p) & (experiment.key$quadrant_96 %in% 'C'))]

        experiment.key$Plate_96_BC[((experiment.key$Plate_384 %in% p) & (experiment.key$quadrant_96 %in% 'D'))]=
        experiment.key$Q4[((experiment.key$Plate_384 %in% p) & (experiment.key$quadrant_96 %in% 'D'))]
    }

    sample_sheet_info=(experiment.key %>% select(Plate_ID,virus_identity,bc_set,quadrant_96, Plate_96_BC, Plate_384, Plate_384_BC, Sample_Well,Sample_ID,index,index2))


    sample.sheet.csv=paste0(rundir, 'SampleSheet.csv')
    sink(sample.sheet.csv)
    cat("[Header]",
    "IEMFileVersion,5",
    "Investigator Name,SwabSeq",
    "Experiment Name, ",
    paste("Date,", date()),
    "Workflow,GenerateFASTQ",
    "Application,FASTQ Only",
    paste("Instrument Type,",instrument_type),
    "Chemistry,Amplicon",
    "[Reads]",
    "26",
    "[Settings]",
    "",
    "[Data]",sep="\n")
    sink()
    write.table(sample_sheet_info, file=sample.sheet.csv, append=T, quote=F, row.names=F, sep=',')


    return(NULL)


}
