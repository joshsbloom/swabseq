data for each and the structure of each subfolder is:
```
runs/
├── v*
│   ├── bcl_link.txt
│   ├── bcls
│   ├── bcls.tar
│   ├── countTable.csv
│   ├── countTable.RDS
│   ├── ExperimentSetup.xlsx
│   ├── SampleSheet.csv
│   └── SwabSeq.xlsx
```

each subdirectory contains:
+ bcl_link.txt a link to google drive .tar file containing bcls and undemuxed fastq.gz files (untar bcl folder in each subdirectory, fastqs are in bcls/out)
+ ExperimentSetup.xlsx an excel file describing experiment setup and goals
+ SwabSeq*.xlsx an excel file formatted for [platemap2samp.py](code/platemap2samp.py) for generating csv sample sheet
+ SampleSheet.csv a properly formatted sample sheet for Illumina instruments, also used to keep track of experiment variables for each sample
+ countTable.csv a table of counts per amplicon
+ countTable.RDS a version 2 RDS file of the counts per amplicon, loaded by scripts in [analysis](../analysis/)

