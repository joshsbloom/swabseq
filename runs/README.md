the structure of each subfolder is:
```
runs/
├── v...
│   ├── bcl_link.txt
│   ├── bcls
│   ├── bcls.tar
│   ├── countTable.csv
│   ├── countTable.RDS
│   ├── ExperimentSetup.xlsx
│   ├── SampleSheet.csv
│   └── SwabSeqv10.xlsx
```

each directory contains:
+ bcl_link.txt a link to google drive .tar file containing bcls and undemuxed fastq.gz files
+ ExperimentSetup.xlsx an excel file describing experiment setup and goals
+ SwabSeq...xlsx an excel file formatted for [platemap2samp.py](code/platemap2samp.py) for generating csv sample sheet
+ SampleSheet.csv a properly formatted sample sheet for Illumina instruments, also used to keep track of experiment variables for each sample
+ countTable.csv a table of counts per amplicon
+ countTable.RDS a version 2 RDS file of the counts per amplicon, loaded by scripts in [analysis](analysis/)





Statistics for each sequencing run 
___

Run | Instrument | Read | Cluster PF (%) | % ≥Q30 | Yield | Error Rate% | Reads PF | Density  | Tiles | Legacy Phas/Prephas (%) | Intensity | PhiX Observerd | PhiX Targeted | Conc. Loaded
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- 
v7 | MiSeq | 1 | 92.71±0.99 | 97.26 | 218.65Mbp | 0.67±0.32 | 8,745,975 | 363±31 | 38 | 0.102/0.044 | 190±34 | 37.28 | 
v8 | MiSeq | 1 | 95.06±0.83 | 98.29 | 475.38Mbp | 0.16±0.02 | 19,015,232 | 776±11 | 38 | 0.084/0.043 | 146±18 | 38.61 |
v9 | MiSeq | 1 | 95.47±0.89 | 98.74 | 346.55Mbp | 0.48±0.08 | 13,862,051 | 563±20 | 38 | 0.095/0.029 | 172±26 | 28.37 |
v10 | MiSeq | 1 | 94.33±0.93 | 97.94 | 440.62Mbp | 0.19±0.04 | 17,624,768 | 727±15 | 38 | 0.103/0.086 | 112±16 | 35.50 |
v11 | MiSeq | 1 | 94.14±1.02 | 98.37 | 423.37Mbp | 0.52±0.13 | 16,934,964 | 703±31 | 38 | 0.090/0.064 | 152±21 | 30.99 | 
v12 | MiSeq | 1 | 93.18±0.96 | 97.35 | 171.26Mbp | 0.35±0.11 | 6,850,514 | 277±24 | 38 | 0.166/0.105 | 147±16 | 44.20 | 
v13 | MiSeq | 1 | 97.16±0.97 | 98.96 | 274.97Mbp | 0.14±0.01 | 10,998,699 | 406±21 | 38 | 0.137/0.101 | 193±30 | 15.26 | 
v14 | MiSeq | 1 | 97.00±0.44 | 99.01 | 463.36Mbp | 0.12±0.01 | 18,534,400 | 742±13 | 38 | 0.138/0.101 | 203±24 | 17.34 | 
v15 | MiSeq | 1 | 75.37±7.07 | 94.56 | 497.81Mbp | 0.83±0.51 | 19,912,204 | 1,103±45 | 38 | 0.027/0.078 | 177±20 | 17.77 |
v16 | MiSeq | 1 | 94.99±1.07 | 98.20 | 598.33Mbp | 0.22±0.34 | 23,933,240 | 973±24 | 38 | 0.032/0.000 | 193±31 | 16.72 | 
v17 | MiSeq | 1 | 93.83±1.71 | 97.63 | 555.45Mbp | 0.13±0.01 | 22,218,176 | 934±12 | 38 | 0.006/0.036 | 175±21 | 16.75 | 
v18 | MiSeq | 1 | 69.46±22.02 | 94.37 | 438.81Mbp | 1.25±0.62 | 17,552,240 | 1,042±57 | 38 | 0.008/0.000 | 192±22 | 11.46 | 
v19 | MiSeq | 1 | 85.26±20.73 | 97.28 | 563.71Mbp | 0.50±0.49 | 22,548,216 | 1,052±47 | 38 | 0.011/0.000 | 189±27 | 13.58 | 
v20 | MiSeq | 1 | 96.09±0.91 | 98.34 | 284.75Mbp | 0.27±0.35 | 11,390,135 | 443±9 | 38 | 0.270/0.020 | 192±31 | 10.94 | 
v21 | NextSeq | 1 | 95.54±0.45 | 98.48 | 0.17Gbp | 0.36±0.10 | 27,807,589 | 33±1 | 72 | 0.119/0.183 | 12741±809 | 8.44 | 
v22 | MiSeq | 1 | 96.35±0.88 | 98.57 | 347.74Mbp | 0.20±0.31 | 13,909,713 | 560±20 | 38 | 0.135/0.107 | 195±34 | 14.84 |
