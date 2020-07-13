# UCLA Amplicon Quantification Software for analyzing amplicon data from SwabSeq Sars-CoV-2 testing
___
see [countAmplicons.R](code/countAmplicons.R) for script to turn BCLs and XLSX description of experiment into tables of counts per amplicon

see [analysis](analysis/) for run-specific analysis scripts 

see [runs](runs/) for BCLs, gzipped fastqs and tables of amplicon counts for each experiment
___

Statistics and notes for each sequencing run 
___

Run | Instrument | % PhiX Observed | % PhiX Targeted | Conc. Loaded (pM) | RT | Indexing Strategy | PCR cycles | Notes | Read | Cluster PF (%) | % ≥Q30 | Yield | Error Rate% | Reads PF | Density  | Tiles | Legacy Phas/Prephas (%) | Intensity 
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- 
[v25](analysis/v25/) |  MiSeq | 36 | 40 | 25 | Taqpath | UDI | 50 | Saliva TBE+/tween; TE,RNASecure,Qiagen Protease; +/- 95 C heating 
[v24](analysis/v24/) |  MiSeq |  42.21 | 40  | 25 | Taqpath | Semi-Combinatorial | 50 | Saliva +/- RNA Secure; NP in VTM, NS, or Amies; dilutions in water |  1  |	92.85±1.17 | 97.79 |  	483.40 Mbp 	| 0.14±0.04 |	19,336,026 | 	829±18 | 	38 |	0.014/0.040| 	183 ±28
v23 |  MiSeq | 37.41  | 40 | 26 | Taqpath | UDI | 50 | NP in VTM diluted; eSwab(aimes); Saliva in cells2cdna; NP in NS diluted; increased lysate vol to 10uL; dilutions in TE | 1 |	92.69±0.86 | 	97.88 | 	570.91 Mbp | 	0.14±0.03 |	22,836,378 | 971±15 | 38 |	0.000/0.000 |	185 ±26
v22 |  MiSeq   | 14.84 | 40   | 25 | Taqpath | Semi-Combinatorial | 40 and 50 | Saliva undiluted; NP into VTM +/- dilution; v14 rerun;                                            | 1 | 96.35±0.88 | 98.57 | 347.74Mbp | 0.20±0.31 | 13,909,713 | 560±20 | 38 | 0.135/0.107 | 195±34
v21 |  NextSeq | 8.44   | 35  | 1.5 | Taqpath | UDI | 40 | NextSeq; mix of v18,v19,v20                                                                                                         | 1 | 95.54±0.45 | 98.48 | 0.17Gbp | 0.36±0.10 | 27,807,589 | 33±1 | 72 | 0.119/0.183 | 12741±809 
v20 |  MiSeq   | 10.94  | 35  | 24 | Taqpath | UDI | 40 | EUA (LoD confirmation);  Saliva undiluted                                                                                            | 1 | 96.09±0.91 | 98.34 | 284.75Mbp | 0.27±0.35 | 11,390,135 | 443±9 | 38 | 0.270/0.020 | 192±31 
v19 |  MiSeq | 13.58   | 35 | 26 | Taqpath | UDI | 40 | EUA (prelimary LoD and Positive and Negative patient samples)                                                                   | 1 | 85.26±20.73 | 97.28 | 563.71Mbp | 0.50±0.49 | 22,548,216 | 1,052±47 | 38 | 0.011/0.000 | 189±27       
v18 |  MiSeq | 11.46   | 30 | 27 | Taqpath | Semi-Combinatorial | 40 | Negative Patient and attempted contrived EUA (high basecalling error rate, diagnostic for amplicon misassignment)| 1 | 69.46±22.02 | 94.37 | 438.81Mbp | 1.25±0.62 | 17,552,240 | 1,042±57 | 38 | 0.008/0.000 | 192±22       
v17 |  MiSeq | 16.75   | 30 | 26 | Taqpath | Semi-Combinatorial| 40 | Negative Patient and attempted contrived EUA (but too high viral input for LoD experiment)                        | 1 | 93.83±1.71 | 97.63 | 555.45Mbp | 0.13±0.01 | 22,218,176 | 934±12 | 38 | 0.006/0.036 | 175±21             
v16 |  MiSeq | 16.72   | 30 | 31.3 | Taqpath | Combinatorial    | 40 | v13/v14 and Simulated Patients; first test no RPP30 unindexed                                                    | 1 | 94.99±1.07 | 98.20 | 598.33Mbp | 0.22±0.34 | 23,933,240 | 973±24 | 38 | 0.032/0.000 | 193±31             
v15 |  MiSeq | 17.77   | 30 | 40.8 | Taqpath | Combinatorial  | 50 | rerun V13 with 384-well plates 1-3; bleach wash                                                                    | 1 | 75.37±7.07 | 94.56 | 497.81Mbp | 0.83±0.51 | 19,912,204 | 1,103±45 | 38 | 0.027/0.078 | 177±20         
v14 |  MiSeq | 17.34 | 30 | 40.8 | Taqpath | Combinatorial | 50 | Saliva in NS diluted; Contrived in HEK extracted; no bleach wash  | 1 | 97.00±0.44 | 99.01 | 463.36Mbp | 0.12±0.01 | 18,534,400 | 742±13 | 38 | 0.138/0.101 | 203±24 
v13 |  MiSeq | 15.26 | 30 | 40.8 | Taqpath | Combinatorial | 50 | Saliva in NS diluted; Contrived in HEK extracted; Titrate Rxn Vols| 1 | 97.16±0.97 | 98.96 | 274.97Mbp | 0.14±0.01 | 10,998,699 | 406±21 | 38 | 0.137/0.101 | 193±30 
v12 |  MiSeq | 44.20 | 40 | 20 | Taqpath | UDI | 50 | MTS in TE diluted and titration of RPP30 unindexed primers                    | 1 | 93.18±0.96 | 97.35 | 171.26Mbp | 0.35±0.11 | 6,850,514 | 277±24 | 38 | 0.166/0.105 | 147±16   
v11 |  MiSeq | 30.99 | 40 | 22 | Taqpath | UDI | 50 | MTS in TE diluted                                                             | 1 | 94.14±1.02 | 98.37 | 423.37Mbp | 0.52±0.13 | 16,934,964 | 703±31 | 38 | 0.090/0.064 | 152±21 
v10 |  MiSeq | 35.50 | 40 | 20 | NEB Luna vs Taqpath | UDI | 50 | NP in NS real-world samples; Saliva in NS diluted                 | 1 | 94.33±0.93 | 97.94 | 440.62Mbp | 0.19±0.04 | 17,624,768 | 727±15 | 38 | 0.103/0.086 | 112±16 
v9  |  MiSeq | 28.37 | 40 | 20 | NEB Luna | UDI | 50 and 60 | MTS in NS diluted; NP in NS real-world samples                        | 1 | 95.47±0.89 | 98.74 | 346.55Mbp | 0.48±0.08 | 13,862,051 | 563±20 | 38 | 0.095/0.029 | 172±26 
v8  |  MiSeq | 38.61 | 40 | 20 | NEB Luna | UDI | 40 and 50 | MTS in NS diluted                                                     | 1 | 95.06±0.83 | 98.29 | 475.38Mbp | 0.16±0.02 | 19,015,232 | 776±11 | 38 | 0.084/0.043 | 146±18 
v7  |  MiSeq | 37.28 | 40 | 16 | NEB Luna | UDI | 40 |  MTS in TE diluted; Saliva in TE with Protease +/- dilution     | 1 | 92.71±0.99 | 97.26 | 218.65Mbp | 0.67±0.32 | 8,745,975 | 363±31 | 38 | 0.102/0.044 | 190±34   

<style>
table {
    width:100%;
}
</style>

