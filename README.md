# UCLA Amplicon Quantification Software for analyzing amplicon data from SwabSeq Sars-CoV-2 testing
___
see [countAmplicons.R](code/countAmplicons.R) for script to turn BCLs and XLSX description of experiment into tables of counts per amplicon

see [analysis](analysis/) for run-specific analysis scripts 

see [runs](runs/) for BCLs, gzipped fastqs and tables of amplicon counts for each experiment
___

Statistics and notes for each sequencing run: 
___
Run | Instrument | Notes | % PhiX Observed | % PhiX Targeted | Conc. Loaded (pM) | Cluster Density | Cluster PF (%) 
--- | --- | --- | --- | --- | --- | --- | --- 
[v49](analysis/v49/) | NextSeqM     | - | 55 | 40 | 1.8 pM | 129K/mm2 | 82%PF
[v48](analysis/v48/) | MiniSeq Fast | - | 63.5 | 40 | 1.8 pM | 131K/mm2 but very low for lanes 2 and 3 | 84%PF
[v47](analysis/v47/) | NextSeqM     | - | 31 | 40 | 1.7 pM | 155K/mm2 but very low for lanes 2 and 3 | 74%PF
[v46](analysis/v46/) | NextSeqM     | - | 13 | 40 | 1.8 pM | 150K/mm2 but varying per lane | 87%PF
[v45](analysis/v45/) | MiniSeq Fast | - | 36 | 40 | 2 pM | 141K/mm2 | 75%PF
[v44](analysis/v44/) | MiniSeq      | - | 31.4 | 40 | 1.8 pM | 147K/mm2 | 90%PF
[v43](analysis/v43/) | NextSeqM     | - | 23.9 | 40 | 1.8 pM | 151K/mm2 | 92%PF
[v42](analysis/v42/) | NextSeqM | - | 31.9 | 40 | 1.75 pM | 113K/mm2 | 90%PF
[v41](analysis/v41/) | NextSeqM | - | 59.6 | 50 | 2 pM |240/mm2 but varying per lane | 70%PF
[v40](analysis/v40/) | NextSeqM | - | 69 | 70 | 1.5  | 93K/mm2 | 94%PF
[v39](analysis/v39/) | NextSeqH | - | 0.4 | 37.5 | 1.5/2.5 pM (phiX/library 50:50) | 126K/mm2 | 42%PF
___

Run | Instrument | Notes | % PhiX Observed | % PhiX Targeted | Conc. Loaded (pM) | RT | Indexing Strategy | PCR cycles | Read | Cluster PF (%) | % ≥Q30 | Yield | Error Rate% | Reads PF | Density  | Tiles | Legacy Phas/Prephas (%) | Intensity 
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- 
[v38](analysis/v38/) | MiSeq | - | 
[v37](analyis/v37/)  | MiSeq  | NS 1:5 | 29.90 | 35 | 25 | Taqpath | - | 50 | 1 | 80.96 ±6.08 |	96.03 |	559.02Mbp | 	0.35±0.32 | 22,360,804 | 	1,135±60 | 	38 |	0.002/0.002 | 	181±20
[v36](analyis/v36/)  | MiSeq  | NS 1:4 and 1:5 | 21.26 | 35 | 25 | Taqpath | - | 50 | 1  |	88.05±3.46 | 	96.10 | 614.70Mbp |	0.21±0.08 |	24,588,136 | 	1,118±29  |	38 |	0.005/0.000 |	149±25 	
[v35](analyis/v35/)  | MiSeq  | ED | 26.89 | 35 | 25 | Taqpath | - | 50 | 1  |	89.50±2.29 | 96.79  |	638.27Mbp |	0.21±0.18 | 25,530,596 | 	1,133±27 | 	38 |	0.000/0.000 |	156±20 	
[v34](analyis/v34/)  | MiSeq  | ED |  31.96 | 35 | 25 | Taqpath | - | 50 | 1 |	90.38±1.97 |	96.92 |	602.46Mbp |	0.26±0.34 |	24,098,420 | 	1,060±29 |	38 |  0.000/0.000 |	170±27 	
[v33](analyis/v33/)  | MiSeq  | - |  33.65 | 35 | 25 | Taqpath | - | 50 | 1 | 	89.05 ±1.85 | 	96.78 |	637.94Mbp |	0.16±0.02 | 25,517,552 |	1,146±15 | 	38 | 0.001/0.000 |	175±26 	
[v32](analyis/v32/)  | MiSeq  | - |   35.56 | 35 | 25 | Taqpath | - | 50 | 1 | 	92.27±1.46 | 	97.46 |	539.25Mbp |	0.23±0.05 |	21,570,046 | 924±14 | 38 | 	0.000 / 0.002 |	195±31
[v31](analyis/v31/)  | MiSeq  | - |  34.78 | 35 | 25 | Taqpath | - | 50 |    1 |	90.15±3.09 | 	96.87 |	564.97Mbp |	0.17±0.04 |	22,598,764 | 996±21 | 	38 |	0.004/0.001 | 	176±11 	
[v30](analyis/v30/)  |  MiSeq | Saliva TBE .5X Tw20; 95C 30 min; Ashe; ED; | 30.38 | 35 | 25 | Taqpath | Semi-Combinatorial | 50 | 1 |	92.63±1.25 | 97.49 |	487.05Mbp |	0.31±0.07 |	19,482,174 |	815±32 |	38  |	0.000/0.000 | 	164 ±22
[v29](analyis/v29/)  |  MiSeq | Saliva TBE .5X Tw20; 95C 30 min; Ashe; ED; | 27.5 | 35 | 25 | Taqpath | Semi-Combinatorial | 50 | 1 | 92.86±1.38 | 97.84 | 	624.03Mbp |	0.20±0.02 |	24,961,348| 	1,057±25  |	38 | 	0.000/0.000 |	198 ±28 	
[v28](analyis/v28/)  |  MiSeq | Saliva TBE .5X Tw20; 95C 30 min; Ashe; ED; | 41 | 35 | 25 | Taqpath | Semi-Combinatorial | 50 |  1 | 	90.11±2.04 	| 97.21 |  	629.71Mbp |  	0.18±0.34 |	25,188,428 | 	1,117±14 |	38 	| 0.000/0.000 |	180±25 	
[v27](analysis/v27/) |  MiSeq | Saliva TE,RNASecure,Qiagen Protease + 95C heating; ATCC after heat; ED; Ashe 2x TBE w/ 0.5% Tw20 7uL  | 18 | 30 | 25 | Taqpath | UDI | 50 | 1 |	94.15±1.28 | 	98.01 	| 556.26Mbp |	0.22±0.04 |  22,250,444 | 	923±23 | 	38 	| 0.009/0.000 | 	208±31 	
[v26](analysis/v26/) |  MiSeq | Saliva TE/RNASecQP; 2X TBE; TBE + Tween ; VTM +95C heat before or after adding ATCC; ED   | 44  | 40 | 25 | Taqpath | Semi-Combinatorial | 50 | 1 |	93.29±1.14 | 97.87 |  550.14Mbp |	0.15±0.02 |	22,005,736 |921±16 | 	38 | 	0.004/0.007 | 	186±29
[v25](analysis/v25/) |  MiSeq | Saliva TBE+/tween; TE,RNASecure,Qiagen Protease; +/- 95 C heating     | 36 | 40 | 25 | Taqpath | UDI | 50 | 
[v24](analysis/v24/) |  MiSeq | Saliva +/- RNA Secure; NP in VTM, NS, or Amies; dilutions in water    | 42 | 40  | 25 | Taqpath | Semi-Combinatorial | 50 |   1  |	92.85±1.17 | 97.79 |  	483.40Mbp 	| 0.14±0.04 |	19,336,026 | 	829±18 | 	38 |	0.014/0.040| 	183 ±28
v23 |                   MiSeq | NP in VTM diluted; eSwab(aimes); Saliva in cells2cdna; NP in NS diluted; increased lysate vol to 10uL; dilutions in TE | 37 | 40 | 26 | Taqpath | UDI | 50 |  1 |	92.69±0.86 | 	97.88 | 	570.91Mbp | 	0.14±0.03 |	22,836,378 | 971±15 | 38 |	0.000/0.000 |	185 ±26
v22 |                   MiSeq | Saliva undiluted; NP into VTM +/- dilution; v14 rerun;                | 15 | 40   | 25 | Taqpath | Semi-Combinatorial | 40 and 50 | 1 | 96.35±0.88 | 98.57 | 347.74Mbp | 0.20±0.31 | 13,909,713 | 560±20 | 38 | 0.135/0.107 | 195±34
v21 |                 NextSeq | NextSeq; mix of v18,v19,v20                                           | 8  | 35  | 1.5 | Taqpath | UDI | 40 | 1 | 95.54±0.45 | 98.48 | 0.17Gbp | 0.36±0.10 | 27,807,589 | 33±1 | 72 | 0.119/0.183 | 12741±809 
v20 |                   MiSeq | EUA (LoD confirmation);  Saliva undiluted                             |  11 | 35  | 24 | Taqpath | UDI | 40 | 1 | 96.09±0.91 | 98.34 | 284.75Mbp | 0.27±0.35 | 11,390,135 | 443±9 | 38 | 0.270/0.020 | 192±31 
v19 |                   MiSeq | EUA (prelimary LoD and Positive and Negative patient samples)         | 14 | 35 | 26 | Taqpath | UDI | 40 |   1 | 85.26±20.73 | 97.28 | 563.71Mbp | 0.50±0.49 | 22,548,216 | 1,052±47 | 38 | 0.011/0.000 | 189±27       
v18 |                   MiSeq | Negative Patient and attempted contrived EUA (high basecalling error rate, diagnostic for amplicon misassignment) | 11 | 30 | 27 | Taqpath | Semi-Combinatorial | 40 | 1 | 69.46±22.02 | 94.37 | 438.81Mbp | 1.25±0.62 | 17,552,240 | 1,042±57 | 38 | 0.008/0.000 | 192±22       
v17 |                   MiSeq | Negative Patient and attempted contrived EUA (but too high viral input for LoD experiment)   | 17 | 30 | 26 | Taqpath | Semi-Combinatorial| 40                      | 1 | 93.83±1.71 | 97.63 | 555.45Mbp | 0.13±0.01 | 22,218,176 | 934±12 | 38 | 0.006/0.036 | 175±21             
v16 |                   MiSeq | v13/v14 and Simulated Patients; first test no RPP30 unindexed         | 17 | 30 | 31.3 | Taqpath | Combinatorial | 40 | 1 | 94.99±1.07 | 98.20 | 598.33Mbp | 0.22±0.34 | 23,933,240 | 973±24 | 38 | 0.032/0.000 | 193±31             
v15 |                   MiSeq | rerun V13 with 384-well plates 1-3; bleach wash                       | 18 | 30 | 40.8 | Taqpath | Combinatorial | 50 |  1 | 75.37±7.07 | 94.56 | 497.81Mbp | 0.83±0.51 | 19,912,204 | 1,103±45 | 38 | 0.027/0.078 | 177±20         
v14 |                   MiSeq | Saliva in NS diluted; Contrived in HEK extracted; no bleach wash      | 17 | 30 | 40.8 | Taqpath | Combinatorial | 50 |  1 | 97.00±0.44 | 99.01 | 463.36Mbp | 0.12±0.01 | 18,534,400 | 742±13 | 38 | 0.138/0.101 | 203±24 
v13 |                   MiSeq | Saliva in NS diluted; Contrived in HEK extracted; Titrate Rxn Vols    | 15 | 30 | 40.8 | Taqpath | Combinatorial | 50 | 1 | 97.16±0.97 | 98.96 | 274.97Mbp | 0.14±0.01 | 10,998,699 | 406±21 | 38 | 0.137/0.101 | 193±30 
v12 |                   MiSeq | MTS in TE diluted and titration of RPP30 unindexed primers            | 44 | 40 | 20 | Taqpath | UDI | 50 |  1 | 93.18±0.96 | 97.35 | 171.26Mbp | 0.35±0.11 | 6,850,514 | 277±24 | 38 | 0.166/0.105 | 147±16   
v11 |                   MiSeq | MTS in TE diluted                                                     | 31 | 40 | 22 | Taqpath | UDI | 50 |  1 | 94.14±1.02 | 98.37 | 423.37Mbp | 0.52±0.13 | 16,934,964 | 703±31 | 38 | 0.090/0.064 | 152±21 
v10 |                   MiSeq | NP in NS real-world samples; Saliva in NS diluted                     | 36 | 40 | 20 | NEB Luna vs Taqpath | UDI | 50 |  1 | 94.33±0.93 | 97.94 | 440.62Mbp | 0.19±0.04 | 17,624,768 | 727±15 | 38 | 0.103/0.086 | 112±16 
v9  |                   MiSeq | MTS in NS diluted; NP in NS real-world samples                        | 28 | 40 | 20 | NEB Luna | UDI | 50 and 60 |  1 | 95.47±0.89 | 98.74 | 346.55Mbp | 0.48±0.08 | 13,862,051 | 563±20 | 38 | 0.095/0.029 | 172±26 
v8  |                   MiSeq | MTS in NS diluted                                                     | 39 | 40 | 20 | NEB Luna | UDI | 40 and 50 |  1 | 95.06±0.83 | 98.29 | 475.38Mbp | 0.16±0.02 | 19,015,232 | 776±11 | 38 | 0.084/0.043 | 146±18 
v7  |                   MiSeq | MTS in TE diluted; Saliva in TE with Protease +/- dilution            | 37 | 40 | 16 | NEB Luna | UDI | 40 |  1 | 92.71±0.99 | 97.26 | 218.65Mbp | 0.67±0.32 | 8,745,975 | 363±31 | 38 | 0.102/0.044 | 190±34   

