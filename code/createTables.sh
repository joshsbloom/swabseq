basedir="/data/Covid/swabseq/"

#populate each run directory with XLSX file describing run or CSV sample sheet and then run the following
# arguments are 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v7/"  -b 195116944  
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v8/"  -b 195154960 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v9/"  -b 195184025 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v10/" -b 195218041 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v11/" -b 195336198 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v12/" -b 195340204 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v13/" -b 195360215 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v14/" -b 195435293 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v15/" -b 195481337 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v16/" -b 195535379 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v17/" -b 195853687 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v18/" -b 195886700 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v19/" -b 195934811 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v20/" -b 196002849 
#nextseq run at TCGB core
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v21/" 
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v22/" -b 196030855 -e T
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v23/" -b 196116944 -e T
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v24/" -b 196133998 -e T
#basespace upload fail for v25, had to pull bcls off of miseq manually
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v25/" -e T
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v26/" -b 196211019 -e T
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v27/" -b 196203090 -e T
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v28/" -b 196244058 -e T
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v29/" -b 196340151 -e T
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v30/" -b 196373205 -e T
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v31/" -b 196423261 -e T
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v32/" -b 196445335 -e T
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v33/" -b 196530346 -e T
Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${basedir}runs/v34/" -b 196554422 -e T
