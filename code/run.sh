#/data/Covid/swabseq/SwabSeqV14/kallisto/kout
# index is in /data/Covid/swabseq/index.idx


#bcl2fastq --runfolder-dir ~/Desktop/test/ --output-dir ~/Desktop/test/out --create-fastq-for-index-reads --use-bases-mask=Y26,I10,I10 --processing-threads 64 --no-lane-splitting --sample-sheet /dev/null
python3.7 ../platemap2samp.py SwabSeqv20.xlsx 
bs download run --id 196002849 -o kallisto/
cd kallisto
bcl2fastq --runfolder-dir . --output-dir out/ --create-fastq-for-index-reads  --ignore-missing-bcl --use-bases-mask=Y26,I10,I10 --processing-threads 64 --no-lane-splitting --sample-sheet /dev/null 















cd kallisto
#kallisto index -i ../../index.idx -k 11 ../../reference_amplicons.fa
bcl2fastq --runfolder-dir . --output-dir out/ --create-fastq-for-index-reads  --ignore-missing-bcl --use-bases-mask=Y26,I10,I10 --processing-threads 64 --no-lane-splitting --sample-sheet /dev/null 
Rscript --vanilla /data/Covid/swabseq/makeWhiteList.R ../SampleSheet.csv whitelist.txt
# swabseq10 specified 10bp indices 
kallisto bus -x swabseq10 --i out/Undetermined_S0 --threads 32 --output-dir kout -i ../../index.idx out/Undetermined_S0_I1_001.fastq.gz \
                                                                                              out/Undetermined_S0_I2_001.fastq.gz \
                                                                                              out/Undetermined_S0_R1_001.fastq.gz

bustools sort -t 32 -o kout/sort.bus kout/output.bus
bustools text -p kout/sort.bus > kout/sort.txt
bustools correct -s -d kout/dump.txt -w whitelist.txt -o kout/sort.correct.bus kout/sort.bus
bustools sort -t 32 -o kout/sort.correct.sort.bus kout/sort.correct.bus
bustools text -p kout/sort.correct.sort.bus > kout/data.txt

#bustools correct -s -d kout/dumpH.txt -w ../hopped.txt -o kout/sort.correctH.bus kout/sort.bus
#bustools sort -t 32 -o kout/sort.correctH.sort.bus kout/sort.correctH.bus
#bustools text -p kout/sort.correctH.sort.bus > kout/dataH.txt

#bustools correct -d kout/dump.txt -w whitelistV12_unique.txt -o kout/sort.correctV12.bus kout/sort.bus
#bustools sort -t 32 -o kout/sort.correctV12.sort.bus kout/sort.correctV12.bus
#bustools text -p kout/sort.correctV12.sort.bus > kout/dataV12.txt
#
#bustools correct -d kout/dump.txt -w impossibleCombos.txt -o kout/sort.correctI.bus kout/sort.bus
#bustools sort -t 32 -o kout/sort.correctI.sort.bus kout/sort.correctI.bus
#bustools text -p kout/sort.correctI.sort.bus > kout/dataI.txt


##w=read.delim('/data/Covid/swabseq/SwabSeqV13/kallisto/whitelist.txt', header=F, stringsAsFactors=F)[,1]
#from 
#x=fread('kout/sort.txt', header=F)
#y=x[order(x[,4], decreasing=T),]
#h=y[!(y$V1 %in% w),]

