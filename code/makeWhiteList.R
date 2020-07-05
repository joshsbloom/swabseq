# short utility to generate whitelist of barcodes from sample sheet for use with kallisto/bustools
args= commandArgs(trailingOnly=T)
filein=args[1]
fileout=args[2]
print(filein)
print(fileout)

ss=read.delim(filein, stringsAsFactors=F, skip=14, sep=',')
ss$mergedIndex=paste0(ss$index, ss$index2)
write.table(ss$mergedIndex, file=fileout, row.names=F, col.names=F, quote=F)


##ssS=ss[seq(1,768,2),]
##ssR=ss[seq(2,768,2),]

#s96index1=split(ss, ss$index2)[[2]]$index
#r96index1=split(ss, ss$index2)[[1]]$index
#
#
#vid=sapply(split(ss, ss$index2), function(x) x$Sample_ID[1])
#s96index2=names(split(ss, ss$index2))[grep('A01-1', vid)]
#r96index2=names(split(ss, ss$index2))[grep('A01-2', vid)]
#
#
#
#impossible1=apply(expand.grid(r96index1, s96index2), 1, paste, collapse='')
#impossible2=apply(expand.grid(s96index1, r96index2), 1, paste, collapse='')
#
#impossibleCombos=c(impossible1,impossible2)
#write.table(impossibleCombos,
#            file='/data/Covid/swabseq/SwabSeqV13/kallisto/impossibleCombos.txt',
#            row.names=F, col.names=F, quote=F)
#
#expand.grid(s96index1, r96index2)

