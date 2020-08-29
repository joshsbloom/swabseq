library(gdata)

#r=i7
#f=i5

primers=list(
's2_r'=read.xls('/home/jbloom/Dropbox/Manuscripts/swabseq/idt_1536.xlsx', sheet=1, header=F, stringsAsFactors=F),
's2_f'=read.xls('/home/jbloom/Dropbox/Manuscripts/swabseq/idt_1536.xlsx', sheet=2, header=F, stringsAsFactors=F),
'rpp30_r'=read.xls('/home/jbloom/Dropbox/Manuscripts/swabseq/idt_1536.xlsx', sheet=3, header=F, stringsAsFactors=F),
'rpp30_f'=read.xls('/home/jbloom/Dropbox/Manuscripts/swabseq/idt_1536.xlsx', sheet=4, header=F, stringsAsFactors=F))


plates=c(rep(1,384),rep(2,384), rep(3,384), rep(4,384))

sonly=lapply(primers, function(x) gsub('^.*_', '', x[,1]))

sonlyp=lapply(sonly, function(x) split(x, plates))


asgrids=lapply(sonlyp, function(x) lapply(x, function(y) t(matrix(y, 24,16))))

asgrids2=lapply(asgrids, function(y) lapply(y ,function(x) cbind(c('',toupper(letters[1:16])), rbind(seq(1,24), x))))

asgrids3=lapply(asgrids2, function(x) do.call('rbind', x))

for(x in names(asgrids3)){
    write.table(asgrids3[[x]], paste0('/home/jbloom/Dropbox/Manuscripts/swabseq/', x, '.csv'), 
                sep='\t', quote=F, row.names=F, col.names=F)
}


