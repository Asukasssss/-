args <- commandArgs(TRUE)
e <- new.env()
load(args[1], envir=e)
for (nm in c('MetaData_M4','MetaData_M5','rna_mapping','master_mapping')) {
 x <- e[[nm]]
 write.table(x, file=file.path(args[2],paste0(nm,'.tsv')), sep='\t', row.names=FALSE, quote=TRUE, na='NA')
 cat('EXPORTED_SERVER_PRIVATE',nm,nrow(x),ncol(x),'\n')
}
