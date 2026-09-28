# Server-only: extract the selected study metadata and LYPLA1, not full atlas export.
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
args=commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2)
out=args[2];dir.create(file.path(out,'private'),recursive=TRUE,showWarnings=FALSE)
message('READ_ATLAS')
obj=readRDS(args[1]);md=slot(obj,'meta.data')
stopifnot('GSE.SRA..Study.' %in% colnames(md));sel=which(as.character(md[['GSE.SRA..Study.']])=='GSE205013');stopifnot(length(sel)>0)
m=md[sel,,drop=FALSE];m$cell_key=rownames(m)
write.table(m,gzfile(file.path(out,'private','atlas_GSE205013_metadata.tsv.gz')),sep='\t',quote=FALSE,row.names=FALSE)
report=c(paste('SELECTED_CELLS',length(sel)),paste('COLUMNS',paste(colnames(m),collapse=' | ')))
for(k in colnames(m)){
 if(k %in% c('Name','orig.ident','cell_key') || grepl('patient|sample|barcode',k,ignore.case=TRUE))next
 if(length(unique(m[[k]]))<=30){report=c(report,paste('COLUMN',k),capture.output(print(table(m[[k]],useNA='ifany'))))}
}
writeLines(report,file.path(out,'public','atlas_schema.txt'))
rna=slot(obj,'assays')[['RNA']];stopifnot(!is.null(rna));counts=GetAssayData(rna,slot='counts')
stopifnot(sum(rownames(counts)=='LYPLA1')==1,identical(colnames(counts),rownames(md)))
ly=as.numeric(counts['LYPLA1',sel]);lib=Matrix::colSums(counts[,sel,drop=FALSE])
stopifnot(all(is.finite(ly)),all(ly>=0),all(ly==floor(ly)),all(lib>0))
d=data.frame(cell_key=rownames(m),LYPLA1_count=ly,total_UMI=lib,LYPLA1_log1p=log1p(10000*ly/lib))
write.table(d,gzfile(file.path(out,'private','atlas_GSE205013_LYPLA1.tsv.gz')),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c('ATLAS_EXTRACTION_DONE','Labels require review before malignancy contrasts; no automatic label promotion.'),file.path(out,'ATLAS_EXTRACTION_DONE.txt'))
message('ATLAS_EXTRACTION_DONE')
