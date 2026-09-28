#!/usr/bin/env Rscript
# Exact barcode join to author's annotation; matrices and cell outputs stay server165.
suppressPackageStartupMessages(library(Matrix))
args=commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
r=normalizePath(args[1]);a=readRDS(file.path(r,'source/hirz_scrna_annotations.rds'));x=readRDS(file.path(r,'source/hirz_scrna_counts.rds'))
stopifnot(!anyDuplicated(rownames(a)))
out=file.path(r,'private/hirz');dir.create(out,showWarnings=FALSE)
genes=Reduce(intersect,lapply(x,rownames));stopifnot(sum(genes=='LYPLA1')==1,!anyDuplicated(genes))
ms=list();os=list();qc=list()
for(s in names(x)){
 m=x[[s]];ids=colnames(m);stopifnot(!anyDuplicated(ids),all(ids%in%rownames(a)))
 cancer=grepl('-T-',s,fixed=TRUE)
 keep=ids[a[ids,'cell1']=='Tumor' & cancer]
 qc[[s]]=data.frame(sample=s,n_all=length(ids),n_author_tumor=sum(a[ids,'cell1']=='Tumor'),n_selected=length(keep))
 if(!length(keep))next
 stopifnot(all(m@x>=0),all(m@x==round(m@x)))
 total=Matrix::colSums(m[,keep,drop=FALSE]);mt=Matrix::colSums(m[grepl('^MT-',rownames(m)),keep,drop=FALSE])/total
 os[[s]]=data.frame(cell=keep,sample=s,author_label=a[keep,'cell1'],library=total,percent_mito=100*mt,UMAP_1=a[keep,'UMAP_1'],UMAP_2=a[keep,'UMAP_2'])
 ms[[s]]=m[genes,keep,drop=FALSE]
}
z=do.call(cbind,ms);o=do.call(rbind,os);stopifnot(identical(colnames(z),o$cell),!anyDuplicated(o$cell))
Matrix::writeMM(z,file.path(out,'counts.mtx'));writeLines(genes,file.path(out,'genes.txt'))
write.table(o,file.path(out,'cells.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,qc),file.path(out,'sample_qc.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
print(do.call(rbind,qc));cat('selected malignant cells',ncol(z),'genes',nrow(z),'\n')
