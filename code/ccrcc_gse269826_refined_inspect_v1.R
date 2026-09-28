library(Matrix)
library(jsonlite)
a<-commandArgs(TRUE);out<-a[1];src<-a[2]
x<-readRDS(src);m<-x@meta.data;y<-x@assays$RNA@counts
stopifnot(identical(rownames(m),colnames(y)),!anyDuplicated(rownames(m)),!anyDuplicated(rownames(y)),all(y@x>=0),all(y@x==floor(y@x)))
full<-readRDS(file.path(out,'private/whole_fresh_lypla1_cache.rds'));ix<-match(rownames(m),rownames(full$metadata));stopifnot(!anyNA(ix))
stopifnot(all(as.character(m$patient)==as.character(full$metadata$patient[ix])),all(as.character(m$sample)==as.character(full$metadata$sample[ix])))
genes<-rownames(full$counts);stopifnot(all(genes%in%rownames(y)));delta<-max(abs(as.matrix(y[genes,,drop=FALSE])-as.matrix(full$counts[,ix,drop=FALSE])))
stopifnot(delta==0)
lib<-full$library[ix];lib2<-Matrix::colSums(y)
cache<-list(metadata=m,counts=full$counts[,ix,drop=FALSE],library=lib,full_metadata=full$metadata[ix,,drop=FALSE],source_genes=full$source_genes,validation=list(cells=nrow(m),refined_genes=nrow(y),full_genes=full$source_genes,exact_barcode_join=TRUE,patient_sample_match=TRUE,panel_counts_max_difference=delta,refined_full_library_max_difference=max(abs(lib-lib2))))
saveRDS(cache,file.path(out,'private/refined_lypla1_cache.rds'))
write_json(cache$validation,file.path(out,'public/refined_join_validation.json'),pretty=TRUE,auto_unbox=TRUE)
write.table(m,file.path(out,'private/refined_metadata.tsv'),sep='\t',quote=FALSE,col.names=NA)
cat('REFINED',dim(y),'\n');print(names(m));for(k in names(m))if(grepl('type|annot|chr3p|scrublet',k)&&length(unique(m[[k]]))<80){cat('\nFIELD',k,'\n');print(table(m[[k]],useNA='ifany'))}
