#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3)
root <- normalizePath(args[1]); donor <- args[2]; mode <- args[3]
stopifnot(mode %in% c('primary','reference_split'))
suppressPackageStartupMessages(library(infercnv))
suppressPackageStartupMessages(library(Matrix))
set.seed(20260925)
options(scipen=100)
input <- file.path(root,'private',donor)
output <- file.path(input,mode)
dir.create(output,showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),file.path(output,'sessionInfo.txt'))
counts <- readMM(file.path(input,'counts.mtx'))
counts <- as(counts,'CsparseMatrix')
rownames(counts) <- readLines(file.path(input,'genes.txt'))
colnames(counts) <- readLines(file.path(input,'cells.txt'))
obj <- CreateInfercnvObject(raw_counts_matrix=counts,
  annotations_file=file.path(input,paste0('annotations_',mode,'.tsv')),
  delim='\t',gene_order_file=file.path(root,'source','gene_order.tsv'),
  ref_group_names='REF_T',min_max_counts_per_cell=c(0,Inf))
params <- list(infercnv_obj=obj,cutoff=0.1,out_dir=output,
  min_cells_per_gene=3,window_length=101,cluster_by_groups=TRUE,
  denoise=TRUE,HMM=FALSE,num_threads=2,analysis_mode='samples',
  plot_steps=FALSE,no_plot=TRUE,save_rds=FALSE,save_final_rds=TRUE,
  write_expr_matrix=FALSE,resume_mode=FALSE)
missing <- setdiff(names(params),names(formals(infercnv::run)))
if(length(missing))stop('Unsupported parameters: ',paste(missing,collapse=','))
writeLines(capture.output(str(params[setdiff(names(params),'infercnv_obj')])),file.path(output,'parameters.txt'))
result <- do.call(infercnv::run,params)
# Final denoised expression ratios are used descriptively, never as malignancy calls.
x <- result@expr.data
go <- result@gene_order
stopifnot(identical(rownames(x),rownames(go)),all(is.finite(x)),all(x>0))
write.table(data.frame(cell=colnames(x),mean_absolute_ratio_deviation=colMeans(abs(x-1)),
  rms_ratio_deviation=sqrt(colMeans((x-1)^2))),file.path(output,'cell_scores.tsv'),
  sep='\t',quote=FALSE,row.names=FALSE)
chr <- as.character(go$chr)
bin <- floor(as.numeric(go$start)/1e7)
keys <- paste(chr,bin,sep=':')
levels <- unique(keys)
bins <- do.call(rbind,lapply(levels,function(key){
  ii <- which(keys==key)
  data.frame(bin_id=key,chr=chr[ii[1]],start=bin[ii[1]]*1e7,
    end=(bin[ii[1]]+1)*1e7,n_genes=length(ii))
}))
means <- t(vapply(levels,function(key){
  ii <- which(keys==key)
  if(length(ii)<10)rep(NA_real_,ncol(x)) else colMeans(log2(x[ii,,drop=FALSE]))
},numeric(ncol(x))))
colnames(means) <- colnames(x)
write.table(bins,file.path(output,'bins.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(bin_id=levels,means,check.names=FALSE),
  gzfile(file.path(output,'bin_log2_ratios.tsv.gz')),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
writeLines(c(paste('infercnv',packageVersion('infercnv')),paste('genes',nrow(x)),
  paste('cells',ncol(x)),paste('min',min(x)),paste('max',max(x))),file.path(output,'DONE'))
