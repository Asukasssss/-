# Server-only, one author unit per process. CNV candidates are not malignant calls.
args <- commandArgs(trailingOnly=TRUE)
unit_dir <- normalizePath(args[1]); run_dir <- normalizePath(args[2])
setwd(unit_dir)
suppressPackageStartupMessages(library(infercnv))
suppressPackageStartupMessages(library(Matrix))
options(scipen=100)
set.seed(20260925)
meta <- read.delim(file.path(unit_dir,'metadata.tsv'),row.names=1,check.names=FALSE)
x <- as(readMM(file.path(unit_dir,'counts.mtx')),'CsparseMatrix')
rownames(x) <- readLines(file.path(unit_dir,'genes.txt'))
colnames(x) <- readLines(file.path(unit_dir,'cells.txt'))
stopifnot(identical(colnames(x),rownames(meta)),!('SLC6A6' %in% rownames(x)))
ref_sets <- list(mixed=c('T_ref','B_ref'),T_only='T_ref',B_only='B_ref')
scores <- list()
for (v in names(ref_sets)) {
  out <- file.path(unit_dir,paste0('infercnv_',v))
  obj <- CreateInfercnvObject(raw_counts_matrix=x,
    annotations_file=file.path(unit_dir,'annotations.tsv'),delim='\t',
    gene_order_file=file.path(run_dir,'private/gene_order_no_chr3.tsv'),
    ref_group_names=ref_sets[[v]],chr_exclude=c('chrX','chrY','chrM','chr3'))
  pars <- list(infercnv_obj=obj,cutoff=0.1,out_dir=out,cluster_by_groups=TRUE,
    denoise=TRUE,HMM=FALSE,num_threads=2,analysis_mode='samples',plot_steps=FALSE,
    no_prelim_plot=TRUE,save_rds=FALSE)
  # Plotting does not affect numeric inference; support versions without no_plot.
  if ('no_plot' %in% names(formals(infercnv::run))) pars$no_plot <- TRUE
  obj <- do.call(infercnv::run,pars)
  a <- obj@expr.data
  stopifnot(!('SLC6A6' %in% rownames(a)))
  order <- obj@gene_order
  stopifnot(all(rownames(a) %in% rownames(order)))
  chrom <- as.character(order[rownames(a),'chr'])
  stopifnot(!any(chrom=='chr3'))
  usable_chr <- names(which(table(chrom)>=50))
  stopifnot(length(usable_chr)>=15)
  # Independent held-out lymphocytes set conservative thresholds.
  held <- intersect(rownames(meta)[grepl('_heldout$',meta$cnv_label)],colnames(a))
  normal_myeloid <- intersect(rownames(meta)[meta$cnv_label=='Myeloid_normal_control'],colnames(a))
  stopifnot(length(held)>=40)
  rms <- sqrt(colMeans((a-1)^2))
  chr_means <- sapply(usable_chr,function(ch) colMeans(a[chrom==ch,,drop=FALSE]-1))
  rownames(chr_means) <- colnames(a)
  q99 <- apply(abs(chr_means[held,,drop=FALSE]),2,quantile,probs=.99,type=8)
  global99 <- max(as.numeric(quantile(rms[held],.99,type=8)),1e-8)
  global95 <- as.numeric(quantile(rms[held],.95,type=8))
  if (length(normal_myeloid)>=20) {
    q99 <- pmax(q99,apply(abs(chr_means[normal_myeloid,,drop=FALSE]),2,quantile,probs=.99,type=8))
    global99 <- max(global99,as.numeric(quantile(rms[normal_myeloid],.99,type=8)))
    global95 <- max(global95,as.numeric(quantile(rms[normal_myeloid],.95,type=8)))
  }
  nchr <- rowSums(sweep(abs(chr_means),2,pmax(q99,1e-8),'>'))
  write.table(data.frame(chr=names(q99),threshold99=q99),file.path(unit_dir,paste0('chr_thresholds_',v,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
  state <- ifelse(rms>global99 & nchr>=2,'CNV_ABNORMAL_CANDIDATE',
    ifelse(rms<=global95 & nchr==0,'REFERENCE_LIKE','UNCERTAIN'))
  z <- data.frame(cell=colnames(a),variant=v,score=rms,n_abnormal_chr=nchr,
    state=state,threshold99=global99,threshold95=global95,check.names=FALSE)
  write.table(z,file.path(unit_dir,paste0('scores_',v,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
  write.table(data.frame(cell=rownames(chr_means),chr_means),file.path(unit_dir,paste0('chromosome_means_',v,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
  scores[[v]] <- z
  saveRDS(obj,file.path(unit_dir,paste0('final_',v,'.rds')))
  rm(obj,a,chr_means);gc()
}
writeLines(capture.output(sessionInfo()),file.path(unit_dir,'sessionInfo.txt'))
writeLines('DONE',file.path(unit_dir,'INFERENCE_DONE'))
