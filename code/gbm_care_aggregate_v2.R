# Server-only sample/category sufficient statistics, preserving author QC cells.
args <- commandArgs(TRUE); out <- args[1]
shard <- if (length(args)>=2) as.integer(args[2]) else 0L
nshards <- if (length(args)>=3) as.integer(args[3]) else 1L
stopifnot(file.exists(file.path(out,'.running')))
meta <- readRDS(file.path(out,'source/celltype_meta_data_2025_01_08.RDS'))
genes <- read.delim(file.path(out,'source/genes_unique.tsv'),check.names=FALSE)$gene
files <- read.delim(file.path(out,'source/CARE_filelist.txt'),check.names=FALSE,comment.char='')
files <- files[files[[1]]=='File',]
dir.create(file.path(out,'private/sample_sums'),showWarnings=FALSE)
for (filename in files$Name[((seq_len(nrow(files))-1L) %% nshards)==shard]) {
  dest <- file.path(out,'private/sample_sums',paste0(filename,'.tsv'))
  if (file.exists(dest) && file.exists(paste0(dest,'.audit'))) next
  src <- file.path(out,'source',filename)
  expected <- files$Size[files$Name==filename]
  waited <- 0
  while (!file.exists(src) || file.info(src)$size!=expected) {
    Sys.sleep(5); waited <- waited+5
    if (waited>1800) stop('count input not available')
  }
  # GEO gzip wraps an already gzip-compressed RDS. Unwrap one layer privately.
  tmp <- file.path(out,'private',paste0('current_counts_shard',shard,'.RDS'))
  a <- gzfile(src,'rb'); b <- file(tmp,'wb')
  repeat {buf <- readBin(a,'raw',n=1048576); if (!length(buf)) break; writeBin(buf,b)}
  close(a);close(b)
  x <- readRDS(tmp); stopifnot(is.matrix(x),!anyDuplicated(colnames(x)),!anyDuplicated(rownames(x)))
  id <- strsplit(filename,'_',fixed=TRUE)[[1]][2]
  m <- meta[meta$ID==id,]; idx <- match(m$CellID,colnames(x))
  stopifnot(nrow(m)>0,!anyNA(idx),!anyDuplicated(m$CellID))
  # Full measured gene library, not candidate-only sum; author QC filter by barcode.
  lib <- colSums(x[,idx,drop=FALSE]); stopifnot(all(is.finite(lib)),all(lib>0))
  gi <- match(genes,rownames(x)); measured <- !is.na(gi)
  cnt <- matrix(NA_real_,length(genes),length(idx)); cnt[measured,] <- x[gi[measured],idx,drop=FALSE]
  stopifnot(all(cnt[measured,]>=0),all(cnt[measured,]==floor(cnt[measured,])))
  z <- log1p(sweep(cnt,2,lib,'/')*10000)
  res <- do.call(rbind,lapply(unique(as.character(m$CellType)),function(ct) {
    use <- m$CellType==ct
    data.frame(sample_id=id,cell_type=ct,gene=genes,n_cells=sum(use),
      sum_expression=rowSums(z[,use,drop=FALSE]),sum_detected=rowSums(cnt[,use,drop=FALSE]>0),measured=measured)
  }))
  write.table(res,paste0(dest,'.partial'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  stopifnot(file.rename(paste0(dest,'.partial'),dest))
  audit <- data.frame(sample_id=id,source_cells=ncol(x),author_qc_cells=nrow(m),genes=nrow(x),candidate_measured=sum(measured),exact_barcode_matches=sum(!is.na(idx)),min_library=min(lib),median_library=median(lib),max_library=max(lib))
  write.table(audit,paste0(dest,'.audit'),sep='\t',quote=FALSE,row.names=FALSE)
  unlink(tmp);rm(x,cnt,z);gc(verbose=FALSE)
  cat('aggregated',which(files$Name==filename),'of',nrow(files),'\n');flush.console()
}
