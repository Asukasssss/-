#!/usr/bin/env Rscript
# Hirz et al. Nat Commun 2023, author-provided Slide-seqV2 counts/annotations.
# All bead-level outputs remain on server165. No reannotation or imputation.
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
root <- normalizePath(args[1], mustWork=TRUE)
src <- file.path(root,"source")
out <- file.path(root,"private")
dir.create(out,showWarnings=FALSE)
x <- readRDS(file.path(src,"hirz_counts.rds"))
a <- readRDS(file.path(src,"hirz_annotations.rds"))
xy <- readRDS(file.path(src,"hirz_coordinates.rds"))
stopifnot(!anyDuplicated(rownames(a)), identical(names(x),names(xy)))
res <- list(); qc <- list()
for (s in names(x)) {
  m <- x[[s]]; ids <- colnames(m)
  stopifnot(!anyDuplicated(ids), !anyDuplicated(rownames(m)),
            sum(rownames(m)=="LYPLA1")==1, all(ids %in% rownames(a)),
            all(ids %in% rownames(xy[[s]])), all(m@x>=0), all(m@x==round(m@x)))
  total <- Matrix::colSums(m)
  genes <- Matrix::colSums(m>0)
  g <- as.numeric(m["LYPLA1",])
  mt <- grepl("^MT-",rownames(m))
  mito <- Matrix::colSums(m[mt,,drop=FALSE])/total
  res[[s]] <- data.frame(bead=ids,section=s,a[ids,,drop=FALSE],
    xy[[s]][ids,,drop=FALSE],total_umi=total,n_genes=genes,mt_fraction=mito,
    LYPLA1_umi=g,LYPLA1_log1p_cp10k=log1p(1e4*g/total),row.names=NULL)
  qc[[s]] <- data.frame(section=s,n_beads=length(ids),n_genes_matrix=nrow(m),
    min_umi=min(total),median_umi=median(total),median_genes=median(genes),
    total_LYPLA1_umi=sum(g),positive_beads=sum(g>0),
    annotation_match=all(ids %in% rownames(a)),coordinates_match=all(ids %in% rownames(xy[[s]])))
}
d <- do.call(rbind,res)
stopifnot(!anyDuplicated(d$bead),nrow(d)==nrow(a),all(is.finite(d$LYPLA1_log1p_cp10k)))
con <- gzfile(file.path(out,"bead_values.tsv.gz"),"wt")
write.table(d,con,sep="\t",quote=FALSE,row.names=FALSE);close(con)
write.table(do.call(rbind,qc),file.path(out,"section_qc.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo.txt"))
print(do.call(rbind,qc)); print(table(d$section,d$cell2))
