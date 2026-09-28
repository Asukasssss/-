# Expression-only object has no per-nucleus malignant/CNV annotation.
root <- commandArgs(TRUE)[1]
x <- readRDS(file.path(root,'source/SS2_expression.rds'))
stopifnot(is.data.frame(x), names(x)[1]=='X1gene_id')
out <- data.frame(object_class=class(x), genes=nrow(x), nuclei=ncol(x)-1,
  object_annotation='gene_id column plus expression columns;no Seurat metadata slots',
  status='NOT_EVALUABLE',reason='need verified cell-type/CNV labels to compare malignant and nonmalignant cells')
write.table(out,file.path(root,'public/SS2_object_inspection.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
