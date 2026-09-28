args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
cat('READ_START',as.character(Sys.time()),'\n')
obj <- readRDS(args[1]);run <- args[2]
cat('CLASS',class(obj),'\n');cat('SLOTS',slotNames(obj),'\n')
cat('ASSAYS',names(obj@assays),'\n')
for(n in names(obj@assays)){
 a<-obj@assays[[n]];cat('ASSAY',n,class(a),dim(a),'\n')
 cat('TARGET',grep('^LYPLA1$|^ENSG00000120992$',rownames(a),value=TRUE),'\n')
}
cat('META_COLUMNS',colnames(obj@meta.data),'\n')
for(n in colnames(obj@meta.data)){
 v<-obj@meta.data[[n]]
 if(length(unique(v))<=80){cat('META',n,'\n');print(table(v,useNA='ifany'))}
}
cat('IMAGES',names(obj@images),'\n')
for(n in head(names(obj@images),2)){
 im<-obj@images[[n]];cat('IMAGE',n,class(im),slotNames(im),'\n')
 if('coordinates'%in%slotNames(im))cat('COORDS',dim(im@coordinates),colnames(im@coordinates),'\n')
 if('scale.factors'%in%slotNames(im))print(im@scale.factors)
 if('image'%in%slotNames(im))cat('IMAGE_DIM',dim(im@image),'\n')
}
saveRDS(obj@meta.data,file.path(run,'private','metadata.rds'))
cat('INSPECTION_COMPLETE','\n')
