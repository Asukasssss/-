args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1]); .libPaths(c(file.path(root,"Rlib"),.libPaths()))
library(edgeR); library(limma);library(jsonlite)
stopifnot(file.exists(file.path(root,"analysis_spec.json")))
spec<-fromJSON(file.path(root,"analysis_spec.json")); stopifnot(spec$primary_dataset=="GSE104966")
raw<-read.delim(gzfile(file.path(root,"source/GSE104966_counts.txt.gz")),row.names=1,check.names=FALSE)
stopifnot(ncol(raw)==24,!anyDuplicated(rownames(raw)),all(is.finite(as.matrix(raw))),all(raw>=0),all(raw==floor(raw)),all(colSums(raw)>0))
meta<-read.delim(file.path(root,"source/author_column_assignment.tsv"),check.names=FALSE)
stopifnot(nrow(meta)==24,all(meta$column==seq_len(24)))
groups<-factor(rep(c("control","shAsns1","shAsns2"),each=4),levels=c("control","shAsns1","shAsns2"))
res<-list();qc<-list();normmeans<-list();validation<-list();coverage<-list()
for(site in c("Tumor","Lung")){
 ix<-which(meta$site==site);stopifnot(length(ix)==12,all(meta$group[ix]==as.character(groups)))
 counts<-as.matrix(raw[,ix]);design<-model.matrix(~0+groups);colnames(design)<-levels(groups)
 y<-DGEList(counts=counts,group=groups);keep<-filterByExpr(y,design=design,min.count=10,min.total.count=15)
 coverage[[site]]<-data.frame(gene=rownames(raw),site=site,tested=keep,reason=ifelse(keep,"tested","low_count_filter"))
 y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y,method="TMM")
 y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
 contrasts<-makeContrasts(shAsns1-control,shAsns2-control,levels=design)
 for(j in 1:2){
  test<-glmQLFTest(fit,contrast=contrasts[,j]);tt<-topTags(test,n=Inf,sort.by="none")$table
  tt$gene<-rownames(tt);tt$site<-site;tt$construct<-paste0("shAsns",j);tt$author_n_treated<-4;tt$author_n_control<-4
  tt$gene_identity<-"author_mouse_gene_symbol;no_human_ortholog_inference";tt$q_within_contrast<-tt$FDR;tt$FDR<-NULL
  tt$ci_lower<-NA_real_;tt$ci_upper<-NA_real_;tt$ci_reason<-"edgeR_QLF_does_not_supply_corresponding_effect_CI;not_fabricated"
  res[[length(res)+1]]<-tt
 }
 qc[[site]]<-data.frame(site=site,group=as.character(groups),column=meta$column[ix],raw_library_total=colSums(counts),nonzero_genes=colSums(counts>0),TMM_factor=y$samples$norm.factors,filtered_library_size=y$samples$lib.size)
 cp<-cpm(y,log=FALSE)
 nm<-sapply(levels(groups),function(g)rowMeans(cp[,groups==g,drop=FALSE]));normmeans[[site]]<-data.frame(gene=rownames(nm),site=site,nm,check.names=FALSE)
 validation[[site]]<-list(input_genes=nrow(raw),tested_genes=sum(keep),low_count_filtered=sum(!keep),design_rank=qr(design)$rank,residual_df=ncol(y)-ncol(design))
 # No sample-level values exported. MDS diagnostics stay server-side.
 pdf(file.path(root,paste0(site,"_MDS_private.pdf")));plotMDS(y,labels=as.character(groups));dev.off()
}
a<-do.call(rbind,res);a$q_global4<-p.adjust(a$PValue,method="BH")
write_tsv<-function(x,name)write.table(x,file.path(root,"public",name),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
# Four separate files keep each public aggregate below the repository size limit.
for(site in c("Tumor","Lung"))for(cn in c("shAsns1","shAsns2"))write_tsv(a[a$site==site & a$construct==cn,],paste0(site,"_",cn,"_all_DE.tsv"))
write_tsv(a[a$gene %in% c("Asns","Gls","Gls2","Glul"),],"glutamine_node_results.tsv")
write_tsv(do.call(rbind,qc),"library_qc.tsv")
write_tsv(do.call(rbind,coverage),"all_gene_test_coverage.tsv")
means<-do.call(rbind,normmeans);write_tsv(means[means$gene %in% c("Asns","Gls","Gls2","Glul"),],"node_group_mean_TMM_CPM.tsv")
summary<-do.call(rbind,lapply(split(a,list(a$site,a$construct),drop=TRUE),function(x)data.frame(site=x$site[1],construct=x$construct[1],tested=nrow(x),q_global_lt005=sum(x$q_global4<.05),up=sum(x$q_global4<.05 & x$logFC>0),down=sum(x$q_global4<.05 & x$logFC<0))))
write_tsv(summary,"contrast_summary.tsv")
concordance<-list()
for(site in c("Tumor","Lung")){
 x<-a[a$site==site & a$construct=="shAsns1",];z<-a[a$site==site & a$construct=="shAsns2",];stopifnot(identical(x$gene,z$gene))
 both<-x$q_global4<.05 & z$q_global4<.05
 concordance[[site]]<-data.frame(site=site,both_q_lt005=sum(both),both_same_direction=sum(both & x$logFC*z$logFC>0),both_opposite_direction=sum(both & x$logFC*z$logFC<0),logFC_spearman=cor(x$logFC,z$logFC,method="spearman"))
}
write_tsv(do.call(rbind,concordance),"construct_concordance.tsv")
validation$status<-"DONE";validation$edgeR<-as.character(packageVersion("edgeR"));validation$limma<-as.character(packageVersion("limma"));validation$R<-R.version.string;validation$global_family_n<-nrow(a);validation$independence<-"author biological replicates; animal identity/pairing unverified; no inferred pairing";validation$metabolomics_analyzed<-FALSE
write_json(validation,file.path(root,"public/validation.json"),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(root,"public/R_sessionInfo.txt"))
print(summary);print(a[a$gene %in% c("Asns","Gls","Gls2","Glul"),c("gene","site","construct","logFC","PValue","q_global4")])
