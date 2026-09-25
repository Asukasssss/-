# KEGG GSEA on frozen DE statistics; no cell/patient matrix required here.
args=commandArgs(TRUE);out=args[1];pub=file.path(out,'public')
stopifnot(dir.exists(file.path(out,'.running')))
suppressPackageStartupMessages({library(fgsea);library(jsonlite)})
spec=fromJSON(file.path(pub,'analysis_spec.json'))
d=read.delim(file.path(pub,'gene_statistics.tsv'),stringsAsFactors=FALSE)
gmt=strsplit(readLines(file.path(out,'kegg_tested.gmt')),'\t',fixed=TRUE)
pathways=lapply(gmt,function(z)z[-c(1,2)]);names(pathways)=vapply(gmt,function(z)z[1],character(1))
titles=setNames(vapply(gmt,function(z)z[2],character(1)),names(pathways))
negative=d$QL_F<0
stopifnot(all(d$QL_F[negative]>-1e-6),all(d$p_value[negative]==1),nrow(d)==11058,!anyDuplicated(d$gene))
rank=setNames(sign(d$log2FC)*sqrt(pmax(0,d$QL_F)),d$gene);rank=rank[order(-rank,names(rank))]
set.seed(spec$GSEA_parameters$seed)
g=as.data.frame(fgseaMultilevel(pathways=pathways,stats=rank,minSize=15,maxSize=500,eps=0,sampleSize=101,nPermSimple=10000,nproc=1))
# fgsea internally computes padj by design; discard it, never select by it.
g$padj=NULL;g$pathway_id=g$pathway;g$pathway=unname(titles[g$pathway_id])
g$leadingEdge=vapply(g$leadingEdge,function(x)paste(x,collapse=';'),character(1))
g$selected=is.finite(g$pval)&g$pval<.05
g$status=ifelse(is.finite(g$pval),'DONE','NOT_EVALUABLE');g$reason=ifelse(is.finite(g$pval),'NOMINAL_P_ONLY','FGSEA_P_UNAVAILABLE')
g=g[order(g$pval,g$pathway_id,na.last=TRUE),]
stopifnot(nrow(g)==length(pathways),all(is.finite(g$NES)))
write.table(g,file.path(pub,'kegg_GSEA.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
capture.output(sessionInfo(),file=file.path(pub,'software_versions.txt'))
v=fromJSON(file.path(pub,'preparation_validation.json'));v$status='PASS';v$GSEA_P05=sum(g$selected);v$GSEA_high_P05=sum(g$selected&g$NES>0);v$GSEA_low_P05=sum(g$selected&g$NES<0);v$GSEA_not_evaluable=sum(!is.finite(g$pval));v$rank_zero_clamped=sum(negative);v$no_new_DE=TRUE;v$FDR_selection=FALSE;v$rank_tied_entries=sum(duplicated(rank));v$software=list(R=R.version.string,fgsea=as.character(packageVersion('fgsea')))
write_json(v,file.path(pub,'validation.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines('DONE',file.path(out,'DONE'));unlink(file.path(out,'.running'),recursive=TRUE)
print(v)
