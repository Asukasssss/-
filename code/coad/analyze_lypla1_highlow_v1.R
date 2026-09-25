# Donor-paired full-transcriptome association, not causal gene perturbation.
args=commandArgs(TRUE);out=args[1];lock_commit=args[2]
stopifnot(nchar(lock_commit)==40,dir.exists(file.path(out,'.running')))
.libPaths(c(file.path(out,'Rlib'),.libPaths()))
suppressPackageStartupMessages({library(edgeR);library(fgsea);library(jsonlite)})
pub=file.path(out,'public');spec=fromJSON(file.path(out,'analysis_spec.json'))
counts=as.matrix(read.delim(file.path(out,'private_pseudobulk_counts.tsv'),row.names=1,check.names=FALSE))
s=read.delim(file.path(out,'private_samples.tsv'),stringsAsFactors=FALSE)
stopifnot(identical(colnames(counts),s$sample),all(counts>=0),all(counts==floor(counts)))
s$patient=factor(s$patient);s$group=factor(s$group,levels=c('low','high'))
design=model.matrix(~patient+group,s)
stopifnot(qr(design)$rank==ncol(design),length(levels(s$patient))>=spec$min_donors)
resume=length(args)>2 && args[3]=='resume'
if(!resume){
y=DGEList(counts=counts)
keep=filterByExpr(y,design=design) & rownames(y)!='LYPLA1'
filter_table=data.frame(gene=rownames(y),status=ifelse(rownames(y)=='LYPLA1','EXCLUDED_GROUPING_GENE',ifelse(keep,'TESTED','LOW_EXPRESSION')))
write.table(filter_table,file.path(pub,'gene_testability.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
y=calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM')
fit=glmQLFit(estimateDisp(y,design,robust=TRUE),design,robust=TRUE)
test=glmQLFTest(fit,coef='grouphigh')
res=topTags(test,n=Inf,sort.by='none')$table;res$gene=rownames(res)
res=res[,c('gene','logFC','logCPM','F','PValue','FDR')]
names(res)=c('gene','log2FC','logCPM','QL_F','p_value','q_value')
res$selected=res$q_value<.05 & abs(res$log2FC)>=.5
res$direction=ifelse(res$log2FC>0,'higher_in_LYPLA1_high','lower_in_LYPLA1_high')
# Repeat the same gene family with a prespecified depth sensitivity term.
s$log_mean_umi=log(s$mean_total_umi)
d2=model.matrix(~patient+log_mean_umi+group,s)
depth_ok=qr(d2)$rank==ncol(d2) && nrow(d2)-ncol(d2)>=3
res$depth_log2FC=NA_real_;res$depth_p_value=NA_real_;res$depth_q_value=NA_real_
if(depth_ok){
  f2=glmQLFit(estimateDisp(y,d2,robust=TRUE),d2,robust=TRUE)
  t2=topTags(glmQLFTest(f2,coef='grouphigh'),n=Inf,sort.by='none')$table
  stopifnot(identical(rownames(t2),res$gene))
  res$depth_log2FC=t2$logFC;res$depth_p_value=t2$PValue;res$depth_q_value=t2$FDR
}
res$depth_direction_agrees=sign(res$log2FC)==sign(res$depth_log2FC)
lc=cpm(y,log=TRUE,prior.count=2)
hi=which(s$group=='high');lo=which(s$group=='low');stopifnot(identical(s$patient[hi],s$patient[lo]))
diff=lc[,hi,drop=FALSE]-lc[,lo,drop=FALSE]
res$n_donors=nlevels(s$patient);res$donors_high_greater=rowSums(diff>0);res$donors_high_lower=rowSums(diff<0)
res$status='DONE';res$reason='EXPLORATORY_DONOR_PAIRED_COEXPRESSION_NOT_CAUSAL'
res=res[order(res$p_value,res$gene),]
write.table(res,file.path(pub,'results.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(res[res$selected,],file.path(pub,'selected_genes.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
saveRDS(list(y=y,fit=fit,test=test,samples=s,design=design),file.path(out,'private_models.rds'))
}else{
  res=read.delim(file.path(pub,'results.tsv'),stringsAsFactors=FALSE)
  stopifnot(!anyDuplicated(res$gene),!('LYPLA1' %in% res$gene),all(res$n_donors==nlevels(s$patient)))
  depth_ok=all(is.finite(res$depth_p_value))
}
gmt=strsplit(readLines(file.path(out,'reactome.gmt')),'\t',fixed=TRUE)
stopifnot(all(vapply(gmt,function(z)grepl('^R-HSA-',z[2]),logical(1))))
pathways=lapply(gmt,function(z)intersect(z[-c(1,2)],res$gene));names(pathways)=vapply(gmt,function(z)z[2],character(1))
titles=setNames(vapply(gmt,function(z)z[1],character(1)),names(pathways));stopifnot(!anyDuplicated(names(pathways)))
sizes=lengths(pathways);eligible=sizes>=15 & sizes<=500;pathways=pathways[eligible]
universe=res$gene;M=length(universe)
ora=list()
for(direction in c('higher_in_LYPLA1_high','lower_in_LYPLA1_high')){
  selected=res$gene[res$selected & res$direction==direction];N=length(selected)
  for(id in names(pathways)){
    pg=pathways[[id]];hit=intersect(pg,selected);k=length(hit);K=length(pg)
    ora[[length(ora)+1]]=data.frame(pathway_id=id,pathway=titles[[id]],direction=direction,background_genes=M,
      selected_genes=N,pathway_tested_genes=K,overlap=k,fold_enrichment=if(N>0)k/N/(K/M) else NA_real_,
      p_value=if(N>0)phyper(k-1,K,M-K,N,lower.tail=FALSE) else 1,overlap_genes=paste(sort(hit),collapse=';'),
      status=if(N>0)'DONE' else 'NOT_EVALUABLE',reason=if(N>0)'EXPLORATORY_ORA' else 'NO_SELECTED_GENES')
  }
}
ora=do.call(rbind,ora);ora$q_value=p.adjust(ora$p_value,'BH');ora=ora[order(ora$q_value,ora$p_value),]
write.table(ora,file.path(pub,'reactome_ORA.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
# Numerical QL likelihood rounding may be slightly below zero with P=1.
# Preserve original F/P/q; clamp only these null ranking values to zero.
negative=res$QL_F<0
stopifnot(all(res$QL_F[negative]>-1e-6),all(res$p_value[negative]==1))
rank=setNames(sign(res$log2FC)*sqrt(pmax(0,res$QL_F)),res$gene);rank=rank[order(-rank,names(rank))]
set.seed(spec$seed)
g=as.data.frame(fgseaMultilevel(pathways=pathways,stats=rank,minSize=15,maxSize=500,eps=0,sampleSize=101,nPermSimple=10000,nproc=1))
g$pathway_id=g$pathway;g$pathway=unname(titles[g$pathway_id]);g$leadingEdge=vapply(g$leadingEdge,function(x)paste(x,collapse=';'),character(1))
g=g[order(g$padj,g$pval),];write.table(g,file.path(pub,'reactome_GSEA.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
summary=list(status='PASS',code_lock_commit=lock_commit,patients=nlevels(s$patient),pseudobulks=nrow(s),genes_source=nrow(counts),
  genes_tested=nrow(res),DE_q05=sum(res$q_value<.05),selected_high=sum(res$selected & res$log2FC>0),selected_low=sum(res$selected & res$log2FC<0),
  depth_sensitivity_evaluable=depth_ok,depth_selected=sum(res$depth_q_value<.05 & abs(res$depth_log2FC)>=.5,na.rm=TRUE),
  primary_selected_depth_direction_agrees=sum(res$selected & res$depth_direction_agrees,na.rm=TRUE),
  pathways_tested=length(pathways),ORA_q05=sum(ora$q_value<.05),GSEA_q05=sum(g$padj<.05,na.rm=TRUE),
  full_rank_design=TRUE,grouping_gene_excluded=TRUE,no_independent_cell_test=TRUE,patient_values_private=TRUE,
  numerical_null_F_clamped_for_ranking=sum(negative),differential_values_reused=resume,
  software=list(R=R.version.string,edgeR=as.character(packageVersion('edgeR')),fgsea=as.character(packageVersion('fgsea'))))
write_json(summary,file.path(pub,'validation.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(pub,'software_versions.txt'))
writeLines('DONE',file.path(out,'DONE'));unlink(file.path(out,'.running'),recursive=TRUE)
print(summary)
