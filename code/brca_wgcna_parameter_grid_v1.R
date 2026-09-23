args=commandArgs(TRUE);R=normalizePath(args[1]);B=normalizePath(args[2]);power=as.integer(args[3]);P=file.path(R,'public');D=file.path(R,'private',paste0('power',power));dir.create(D,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(file.path(B,'Rlib'),.libPaths()));suppressPackageStartupMessages(library(WGCNA));suppressPackageStartupMessages(library(jsonlite));options(stringsAsFactors=FALSE);allowWGCNAThreads(4)
v=readRDS(file.path(B,'private','preservation_input.rds'));X=v$X;Y=v$Y;basecols=v$colors
grid=read.delim(file.path(P,'parameter_grid.tsv'));grid=grid[grid$power==power,]
genes=strsplit('ABHD12 CEPT1 CHKA CHKB CHPT1 ENPP2 ETNK1 ETNK2 GDPD5 GPCPD1 LPCAT1 LPCAT2 LPCAT3 LPCAT4 LYPLA1 LYPLA2 PCYT1A PCYT1B PCYT2 SLC22A2',' ')[[1]]
saveTab=function(x,n)write.table(x,file.path(P,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
for(j in seq_len(nrow(grid))){
 cfg=grid[j,];id=cfg$config;set.seed(20260923);cat('START',id,'\n');flush.console()
 tomBase=file.path(D,'TOM');tom=file.path(D,'TOM-block.1.RData');load=file.exists(tom)
 net=blockwiseModules(X,power=cfg$power,networkType='signed',TOMType='signed',corType='bicor',maxPOutliers=.05,pearsonFallback='individual',minModuleSize=cfg$minModuleSize,mergeCutHeight=cfg$mergeCutHeight,deepSplit=cfg$deepSplit,maxBlockSize=ncol(X)+1,numericLabels=FALSE,pamRespectsDendro=FALSE,saveTOMs=TRUE,saveTOMFileBase=tomBase,loadTOM=load,TOMFiles=if(load)tom else NULL,randomSeed=20260923,nThreads=4,verbose=2)
 cols=net$colors;names(cols)=colnames(X);me=moduleEigengenes(X,cols,excludeGrey=TRUE)$eigengenes;k=bicor(X,me,use='p',maxPOutliers=.05,pearsonFallback='individual')
 common=intersect(names(cols),colnames(Y));ey=moduleEigengenes(Y[,common,drop=FALSE],cols[common],excludeGrey=TRUE)$eigengenes;ky=bicor(Y[,common,drop=FALSE],ey,use='p',maxPOutliers=.05,pearsonFallback='individual')
 out=data.frame(config=id,gene=colnames(X),module=unname(cols),module_size=NA_integer_,kME=NA_real_,rank=NA_real_,rank_fraction=NA_real_,Pal_kME=NA_real_,Pal_rank_fraction=NA_real_)
 for(m in setdiff(unique(cols),'grey')){
  ix=which(cols==m);out$module_size[ix]=length(ix);out$kME[ix]=k[ix,paste0('ME',m)];out$rank[ix]=rank(-out$kME[ix],ties.method='min');out$rank_fraction[ix]=out$rank[ix]/length(ix)
  cg=intersect(names(cols)[ix],common);yi=match(cg,out$gene);out$Pal_kME[yi]=ky[cg,paste0('ME',m)];out$Pal_rank_fraction[yi]=rank(-out$Pal_kME[yi],ties.method='min')/length(cg)
 }
 saveTab(out,paste0(id,'_all_genes.tsv'));saveRDS(net,file.path(D,paste0(id,'_network.rds')))
 q=merge(data.frame(gene=genes),out,by='gene',all.x=TRUE,sort=FALSE);q$config=id;q$status=ifelse(is.na(q$module),'NOT_EVALUABLE',ifelse(q$module=='grey','NOT_EVALUABLE','DONE'));q$reason=ifelse(is.na(q$module),'not in frozen 6000 gene input',ifelse(q$module=='grey','unassigned grey',''))
 q$baseline_module_jaccard=NA_real_;for(i in seq_len(nrow(q))){g=q$gene[i];if(q$status[i]=='DONE'&&g%in%names(basecols)&&basecols[g]!='grey'){a=names(cols)[cols==cols[g]];b=names(basecols)[basecols==basecols[g]];q$baseline_module_jaccard[i]=length(intersect(a,b))/length(union(a,b))}}
 saveTab(q,paste0(id,'_topic20.tsv'))
 saveTab(data.frame(config=id,power=cfg$power,deepSplit=cfg$deepSplit,mergeCutHeight=cfg$mergeCutHeight,minModuleSize=cfg$minModuleSize,modules=sum(unique(cols)!='grey'),grey_genes=sum(cols=='grey'),input_genes=ncol(X),baseline_exact_labels=if(id=='baseline')identical(unname(cols),unname(basecols))else NA),paste0(id,'_summary.tsv'))
 cat('DONE',id,'\n');flush.console()
}
capture.output(sessionInfo(),file=file.path(P,paste0('R_sessionInfo_power',power,'.txt')))
