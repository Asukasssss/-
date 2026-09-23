args=commandArgs(TRUE);R=normalizePath(args[1]);P=file.path(R,'public');D=file.path(R,'private')
.libPaths(c(file.path(R,'Rlib'),.libPaths()))
suppressPackageStartupMessages(library(WGCNA));suppressPackageStartupMessages(library(jsonlite))
options(stringsAsFactors=FALSE);allowWGCNAThreads(nThreads=4);set.seed(20260923)
saveTab=function(x,n)write.table(x,file.path(P,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
readExpr=function(n)as.data.frame(read.delim(file.path(D,paste0(n,'.tsv')),row.names=1,check.names=FALSE))
X=readExpr('Wu2021_expr');meta=read.delim(file.path(D,'Wu2021_metacells.tsv'))
stopifnot(identical(rownames(X),meta$metacell))
inputGenes=colnames(X);gsg=goodSamplesGenes(X,verbose=0);stopifnot(all(gsg$goodSamples));X=X[,gsg$goodGenes,drop=FALSE]
saveTab(data.frame(gene=inputGenes,pass=gsg$goodGenes),'network_quality_genes.tsv')
powers=c(1:10,seq(12,30,2));sf=pickSoftThreshold(X,powerVector=powers,networkType='signed',corFnc='bicor',corOptions=list(use='p',maxPOutliers=.05),verbose=2)
fit=sf$fitIndices;saveTab(fit,'soft_power.tsv')
ok=which(fit$SFT.R.sq>=.8 & fit$slope<0 & fit$mean.k.>=5)
fallback=length(ok)==0
if(!fallback) power=fit$Power[min(ok)] else {ok=which(fit$mean.k.>=5 & fit$slope<0);stopifnot(length(ok)>0);power=fit$Power[ok[which.max(fit$SFT.R.sq[ok])]]}
cat('POWER',power,'FALLBACK',fallback,'\n');flush.console()
net=blockwiseModules(X,power=power,networkType='signed',TOMType='signed',corType='bicor',maxPOutliers=.05,pearsonFallback='individual',minModuleSize=30,mergeCutHeight=.25,deepSplit=2,maxBlockSize=ncol(X)+1,numericLabels=FALSE,pamRespectsDendro=FALSE,saveTOMs=FALSE,randomSeed=20260923,nThreads=4,verbose=2)
cols=net$colors;names(cols)=colnames(X);ME=moduleEigengenes(X,colors=cols,excludeGrey=TRUE)$eigengenes
k=bicor(X,ME,use='p',maxPOutliers=.05,pearsonFallback='individual')
members=data.frame(gene=colnames(X),module=unname(cols),kME=NA_real_,rank_in_module=NA_integer_)
for(m in setdiff(unique(cols),'grey')) {ix=which(cols==m);members$kME[ix]=k[ix,paste0('ME',m)];members$rank_in_module[ix]=rank(-members$kME[ix],ties.method='min')}
saveTab(members,'gene_modules.tsv');saveRDS(list(net=net,ME=ME,expression=X),file.path(D,'network.rds'))
mods=as.data.frame(table(cols));colnames(mods)=c('module','n_genes');mods$module=as.character(mods$module)
ly=members[members$gene=='LYPLA1',];saveTab(ly,'LYPLA1_membership.tsv')
target=if(nrow(ly))ly$module else 'FILTERED';cat('LYPLA1_MODULE',target,'\n');flush.console()
if(!target%in%c('grey','FILTERED')) {
 ix=which(cols==target);g=colnames(X)[ix];saveTab(members[ix,][order(-members$kME[ix]),],'LYPLA1_module_members.tsv')
 lc=bicor(X[,'LYPLA1'],X[,g,drop=FALSE],use='p',maxPOutliers=.05)
 saveTab(data.frame(gene=g,bicor_with_LYPLA1=as.numeric(lc),signed_adjacency=((1+as.numeric(lc))/2)^power),'LYPLA1_module_edges_all.tsv')
 # Fixed gene set sensitivity, eigengene excludes LYPLA1 to avoid self inclusion.
 loo=list();for(donor in unique(meta$donor)){
   z=X[meta$donor!=donor,g,drop=FALSE];other=setdiff(g,'LYPLA1');e=moduleEigengenes(z[,other,drop=FALSE],rep('module',length(other)))$eigengenes[,1]
   loo[[length(loo)+1]]=data.frame(kME_excluding_self=bicor(z$LYPLA1,e,maxPOutliers=.05),remaining_donors=length(unique(meta$donor))-1)
 }
 saveTab(do.call(rbind,loo),'LYPLA1_leave_one_donor_out.tsv')
 # Plot-ready donor-equal subtype summaries; no metacell-level significance tests.
 donorME=aggregate(ME,by=list(donor=meta$donor),FUN=mean);dm=unique(meta[,c('donor','subtype','treatment')]);donorME=merge(donorME,dm,by='donor')
 write.table(donorME,file.path(D,'donor_eigengenes.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 sr=list();for(st in unique(donorME$subtype))for(m in names(ME)){
    v=donorME[donorME$subtype==st,m];sr[[length(sr)+1]]=data.frame(subtype=st,module=sub('^ME','',m),n_donors=length(v),mean_eigengene=mean(v),sd_eigengene=sd(v))
 };saveTab(do.call(rbind,sr),'module_subtype_descriptive.tsv')
}
# External fixed modules: Pal kME is descriptive, not independent-cell significance.
cat('WU_NETWORK_DONE_WAITING_PREPARATION\n');flush.console()
while(!file.exists(file.path(P,'preparation_validation.json')))Sys.sleep(10)
Y=readExpr('Pal2021_reprocessed_expr');common=intersect(colnames(X),colnames(Y));Y=Y[,common,drop=FALSE]
ey=moduleEigengenes(Y,colors=cols[common],excludeGrey=TRUE)$eigengenes;ky=bicor(Y,ey,use='p',maxPOutliers=.05,pearsonFallback='individual')
members$Pal_kME=NA_real_;for(m in setdiff(unique(cols),'grey')){g=intersect(names(cols)[cols==m],common);members$Pal_kME[match(g,members$gene)]=ky[g,paste0('ME',m)]}
saveTab(members,'gene_modules_external_kME.tsv')
for(m in mods$module){ix=which(cols==m);mods$median_Wu_kME[mods$module==m]=median(members$kME[ix],na.rm=TRUE);mods$median_Pal_kME[mods$module==m]=median(members$Pal_kME[ix],na.rm=TRUE)}
saveTab(mods,'module_summary.tsv')
write_json(list(WGCNA=as.character(packageVersion('WGCNA')),R=R.version.string,power=power,power_fallback=fallback,network_genes=ncol(X),Wu_metacells=nrow(X),Pal_metacells=nrow(Y),modules_excluding_grey=sum(mods$module!='grey'),LYPLA1_module=target),file.path(P,'network_summary.json'),auto_unbox=TRUE,pretty=TRUE)
capture.output(sessionInfo(),file=file.path(P,'R_sessionInfo.txt'))
saveRDS(list(colors=cols,X=X,Y=Y,members=members),file.path(D,'preservation_input.rds'))
pdf(file.path(P,'dendrogram.pdf'),width=16,height=6);plotDendroAndColors(net$dendrograms[[1]],cols[net$blockGenes[[1]]],'Modules',dendroLabels=FALSE,hang=.03,addGuide=TRUE,guideHang=.05);dev.off()
cat('NETWORK_DONE\n')
