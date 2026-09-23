args=commandArgs(TRUE);R=normalizePath(args[1]);P=file.path(R,'public');D=file.path(R,'private');.libPaths(c(file.path(R,'Rlib'),.libPaths()))
suppressPackageStartupMessages(library(WGCNA));suppressPackageStartupMessages(library(org.Hs.eg.db));suppressPackageStartupMessages(library(GO.db));suppressPackageStartupMessages(library(jsonlite))
options(stringsAsFactors=FALSE);allowWGCNAThreads(nThreads=4)
saveTab=function(x,n)write.table(x,file.path(P,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
v=readRDS(file.path(D,'preservation_input.rds'));genes=names(v$colors)
ann=AnnotationDbi::select(org.Hs.eg.db,keys=intersect(genes,AnnotationDbi::keys(org.Hs.eg.db,keytype='SYMBOL')),keytype='SYMBOL',columns=c('GOALL','ONTOLOGYALL'))
ann=unique(ann[!is.na(ann$GOALL)&ann$ONTOLOGYALL=='BP',c('SYMBOL','GOALL')]);bg=unique(ann$SYMBOL);sets=split(ann$SYMBOL,ann$GOALL);sets=lapply(sets,unique);sets=sets[lengths(sets)>=10&lengths(sets)<=500]
terms=AnnotationDbi::select(GO.db,keys=names(sets),columns='TERM',keytype='GOID');termname=setNames(terms$TERM,terms$GOID)
res=list();mods=setdiff(unique(v$colors),'grey')
for(m in mods){g=intersect(genes[v$colors==m],bg);n=length(g);if(n<1)next
 for(id in names(sets)){hit=intersect(g,sets[[id]]);k=length(hit);K=length(sets[[id]]);res[[length(res)+1]]=data.frame(module=m,GO=id,description=termname[[id]],overlap=k,module_annotated_genes=n,term_universe_genes=K,universe_genes=length(bg),p=phyper(k-1,K,length(bg)-K,n,lower.tail=FALSE),genes=paste(hit,collapse=';'))}
}
en=do.call(rbind,res);en$q_all_module_terms=p.adjust(en$p,'BH');saveTab(en,'GO_BP_all_modules.tsv')
ly=v$members[v$members$gene=='LYPLA1',];if(nrow(ly)&&ly$module!='grey')saveTab(en[en$module==ly$module,][order(en$q_all_module_terms[en$module==ly$module]),],'LYPLA1_module_GO_BP.tsv')
write_json(list(database_version=as.character(packageVersion('org.Hs.eg.db')),GO_db_version=as.character(packageVersion('GO.db')),universe=length(bg),terms=length(sets),module_term_tests=nrow(en),family='all tested non-grey module x GO BP term pairs,including zero overlap',min_term=10,max_term=500),file.path(P,'enrichment_spec.json'),auto_unbox=TRUE,pretty=TRUE)
cat('ENRICHMENT_DONE',nrow(en),'\n');flush.console()
# Formal preservation evaluated on donor aggregates, NOT metacells as replicates.
readD=function(co)as.data.frame(read.delim(file.path(D,paste0(co,'_donor_expr.tsv')),row.names=1,check.names=FALSE))
a=readD('Wu2021');b=readD('Pal2021_reprocessed');common=intersect(genes,intersect(colnames(a),colnames(b)));common=common[apply(a[,common],2,sd)>0&apply(b[,common],2,sd)>0]
a=a[,common,drop=FALSE];b=b[,common,drop=FALSE];cols=v$colors[common]
mp=modulePreservation(list(Wu=list(data=a),Pal=list(data=b)),list(Wu=cols),referenceNetworks=1,nPermutations=100,networkType='signed',corFnc='bicor',corOptions="use='p',maxPOutliers=0.05,pearsonFallback='individual'",randomSeed=20260923,maxModuleSize=1000,maxGoldModuleSize=1000,savePermutedStatistics=FALSE,verbose=2,parallelCalculation=FALSE)
saveRDS(mp,file.path(D,'module_preservation.rds'))
z=mp$preservation$Z$ref.Wu$inColumnsAlsoPresentIn.Pal;o=mp$preservation$observed$ref.Wu$inColumnsAlsoPresentIn.Pal
saveTab(data.frame(module=rownames(z),z,check.names=FALSE),'preservation_Z.tsv');saveTab(data.frame(module=rownames(o),o,check.names=FALSE),'preservation_observed.tsv')
write_json(list(reference_donors=nrow(a),external_donors=nrow(b),common_genes=length(common),permutations=100,max_module_genes=1000,statistical_unit='donor aggregates',purpose='exploratory preservation of metacell-discovered modules across donor-level expression',not_proof='not causal or therapeutic validation'),file.path(P,'preservation_spec.json'),auto_unbox=TRUE,pretty=TRUE)
cat('FOLLOWUP_DONE\n')
