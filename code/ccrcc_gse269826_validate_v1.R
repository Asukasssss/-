# Independent base-R check of exported cell arithmetic, paired aggregation and exact null.
args<-commandArgs(TRUE);out<-args[1];library(jsonlite)
d<-read.delim(file.path(out,'private/cells_analysis.tsv'),check.names=FALSE);sel<-fromJSON(file.path(out,'source/ccrcc_gse269826_selection.json'))
d$condition<-ifelse(d$tissue=='ccRCC' & d$annotation %in% sel$tumor_labels,'Tumor',ifelse(d$tissue=='Normal' & d$annotation %in% sel$normal_labels,'Normal','Exclude'))
d<-d[d$condition!='Exclude',];recalc<-log1p(10000*d$counts/d$library);stopifnot(max(abs(recalc-d$expression))<1e-12)
meanx<-aggregate(recalc,list(patient=d$patient,condition=d$condition),mean);names(meanx)[3]<-'value'
n<-aggregate(rep(1,nrow(d)),list(patient=d$patient,condition=d$condition),sum);names(n)[3]<-'n';meanx<-merge(meanx,n);meanx<-meanx[meanx$n>=20,]
t<-meanx[meanx$condition=='Tumor',];nn<-meanx[meanx$condition=='Normal',];z<-merge(t,nn,by='patient',suffixes=c('_tumor','_normal'));delta<-z$value_tumor-z$value_normal
sgn<-as.matrix(expand.grid(rep(list(c(-1,1)),length(delta))));p<-mean(abs(sgn%*%delta/length(delta))>=abs(mean(delta))-1e-12)
r<-read.delim(file.path(out,'public/results.tsv'));rr<-r[r$mode=='primary',];stopifnot(length(delta)==rr$n,abs(mean(delta)-rr$effect)<1e-12,abs(p-rr$p_value)<1e-12)
a<-list(status='PASS',independent_language='R',gene_expression_from_raw_counts=TRUE,paired_n=length(delta),effect=mean(delta),p_exact=p,p_wilcoxon=wilcox.test(delta,exact=TRUE)$p.value,patient_identity='author field; no sample order inference')
write_json(a,file.path(out,'public/independent_numeric_validation.json'),pretty=TRUE,auto_unbox=TRUE,digits=16);print(a)
