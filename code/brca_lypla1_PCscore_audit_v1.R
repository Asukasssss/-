# Independent numerical audit; only aggregate audit output is public.
a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==1);root<-a[1]
f<-read.delim(file.path(root,'public','frozen_score_members.tsv'),check.names=FALSE)
genes<-f$gene[tolower(f$used_in_all3_cohorts)=='true']
stopifnot(length(genes)>=5,!any(c('LYPLA1','LYPLA2')%in%genes))
score_error<-c()
for(co in c('Wu2021','Pal2021_reprocessed','CAMP')) {
 e<-read.delim(file.path(root,'private',paste0(co,'_expression.tsv')),row.names=1,check.names=FALSE)
 p<-read.delim(file.path(root,'private',paste0(co,'_scored_expression.tsv')),row.names=1,check.names=FALSE)
 score<-rowMeans(scale(e[,genes,drop=FALSE]))
 err<-max(abs(score-p$PC_score));stopifnot(err<1e-10)
 score_error[co]<-err
}
d<-read.delim(file.path(root,'private','audit_test_inputs.tsv'),check.names=FALSE)
r<-read.delim(file.path(root,'public','results.tsv'),check.names=FALSE)
out<-list()
for(i in seq_len(nrow(r))){
 z<-d[d$test_id==r$test_id[i],];stopifnot(nrow(z)==r$n[i],!anyDuplicated(z$source_unit))
 if(grepl('ER_PC12',r$test_id[i])){
   # QR residuals from a differently implemented base-R linear model.
   z$PC1rank<-rank(z$PC1);z$PC2rank<-rank(z$PC2)
   xx<-resid(lm(rank(x)~ER+PC1rank+PC2rank,data=z))
   yy<-resid(lm(rank(y)~ER+PC1rank+PC2rank,data=z))
   rho<-cor(xx,yy)
 }else{rho<-cor(z$x,z$y,method='spearman')}
 stopifnot(abs(rho-r$effect[i])<1e-10)
 stopifnot(abs((r$permutation_extreme[i]+1)/(r$permutation_n[i]+1)-r$p_value[i])<1e-12)
 out[[i]]<-data.frame(test_id=r$test_id[i],n=nrow(z),R_rho=rho,absolute_difference=abs(rho-r$effect[i]),P_formula_matches=TRUE)
}
stopifnot(max(abs(p.adjust(r$p_value,'BH')-r$q_value))<1e-12)
write.table(do.call(rbind,out),file.path(root,'public','independent_R_audit.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(c('PASS: common-member scores recomputed in all3 cohorts; 8 correlations; BH8; plusone P formula.',
 paste('Maximum score difference:',max(score_error)),paste('R version:',getRversion()),
 'Permutation draws and bootstrap intervals were not independently repeated.'),file.path(root,'public','independent_R_audit.txt'))
cat('PASS: score and correlation audit\n')
