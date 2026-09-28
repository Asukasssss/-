# Independent R recomputation of matched-patient effects, Wilcoxon P, and BH q.
args <- commandArgs(TRUE)
root <- args[1]
d <- read.delim(file.path(root, 'private/combined_patient_profiles.tsv'), check.names=FALSE)
r <- read.delim(file.path(root, 'public/LYPLA1_paired_results.tsv'), check.names=FALSE)
maxeffect <- 0; maxp <- 0
rp <- rep(NA_real_, nrow(r))
for (i in seq_len(nrow(r))) {
  z <- d[d$cohort==r$cohort[i] & d$partition==r$partition[i] & d$n_cells>=r$min_cells[i], ]
  a <- z[z$category=='Malignant', c('patient','mean_expression')]
  b <- z[z$category==r$comparator[i], c('patient','mean_expression')]
  ab <- merge(a,b,by='patient',suffixes=c('_a','_b'))
  stopifnot(nrow(ab)==r$n[i])
  delta <- ab$mean_expression_a-ab$mean_expression_b
  if (length(delta)) maxeffect <- max(maxeffect, abs(mean(delta)-r$effect[i]))
  if (r$status[i]=='DONE') {
    rp[i] <- if(all(delta==0)) 1 else wilcox.test(delta, alternative='two.sided',
      exact=length(delta)<=50 && !any(delta==0), correct=FALSE)$p.value
    maxp <- max(maxp,abs(rp[i]-r$p_value[i]))
  } else stopifnot(is.na(r$p_value[i]), is.na(r$q_value[i]))
}
maxq <- 0
for (f in unique(r$test_family)) {
  k <- which(r$test_family==f & r$status=='DONE')
  if(length(k)) maxq <- max(maxq,max(abs(p.adjust(rp[k],method='BH')-r$q_value[k])))
}
stopifnot(maxeffect<1e-10,maxp<1e-10,maxq<1e-10)
write.table(data.frame(check=c('paired_effect','Wilcoxon_P','BH_q'), max_absolute_error=c(maxeffect,maxp,maxq),
  status='PASS',R_version=R.version.string),file.path(root,'public/independent_R_validation.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(data.frame(maxeffect,maxp,maxq))
