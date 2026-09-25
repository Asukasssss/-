# Independent base-R checks; private patient differences stay on server165.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
root <- args[1]
d <- read.delim(file.path(root,'private','paired_differences.tsv'),check.names=FALSE)
r <- read.delim(file.path(root,'public','results.tsv'),check.names=FALSE)
c <- read.delim(file.path(root,'public','joint_directions.tsv'),check.names=FALSE)
answer <- list()
for (i in seq_len(nrow(r))) {
  x <- d[d$metabolite_name==r$metabolite_name[i],]
  if (r$analysis_type[i]=='author_available') x <- x[tolower(x$both_author_available)=='true',]
  stopifnot(nrow(x)==r$n[i], !anyDuplicated(x$case_row))
  rho <- cor(x$RNA_delta,x$metabolite_delta,method='spearman')
  stopifnot(abs(rho-r$effect[i])<1e-12)
  cs <- c[c$metabolite_name==r$metabolite_name[i] & c$analysis_type==r$analysis_type[i],]
  state <- function(v) ifelse(v>0,'UP',ifelse(v<0,'DOWN','EQUAL'))
  for (j in seq_len(nrow(cs))) {
    actual <- sum(state(x$RNA_delta)==cs$RNA_direction[j] & state(x$metabolite_delta)==cs$metabolite_direction[j])
    stopifnot(actual==cs$count[j])
  }
  stopifnot(abs((r$permutation_extreme[i]+1)/(r$permutation_n[i]+1)-r$p_value[i])<1e-14)
  answer[[i]] <- data.frame(metabolite_name=r$metabolite_name[i],analysis_type=r$analysis_type[i],
                           n=nrow(x),R_rho=rho,absolute_rho_difference=abs(rho-r$effect[i]),
                           all9_direction_counts_match=TRUE,permutation_P_formula_matches=TRUE)
}
stopifnot(max(abs(p.adjust(r$p_value,'BH')-r$q_value))<1e-12)
write.table(do.call(rbind,answer),file.path(root,'public','independent_R_audit.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(c('PASS: rho recomputed in R for four models; all9 direction counts per model; BH4; plus-one P formula.',
             paste('R version:',getRversion()),
             'Permutation draws and bootstrap intervals were not independently repeated.'),file.path(root,'public','independent_R_audit.txt'))
cat('PASS: independent R audit\n')
