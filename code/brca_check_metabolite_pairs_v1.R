# Independent signed-rank and BH check; difference table stays on server165.
args <- commandArgs(trailingOnly=TRUE); out <- args[1]
d <- read.delim(file.path(out,'private','differences_for_R.tsv'),check.names=FALSE)
groups <- split(d,interaction(d$metabolite_name,d$family,drop=TRUE))
rows <- lapply(groups,function(z) {
 x <- z$difference
 p <- if(length(x)<8) NA_real_ else if(all(x==0)) 1 else wilcox.test(x,mu=0,alternative='two.sided',exact=FALSE,correct=TRUE)$p.value
 data.frame(metabolite_name=z$metabolite_name[1],test_family=z$family[1],n=length(x),higher=sum(x>0),lower=sum(x<0),equal=sum(x==0),p=p)
})
r <- do.call(rbind,rows);r$q <- NA_real_
for(f in unique(r$test_family)){ ix <- which(r$test_family==f & is.finite(r$p));r$q[ix] <- p.adjust(r$p[ix],method='BH') }
write.table(r,file.path(out,'private','independent_R_checks.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('R checks completed:',nrow(r),'rows\n')
