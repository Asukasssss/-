# Independent R formula check of the fixed paired mean comparison; server only.
args <- commandArgs(trailingOnly=TRUE)
out <- args[1]
root <- '/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716'
map <- read.delim(file.path(out,'private','paired_RNA_map.tsv'),check.names=FALSE)
rna <- read.csv(file.path(root,'data/candidates/camp_primary_tissue_multicancer/gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed/GSE37751.hugene10st.gene_symbol.csv'),row.names=1,check.names=FALSE)
genes <- read.delim(file.path(root,'results/BRCA_117_screen_20260911/tables/BRCA_117_RNA_differences.tsv'))$gene
stopifnot(nrow(map)==45,length(genes)==117)
rows <- lapply(genes,function(g) {
  if(!g %in% rownames(rna)) return(data.frame(gene=g,effect=NA,p=NA,lower=NA,upper=NA))
  x <- as.numeric(rna[g,map$RNAID_tumor]); y <- as.numeric(rna[g,map$RNAID_normal])
  good <- is.finite(x)&is.finite(y)
  if(sum(good)<8 || sd(x[good]-y[good])==0) return(data.frame(gene=g,effect=NA,p=NA,lower=NA,upper=NA))
  z <- t.test(x[good],y[good],paired=TRUE)
  data.frame(gene=g,effect=unname(z$estimate),p=z$p.value,lower=z$conf.int[1],upper=z$conf.int[2])
})
d <- do.call(rbind,rows);d$q <- NA;ok <- is.finite(d$p);d$q[ok] <- p.adjust(d$p[ok],method='BH')
write.table(d,file.path(out,'private','R_paired_check.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('Independent paired R check completed:',sum(ok),'genes\n')
