args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(args[1]);library(limma);library(jsonlite)
spec<-fromJSON(file.path(root,"analysis_spec.json"));sets<-fromJSON(file.path(root,"source/frozen_program_sets.json"))
old<-spec$previous_public;out<-list();members<-list()
for(site in c("Tumor","Lung"))for(con in c("shAsns1","shAsns2")){
 a<-read.delim(file.path(old,paste0(site,"_",con,"_all_DE.tsv")),check.names=FALSE)
 a<-a[a$gene!="Asns",];stopifnot(!anyDuplicated(a$gene),all(is.finite(a$F)),all(a$F>=0))
 stat<-sign(a$logFC)*sqrt(a$F);names(stat)<-a$gene
 ix<-lapply(sets,function(g)which(a$gene %in% setdiff(g,"Asns")))
 stopifnot(all(lengths(ix)>=10))
 test<-cameraPR(stat,index=ix,use.ranks=TRUE,inter.gene.cor=0.01,sort=FALSE)
 for(name in names(sets)){
  i<-ix[[name]];t<-test[name,]
  out[[length(out)+1]]<-data.frame(cancer="BRCA",cohort=paste0("GSE104966_",site),stage_id="05_FUNCTION",run_id=basename(root),analysis_version="asns117_programs_v1",analysis_type="cameraPR_signed_sqrtQLF_ranks",metabolite_key=NA,metabolite_name=NA,gene=NA,unit="author_biological_replicate",n=4,n_reference=4,effect_type="descriptive_median_log2FC_in_set",effect=median(a$logFC[i]),ci_lower=NA,ci_upper=NA,p_value=t$PValue,q_value=NA,test_family="three_programs_four_contrasts12",family_n_evaluable=12,status="DONE",reason="competitive_rank_test;fixed_inter_gene_cor0.01;post_selection_exploratory",source_id="GSE104966;MSigDB2025.1.Mm",program=name,site=site,construct=con,set_genes_total=length(unique(setdiff(sets[[name]],"Asns"))),set_genes_tested=length(i),background_genes=nrow(a),camera_direction=t$Direction,fraction_logFC_positive=mean(a$logFC[i]>0),fraction_logFC_negative=mean(a$logFC[i]<0))
  missing<-setdiff(setdiff(sets[[name]],"Asns"),a$gene)
  members[[length(members)+1]]<-data.frame(program=name,site=site,construct=con,mouse_gene=a$gene[i],logFC=a$logFC[i],original_q_global4=a$q_global4[i],status="DONE")
  if(length(missing))members[[length(members)+1]]<-data.frame(program=name,site=site,construct=con,mouse_gene=missing,logFC=NA,original_q_global4=NA,status="NOT_EVALUABLE")
 }
}
b<-do.call(rbind,out);b$q_value<-p.adjust(b$p_value,"BH");write.table(b,file.path(root,"public/program_results.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
write.table(do.call(rbind,members),file.path(root,"public/program_gene_detail.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
write_json(list(status="DONE",program_tests=nrow(b),limma=as.character(packageVersion("limma")),inter_gene_correlation="assumed0.01;not_estimated_from_samples",target_Asns_removed_from_sets_and_background=TRUE,old_edgeR_fits_reused=TRUE),file.path(root,"public/program_validation.json"),pretty=TRUE,auto_unbox=TRUE)
print(b[,c("site","construct","program","set_genes_tested","effect","camera_direction","p_value","q_value")])
