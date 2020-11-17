get_thresh<-function(pip, alpha=0.05){
    lfdr = sort(1-pip)
    fdr = cumsum(lfdr)/1:length(lfdr)
    thresh = 1 - lfdr[max(which(fdr<=alpha))]
    return(thresh)
}

truth = read.table("results/sim.coloc.gene.list")$V1


total_sites = length(truth)



dd = read.table("coloc_out/coloc.gene.default_prior.out",head=T)
dd_thresh = get_thresh(dd$PPH4)
cat("coloc with estimated enrichment prior (FDR 5% level)\n\n");
cat("PPH4 threshold: ",dd_thresh, " \n")
dd_rej = dd$Gene[dd$PPH4>= dd_thresh]
rej_count = length(dd_rej)
true_rej = intersect(dd_rej,truth)
false_rej_count = length(dd_rej)-length(true_rej)

cat("FDR = ",round(false_rej_count/length(dd_rej),3),"\n") 
cat("Power = ", round((rej_count - false_rej_count)/total_sites, 3), "\n\n")


dd = read.table("coloc_out/coloc.gene.true_prior.out",head=T)
dd_thresh = get_thresh(dd$PPH4)
cat("coloc with TRUE enrichment prior (FDR 5% level)\n\n");
cat("PPH4 threshold: ",dd_thresh, " \n")
dd_rej = dd$Gene[dd$PPH4>= dd_thresh]
rej_count = length(dd_rej)
true_rej = intersect(dd_rej,truth)
false_rej_count = length(dd_rej)-length(true_rej)

cat("FDR = ",round(false_rej_count/length(dd_rej),3),"\n")
cat("Power = ", round((rej_count - false_rej_count)/total_sites, 3), "\n\n")
