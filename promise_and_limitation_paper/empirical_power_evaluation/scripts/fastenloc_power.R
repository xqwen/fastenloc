get_thresh<-function(pip, alpha=0.05){
    lfdr = sort(1-pip)
    fdr = cumsum(lfdr)/1:length(lfdr)
    thresh = 1 - lfdr[max(which(fdr<=alpha))]
    return(thresh)
}

de_summary = read.table("results/fastenloc.est_prior.summary",head=T)
dt_summary = read.table("results/fastenloc.true_prior.summary",head=T)

de_truth = de_summary$sig_cluster

total_sites = length(de_truth)



de = read.table("fastenloc_out/sim.enloc.sig.out",head=T)
de_thresh = get_thresh(de$RCP)
cat("fastENLOC with estimated enrichment prior (FDR 5% level)\n\n");
cat("RCP threshold: ",de_thresh, " \n")
de_rej = de$Signal[de$RCP>= de_thresh]
rej_count = length(de_rej)
true_rej = intersect(de_rej,de_truth)
false_rej_count = length(de_rej)-length(true_rej)

cat("FDR = ",round(false_rej_count/length(de_rej),3),"\n") 
cat("Power = ", round((rej_count - false_rej_count)/total_sites, 3), "\n\n")



dt = read.table("fastenloc_out/sim.true_prior.enloc.sig.out",head=T)
dt_thresh = get_thresh(dt$RCP)
cat("fastENLOC with TRUE enrichment prior (FDR 5% level)\n\n");
cat("RCP threshold: ",dt_thresh, " \n")
dt_rej = dt$Signal[dt$RCP>= dt_thresh]
rej_count = length(dt_rej)
true_rej = intersect(dt_rej,de_truth)
false_rej_count = length(dt_rej)-length(true_rej)

cat("FDR = ",round(false_rej_count/length(dt_rej),3),"\n")
cat("Power = ", round((rej_count - false_rej_count)/total_sites, 3), "\n\n")



