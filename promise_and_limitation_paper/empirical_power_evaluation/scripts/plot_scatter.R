get_thresh<-function(pip, alpha=0.05){
    lfdr = sort(1-pip)
    fdr = cumsum(lfdr)/1:length(lfdr)
    thresh = 1 - lfdr[max(which(fdr<=alpha))]
    return(thresh)
}


d1 = read.table("results/fastenloc.est_prior.summary", head=T)
d2 = read.table("fastenloc_out/sim.enloc.sig.out",head=T)

eff_eqtl = d1$eQTL_eff 
eff_gwas = d1$GWAS_eff
spip_gwas = d1$GWAS_spip
spip_eqtl = d1$eQTL_spip
rcp = d1$RCP


rcp_cutoff = get_thresh(d2$RCP)
gwas_cutoff = get_thresh(d2$CPIP_qtl)
eqtl_cutoff = get_thresh(d2$CPIP_gwas_marginal)

pdf(file="results/all_class.pdf", height=9.5, width=9.5, bg="white", pointsize=18)

plot(eff_eqtl~eff_gwas, xlim=c(-4,4), ylim=c(-4,4),pch=16,col="gray",xlab = "GWAS effect", ylab="eQTL effect")
points(eff_eqtl[spip_gwas>= gwas_cutoff&spip_eqtl >= eqtl_cutoff] ~ eff_gwas[spip_gwas>= gwas_cutoff&spip_eqtl >= eqtl_cutoff],pch=16, col="cyan")
points(eff_eqtl[rcp>=rcp_cutoff]~ eff_gwas[rcp>=rcp_cutoff], pch=16,col="red")

abline(v=0, lty=2)
abline(h=0, lty=2)

legend("bottomrigh", pch = c(16, 16, 16), col=c("gray", "cyan", "red"), legend=c("Class I False Negative", "Class II False Negative", "Detected"))
dev.off()


