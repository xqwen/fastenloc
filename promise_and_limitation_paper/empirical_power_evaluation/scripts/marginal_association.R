get_thresh<-function(pip, alpha=0.05){
    lfdr = sort(1-pip)
    fdr = cumsum(lfdr)/1:length(lfdr)
    thresh = 1 - lfdr[max(which(fdr<=alpha))]
    return(thresh)
}

d_truth = read.table("sim_data/sim.truth.rst",head=T)
attach(d_truth)


d_gwas = read.table("summary/sim.gwas_analysis.summary")
gwas_thresh = get_thresh(d_gwas$V2)
gwas_total = length(which(gwas_eff!=0))
gwas_fdr = length(which(d_gwas$V2>=gwas_thresh & d_gwas$V3 == 0))/length(which(d_gwas$V2>=gwas_thresh))
gwas_power =  length(which(d_gwas$V2>=gwas_thresh & d_gwas$V3 == 1))/gwas_total

cat(paste0("GWAS (FDR 5%):   SPIP Cutoff  ", round(gwas_thresh,3),  " Realized FDR ", round(gwas_fdr,3), "   Power ", round(gwas_power,3), " (", length(which(d_gwas$V2>=gwas_thresh & d_gwas$V3 == 1)),"/",gwas_total,")"),'\n')





d_eqtl = read.table("summary/sim.eqtl_analysis.summary")
eqtl_thresh = get_thresh(d_eqtl$V2)
eqtl_total = length(which(eqtl_eff!=0))
eqtl_fdr = length(which(d_eqtl$V2>=eqtl_thresh & d_eqtl$V3 == 0))/length(which(d_eqtl$V2>=eqtl_thresh))
eqtl_power =  length(which(d_eqtl$V2>=eqtl_thresh & d_eqtl$V3 == 1))/eqtl_total

cat(paste0("eQTL (FDR 5%):   SPIP Cutoff ", round(eqtl_thresh,3), " Realized FDR ", round(eqtl_fdr,3), "   Power ", round(eqtl_power,3), " (", length(which(d_eqtl$V2>=eqtl_thresh & d_eqtl$V3 == 1)),"/",eqtl_total,")"), '\n')


d = read.table("results/fastenloc.est_prior.summary", head=T)
cat(paste0("GWAS power in colocalized sites:  ", round(length(d$GWAS_spip[d$GWAS_spip >= gwas_thresh])/length(d$GWAS_spip),3)," (",length(d$GWAS_spip[d$GWAS_spip >= gwas_thresh]),"/",length(d$GWAS_spip),")\n"))
cat(paste0("eQTL power in colocalized sites:  ", round(length(d$eQTL_spip[d$eQTL_spip >= eqtl_thresh])/length(d$eQTL_spip),3)," (",length(d$eQTL_spip[d$eQTL_spip >= eqtl_thresh]),"/",length(d$eQTL_spip),")\n"))
