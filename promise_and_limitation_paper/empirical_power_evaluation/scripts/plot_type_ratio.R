d = read.table("results/error_class.dat", head=T)
attach(d)
ratio = eqtl_ratio*gwas_ratio
pdf(file="results/ratio_detect_eqtl.pdf", height=4, width=7, bg="white", pointsize=18)
hist(eqtl_ratio[type==3],breaks=10, xlim=c(0,1), xlab="causal-to-lead PIP ratio", main= "eQTL (detected)")
dev.off()

pdf(file="results/ratio_detect_gwas.pdf", height=4, width=7, bg="white", pointsize=18)
hist(gwas_ratio[type==3],breaks=10, xlim=c(0,1), xlab="causal-to-lead PIP ratio", main= "GWAS (detected)")
dev.off()

pdf(file="results/ratio_c2fn_eqtl.pdf", height=4, width=7, bg="white", pointsize=18)
hist(eqtl_ratio[type==2],breaks=10, xlim=c(0,1), xlab="causal-to-lead PIP ratio", main= "eQTL (class II FN)")
dev.off()

pdf(file="results/ratio_c2fn_gwas.pdf", height=4, width=7, bg="white", pointsize=18)
hist(gwas_ratio[type==2],breaks=10, xlim=c(0,1), xlab="causal-to-lead PIP ratio", main= "GWAS (class II FN)")
dev.off()


pdf(file="results/ratio_c2fn_combined.pdf", height=4, width=7, bg="white", pointsize=18)
hist(ratio[type==2], breaks=10, xlim=c(0,1), xlab="combined PIP ratio",  main = "Class II FN")
dev.off()


pdf(file="results/ratio_detect_combined.pdf", height=4, width=7, bg="white", pointsize=18)
hist(ratio[type==3], breaks=10, xlim=c(0,1), xlab="combined PIP ratio", main="Detected")
dev.off()
