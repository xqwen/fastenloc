truth  = read.table("sim_data/true.coloc.gene")$V1


d = read.table("fastenloc_out/fin_vs_fin.enloc.gene.out")
lfdr = sort(1-d$V2)
FDR = cumsum(lfdr)/1:length(lfdr)
thresh = 1 - lfdr[max(which(FDR<=0.05))]
rej = d$V1[d$V2>thresh]
int = intersect(rej,truth)
(length(rej)-length(int))/length(rej)
length(int)/length(truth)
length(rej)
print("====")


d = read.table("fastenloc_out/fin_vs_gbr.enloc.gene.out")
lfdr = sort(1-d$V2)
FDR = cumsum(lfdr)/1:length(lfdr)
thresh = 1 - lfdr[max(which(FDR<=0.05))]
rej = d$V1[d$V2>thresh]
int = intersect(rej,truth)
(length(rej)-length(int))/length(rej)
length(int)/length(truth)
length(rej)
print("====")



d = read.table("fastenloc_out/fin_vs_tsi.enloc.gene.out")
lfdr = sort(1-d$V2)
FDR = cumsum(lfdr)/1:length(lfdr)
thresh = 1 - lfdr[max(which(FDR<=0.05))]
rej = d$V1[d$V2>thresh]
int = intersect(rej,truth)
(length(rej)-length(int))/length(rej)
length(int)/length(truth)
length(rej)
print("====")

d = read.table("fastenloc_out/fin_vs_ceu.enloc.gene.out")
lfdr = sort(1-d$V2)
FDR = cumsum(lfdr)/1:length(lfdr)
thresh = 1 - lfdr[max(which(FDR<=0.05))]
rej = d$V1[d$V2>thresh]
int = intersect(rej,truth)
(length(rej)-length(int))/length(rej)
length(int)/length(truth)
length(rej)
print("====")



d = read.table("fastenloc_out/fin_vs_yri.enloc.gene.out")
lfdr = sort(1-d$V2)
FDR = cumsum(lfdr)/1:length(lfdr)
thresh = 1 - lfdr[max(which(FDR<=0.05))]
rej = d$V1[d$V2>thresh]
int = intersect(rej,truth)
(length(rej)-length(int))/length(rej)
length(int)/length(truth)
length(rej)
print("====")
