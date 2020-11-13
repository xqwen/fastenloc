truth  = read.table("true.coloc.gene")$V1


d = read.table("fin_vs_fin.truep.enloc.gene.out")
lfdr = sort(1-d$V2)
FDR = cumsum(lfdr)/1:length(lfdr)
thresh = 1 - lfdr[max(which(FDR<=0.05))]
rej = d$V1[d$V2>thresh]
int = intersect(rej,truth)
(length(rej)-length(int))/length(rej)
length(int)/length(truth)
length(rej)
print("====")


d = read.table("fin_vs_gbr.truep.enloc.gene.out")
lfdr = sort(1-d$V2)
FDR = cumsum(lfdr)/1:length(lfdr)
thresh = 1 - lfdr[max(which(FDR<=0.05))]
rej = d$V1[d$V2>thresh]
int = intersect(rej,truth)
(length(rej)-length(int))/length(rej)
length(int)/length(truth)




d = read.table("fin_vs_tsi.truep.enloc.gene.out")
lfdr = sort(1-d$V2)
FDR = cumsum(lfdr)/1:length(lfdr)
thresh = 1 - lfdr[max(which(FDR<=0.05))]
rej = d$V1[d$V2>thresh]
int = intersect(rej,truth)
(length(rej)-length(int))/length(rej)
length(int)/length(truth)



d = read.table("fin_vs_ceu.truep.enloc.gene.out")
lfdr = sort(1-d$V2)
FDR = cumsum(lfdr)/1:length(lfdr)
thresh = 1 - lfdr[max(which(FDR<=0.05))]
rej = d$V1[d$V2>thresh]
int = intersect(rej,truth)
(length(rej)-length(int))/length(rej)
length(int)/length(truth)



d = read.table("fin_vs_yri.truep.enloc.gene.out")
lfdr = sort(1-d$V2)
FDR = cumsum(lfdr)/1:length(lfdr)
thresh = 1 - lfdr[max(which(FDR<=0.05))]
rej = d$V1[d$V2>thresh]
int = intersect(rej,truth)
(length(rej)-length(int))/length(rej)
length(int)/length(truth)

