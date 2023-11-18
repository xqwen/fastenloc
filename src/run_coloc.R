library(coloc)
args = commandArgs(trailingOnly=T)

data = read.table(args[1])
loci = unique(data[,1])


run_coloc <-function(d) {
    snp = d[,2]
    beta1 = d[,3]
    varbeta1 = d[,4]^2

    beta2 = d[,5]
    varbeta2 = d[,6]^2

    d1 = list()
    d1$beta = beta1
    d1$varbeta = varbeta1
    d1$sdY = 1
    d1$type = "quant"
    d1$snp = snp

    d2 = list()
    d2$beta = beta2
    d2$varbeta = varbeta2
    d2$sdY = 1
    d2$type = "quant"
    d2$snp = snp

    rst = coloc.abf(dataset1=d1, dataset2=d2)

    cat(d[1,1],"\t RCP = ",rst$summary[6], "\t LCP = ", rst$summary[5]+rst$summary[6],"\n")  
    


}


sapply(loci, function(x) run_coloc(data[which(data[,1]==x),]))




