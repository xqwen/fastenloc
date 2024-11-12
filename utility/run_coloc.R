library(coloc)


args = commandArgs(trailingOnly=TRUE)
enloc_summary = args[1]

d = read.table(file=enloc_summary)

run_coloc_quant_locus <-function(locus){

	ds = d[d$V1 == locus,]
	dse= list()
	dsg = list()
	
	# dse is eQTL data, dsg is GWAS data
	
	dse$snp = dsg$snp = ds$V2
	
	dse$type = dsg$type = "quant"
	dse$sdY = dsg$sdY = 1
	
	dse$beta = ds$V3
	dse$varbeta = ds$V4^2
	dsg$beta = ds$V5
	dsg$varbeta = ds$V6^2
	
	rst = coloc.abf(dataset1 = dse, dataset2=dsg)
	return(c(locus,rst$summary[6], rst$summary[6]+rst$summary[5]))

}



out=t(sapply(unique(d$V1), function(x) run_coloc_quant_locus(x)))
write(file="", t(out), ncol=3)
