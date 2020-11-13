fastenloc:
	if [ ! -d "fastenloc_out" ]; then mkdir fastenloc_out; fi
	fastenloc -g sim_data/gwas.pip.gz -e sim_data/eqtl.vcf.gz -prefix fastenloc_out/sim -tv 30000000 --all
	fastenloc -g sim_data/gwas.pip.gz -e sim_data/eqtl.vcf.gz -prefix fastenloc_out/sim.true_prior -tv 30000000 --all -a0 -7.315 -a1 4.001 

