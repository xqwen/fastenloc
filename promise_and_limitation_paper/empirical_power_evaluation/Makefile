fastenloc:
	if [ ! -d "fastenloc_out" ]; then mkdir fastenloc_out; fi
	fastenloc -g sim_data/gwas.pip.gz -e sim_data/eqtl.vcf.gz -prefix fastenloc_out/sim -tv 30000000 --all
	fastenloc -g sim_data/gwas.pip.gz -e sim_data/eqtl.vcf.gz -prefix fastenloc_out/sim.true_prior -tv 30000000 --all -a0 -7.315 -a1 4.001 
analysis:
	if [ ! -d "results" ]; then mkdir results; fi
	perl scripts/get_true_coloc_gene.pl > results/sim.coloc.gene.list
	perl scripts/summarize_fastenloc_estp.pl > results/fastenloc.est_prior.summary
	perl scripts/summarize_fastenloc_truep.pl > results/fastenloc.true_prior.summary
	echo "Power Analysis for fastENLOC"
	Rscript scripts/fastenloc_power.R
	Rscript scripts/coloc_power.R
	echo "Generate plots"
	Rscript scripts/plot_scatter.R 
	if [ ! -d "dap_eqtl" ]; then exit; fi
	perl scripts/classify_snp.pl > results/error_class.dat
	Rscript scripts/plot_type_ratio.R

