fastenloc:
	if [ ! -d "fastenloc_out" ]; then mkdir fastenloc_out; fi
	fastenloc -g sim_data/gwas_fin.ld_sim.pip.gz -e sim_data/eqtl.fin.vcf.gz -prefix fastenloc_out/fin_vs_fin -tv 4888263 --all
	fastenloc -g sim_data/gwas_gbr.ld_sim.pip.gz -e sim_data/eqtl.fin.vcf.gz -prefix fastenloc_out/fin_vs_gbr -tv 4888263 --all
	fastenloc -g sim_data/gwas_ceu.ld_sim.pip.gz -e sim_data/eqtl.fin.vcf.gz -prefix fastenloc_out/fin_vs_ceu -tv 4888263 --all
	fastenloc -g sim_data/gwas_tsi.ld_sim.pip.gz -e sim_data/eqtl.fin.vcf.gz -prefix fastenloc_out/fin_vs_tsi -tv 4888263 --all
	fastenloc -g sim_data/gwas_yri.ld_sim.pip.gz -e sim_data/eqtl.fin.vcf.gz -prefix fastenloc_out/fin_vs_yri -tv 4888263 --all
analysis:
	perl scripts/parse_enloc_gene.pl fastenloc_out/fin_vs_fin.enloc.sig.out > fastenloc_out/fin_vs_fin.enloc.gene.out
	perl scripts/parse_enloc_gene.pl fastenloc_out/fin_vs_gbr.enloc.sig.out > fastenloc_out/fin_vs_gbr.enloc.gene.out
	perl scripts/parse_enloc_gene.pl fastenloc_out/fin_vs_ceu.enloc.sig.out > fastenloc_out/fin_vs_ceu.enloc.gene.out
	perl scripts/parse_enloc_gene.pl fastenloc_out/fin_vs_tsi.enloc.sig.out > fastenloc_out/fin_vs_tsi.enloc.gene.out
	perl scripts/parse_enloc_gene.pl fastenloc_out/fin_vs_yri.enloc.sig.out > fastenloc_out/fin_vs_yri.enloc.gene.out
