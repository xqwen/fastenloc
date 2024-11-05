#  Running fastENLOC with probabilistic fine-mapping input

The command line syntax to run fastENLOC with the probabilistic fine-mapping input is

```
fastenloc -eqtl eqtl_fm_annotation -gwas gwas_fm_annotation [-total_variants total_snp] [-t tissue_name] 
```

Required:

+ ``-eqtl``: specify the prepared eQTL annotation file
+ ``-gwas``: specify the prepared GWAS PIP file
+ ``-tissue``: specify the tissue name (only required if multi-tissue eQTL annotation is used)

Recommended:

+ ``-total_variants``: specify the number of total SNPs interrogated in GWAS data. This option is highly important if GWAS input does not contain all SNPs interrogated (e.g., in some cases, only fine-mapped gnomic regions are included). fastENLOC also assumes that all annotated eQTL SNPs from the eQTL annotation file are part of the GWAS analysis.   

