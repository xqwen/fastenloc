## Quick start

### Compile the source code

Run ``make`` to compile the source code in this directory. The resulting executable is named ``fastenloc``. Select the appropriate Makefile for your operating system. 

### Sample data download

+ GTEx v8 eQTL annotation: [\[with hg38 position ID\]](https://drive.google.com/open?id=1kfH_CffxyCtZcx3z7k63rIARNidLv1_P)    [\[with rs ID\]](https://drive.google.com/open?id=1rSaHenk8xOFtQo7VuDZevRkjUz6iwuj0)
+ Height GWAS z-scores in torus format: [download](https://drive.google.com/open?id=1kxZge6NQ8_8oJjVhkO4lKdmZiG2jbu1m)
+ Processing GWAS data with ``torus``

```
torus -d Height.torus.zval.gz --load_zval -dump_pip Height.gwas.pip
gzip Height.gwas.pip
```

### GWAS data input

The required GWAS input format is

```
SNP_ID  LOCUS_ID  Posterior_Inclusion_Prob
```

The PIP can be generated from any Bayesian fine-mapping algorithms (DAP-G, CAVIAR, etc.). The sample GWAS data is generated from ```torus``` based on the z-scores and the locus partition algorithm by Berisa and Pickrell.   



### Running colocalization

Command line syntax

```
fastenloc -e eqtl_annotation_gzipped -g gwas_data_gzipped -t tissue_name [-thread n] [-prefix prefix_name] [-total_variant total_snp] [-s shrinkage] [--enrich] 
```
For example to run colocalization analysis of Height GWAS and GTEx (v8) whole blood eQTLs, run
```
fastenloc -e gtex_v8.eqtl_annot.vcf.gz -g Height.gwas.pip -t Whole_Blood -prefix Height_Blood
```

#### Command line arguments

+ ```-e eqtl_input``` (required): specify the eQTL annotation input file.
                                  
+ ```-g gwas_input``` (required): specify the GWAS PIP input file.

+ ```-t tissue_name``` (required for multi-tissue eQTL input): specify the tissue name for eQTLs.

+ ```-thread n``` (optional): specify the number of simultaneous threads for multiple imputation in enrichment analysis. The default is to use a single thread.

+ ```-prefix prefix_name``` (optional): specify the alternative prefix for the output files. By default, all outputs have prefix "enloc".

+ ```-total_variant total_gwas_snps``` (optional): specify the number of total GWAS SNPs for analysis. It is important to specify this number to obtain calibrated enrichment estimates.

+ ```-s shrinkage``` (optional): the shrinkage prior for the enrichment parameter. The prior is assumed be N(0, 1/shrinkage). By default, we set shrinkage = 1. if shrinkage is set to 0 , the prior becomes flat, and there is no shrinkage effect.

+ ```--enrich``` (optional): run enrichment analysis only.


                                                               


### Output

1. Enrichment analysis result ``prefix.enloc.enrich.rst``: estimated enrichment parameters and standard errors. If shrinkage parameter is non-zero, both the unshrunk and the shrunk enrichment estimates will be output. Note that, only the shrunk estimate will be used in the colocalization computation.
2. Signal-level colocalization result ``prefix.enloc.sig.out``:  the main output from the colocalization analysis with the following format
    + column 1: signal cluster name (from eQTL analysis)
    + column 2: number of member SNPs
    + column 3: cluster PIP of eQTLs
    + column 4: cluster PIP of GWAS hits (***without*** eQTL prior)
    + column 5: cluster PIP of GWAS hits (***with*** eQTL prior)
    + column 6: regional colocalization probability (RCP)
3. SNP-level colocalization result ``prefix.enloc.snp.out``: SNP-level colocalization output with the following format
    + column 1: signal cluster name
    + column 2: SNP name
    + column 3: SNP-level PIP of eQTLs
    + column 4: SNP-level PIP of GWAS (without eQTL prior)
    + column 5: SNP-level PIP of GWAS (with eQTL prior)
    + column 6: SNP-level colocalization probability


For more detailed instructions, please refer to the [tutorial](../tutorial/).

