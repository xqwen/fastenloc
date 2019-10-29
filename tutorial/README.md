# fastENLOC User's Guide

To run fastENLOC, one needs to prepare probabilistic eQTL annotations generated from software package ``[DAP-G](https://github.com/xqwen/dap/)`` and posterior probabilities from analyzing GWAS data. This document illustrates the details on each step.


## 1. Preparing probabilistic eQTL annotation

### Use pre-computed GTEx multi-tissue eQTL annotation

Simply download appropriate vcf files below. Note that the variant IDs have to match in the GWAS and eQTL annotations. 

+  [Multi-tissue eQTL annotation with hg38 position ID](https://drive.google.com/open?id=1kfH_CffxyCtZcx3z7k63rIARNidLv1_P)
+  [Multi-tissue eQTL annotation with rs ID](https://drive.google.com/open?id=1rSaHenk8xOFtQo7VuDZevRkjUz6iwuj0)


### Derive annotations based on your own eQTL data

The required annotations can be obtained from Bayeisan multi-SNP fine-mapping analysis using software package [``DAP``](https://github.com/xqwen/dap/). This [document](https://github.com/xqwen/dap/tree/master/gtex_v8_analysis) provides a summary on the processing of GTEx v8 data. 

Once the fine-mapping by DAP is completed

1. make sure each result file from a gene is named by the corresponding gene name, e.g., ``GENE\_NAME.postfix`` (postfix can be arbitrary but it is required)
2. Put all fine-mapping results into a single empty directory (``dap_rst_dir``)
3. provide a vcf file to annotate all SNP's positions (``snp_vcf_file``).
4. Run script ``summarize_dap2enloc.pl`` to generate annotation vcf file

```
   summarize_dap2enloc.pl -dir dap_rst_dir -vcf snp_vcf_file [-tissue tissue_name] | gzip - > fastenloc.eqtl.annotation.vcf.gz
```

Note that you can (but don't have to) specify the eQTL dataset name through ``-tissue`` option. This option may be useful if you consider merge eQTL annotations from multiple data sets.


## 2. Preparing GWAS PIP input

The required GWAS PIP input file should be a gzipped text file with the following format 

+ column 1: SNP ID
+ column 2: Signal ID
+ column 3: SNP PIP

Additional columns are allowed, but the information in those columns will not be used. 





### Use summary-level GWAS statistics

The minimum requirement for GWAS data are single-SNP testing summary statistics, e.g. the z-scores for each SNPs. All the input data for a given GWAS trait should be organized into a single gzipped text file with the following simple format

```
    chr1_13550_G_A_b38  Loc1  -0.18304155192313465
    chr1_14671_G_C_b38  Loc1  0.039355713780289986
    chr1_14677_G_A_b38  Loc1  -0.33373972601437063
    chr1_14933_G_A_b38  Loc1  -0.04685911993306961
    chr1_16841_G_T_b38  Loc1  0.1852334931031181
    chr1_17005_A_G_b38  Loc1  -0.01500239015749479
    chr1_17147_G_A_b38  Loc1  0.19543816258120805
    chr1_17407_G_A_b38  Loc1  0.004952064233334381
```
The first and second columns represent the ID the LD block of the corresponding SNP, respectively. The last column i
ndicates the z-scores. The LD blocks are defined based on the results of [Berisa and Pickrell, 2015](http://bioinfor
matics.oxfordjournals.org/content/32/2/283). The segmentation of the LD blocks are population-specific, the detailed
 information can be found [here](https://bitbucket.org/nygcresearch/ldetect-data).

For reference, we provide a complete sample data from [Height GWAS](https://drive.google.com/open?id=1kxZge6NQ8_8oJjVhkO4lKdmZiG2jbu1m) The data set is orignially from UK Biobank with additional SNPs imputed to match the GTEx SNP panel.

To convert the z-scores to PIPs, run 

```
torus -d Height.torus.zval.gz --load_zval -dump_pip Height.gwas.pip
gzip Height.gwas.pip
```
The resulting ``Height.gwas.pip.gz`` is ready to be used in fastENLOC


### Use fine-mapped GWAS data

We highly recommend to perform multi-SNP fine-mapping analysis on your GWAS data whenever is possible. Importantly, SNPs should be grouped into fine-mapped loci, and corresponding locus names should be used as the signal IDs. For this reason, we recommend fine-mapping methods that can automatically construct Bayesian credible sets, e.g., [DAP](https://github.com/dap/) and [susieR](https://github.com/stephenslab/susieR). 


**Importantly, GWAS input file does not have to include all SNPs interrogated. If an eQTL SNP is missing from the GWAS PIP input, fastENLOC assumes the corresponding PIP is close to 0**. It is allowed to only include important GWAS signals in the GWAS input file.



## 3. Running fastENLOC

The command line syntax to run fastENLOC is 
```
fastenloc -eqtl eqtl_annotation_gzipped -gwas gwas_data_gzipped [-total_variants total_snp] [-t tissue_name] [-thread n] [-prefix prefix_name] [-pc pseudo_count]
```

Required:

+ ``-eqtl``: specify the prepared eQTL annotation file
+ ``-gwas``: specify the prepared GWAS PIP file
+ ``-tissue``: specify the tissue name (only required if multi-tissue eQTL annotation is used)

Recommended:

+ ``-total_varaints``: specify the number of total SNPs interrogated. This option is important if GWAS or eQTL file does not contain all SNPs interrogated. 

Optional:

+ ``-thread``: number of parallel threads for analysis. By default, a single thread is used
+ ``-pc``: specify the shrinkage prior/pseudo count in the enrichment analysis. By default, pc = 1. increase pc values if the enrichment estimates have large estimated variance. 
+ ``-prefix``: specify the prefix for the output files.




## Output from  fastENLOC

1. Enrichment analysis result ``prefix.enloc.enrich.rst``: estimated enrichment parameters and standard errors.
2. Signal-level colocalization result ``prefix.enloc.sig.out``:  the main output from the colocalization analysis wi
th the following format
    + column 1: signal cluster name (from eQTL analysis)
    + column 2: number of member SNPs
    + column 3: cluster PIP of eQTLs
    + column 4: cluster PIP of GWAS hits (***without*** eQTL prior)
    + column 5: cluster PIP of GWAS hits (***with*** eQTL prior)
    + column 6: regional colocalization probability (RCP)
3. SNP-level colocalization result ``prefix.enloc.snp.out``: SNP-level colocalization output with the following form
at
    + column 1: signal cluster name
    + column 2: SNP name
    + column 3: SNP-level PIP of eQTLs
    + column 4: SNP-level PIP of GWAS (without eQTL prior)
    + column 5: SNP-level PIP of GWAS (with eQTL prior)
    + column 6: SNP-level colocalization probability

Use 
```
 sort -grk6 prefix.enloc.sig.out 
```
to get a sorted list of colocalization signals.
