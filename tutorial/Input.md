

# Input Format


FastENLOC analysis requires summary results from genetic association analyses of a molecular trait and a complex trait to perform colocalization analysis. We refer to these summary results as eQTL and GWAS information for the molecular trait and complex trait analyses, respectively.

The current version of fastENLOC supports two types of summary information:

1. **Probabilistic Fine-mapping Input**: probabilistic association information from multi-SNP fine-mapping analyses of both traits (supported fine-mapping software includes SuSiE and DAP-G).

2. **Single-SNP Summary Statistics Input**: summary statistics from single-SNP association analyses, specifically the estimated effect size, beta, and its standard error, se(beta), for each SNP in each trait.

The first type of input is generally preferred, as it tends to yield more accurate colocalization analysis results. However, the computational cost to obtain the second type of input is considerably lower.



## Single-SNP summary statistics input

The single-SNP summary statistics input requires a text file, which can be either compressed or uncompressed. The file should follow the format of a data matrix:
```
Locus_id   SNP_id   beta_eqtl   se_eqtl    beta_gwas   se_gwas
```
The ``Locus_id`` can be used to represent a gene in an eQTL study. While the signs of ``beta_eqtl`` and ``beta_gwas`` are not critical for the colocalization analysis, it is strongly recommended that both are defined using the same reference allele.

An example input file can be downloaded from [here]().


## Probabilistic fine-mapping input

To run fastENLOC, one needs to prepare probabilistic eQTL annotations generated from software package [``DAP-G``](https://github.com/xqwen/dap/) and posterior probabilities from analyzing GWAS data. This document illustrates the details on each step.


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

