# fastENLOC: fast enrichment estimation aided colocalization analysis

**Version 3.1** 


This repository contains the software implementation of fastENLOC, which enables integrative genetic association analysis of molecular QTL data and GWAS data. The statistical model and the key computational procedures are described in \[1\], \[2\], \[3\], and \[4\].

For questions/comments regarding the software package, please contact Xiaoquan (William) Wen (xwen at umich dot edu).


## New features in version 3.1

- Unified VCF input format for both complex- and molecular-trait data
- Approximate computation with only summary-level data
- Improved computational methods for Bayesian fine-mapping results using signal clusters/credible sets
- Improved multiple imputation procedure with diagnostic output 
- Updated tutorial and documentation
- Utility to convert fine-mapping results from other software tools (e.g., SuSiE)

## License

Software distributed under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [LICENSE](http://www.gnu.org/licenses/gpl-3.0.en.html) for more details.


## Quick start

A detailed tutorial is provided in [``tutorial``](https://github.com/xqwen/fastenloc/tree/master/tutorial/) directory.


### Running with fine-mapping results

For best results, we recommend performing multi-SNP probabilistic fine-mapping before the fastENLOC analysis.
The fine-mapping results should be organized to two VCF files, one for each trait. (We provide utilities to convert DAP-G and SuSiE fine-mapping results to compatible fastENLOC VCF format.) 
To run fastENLOC analysis, simply issue command
```
fastENLOC -g gwas.vcf.gz -e eqtl.vcf.gz 
```
For GWAS file using the legacy fastENLOC GWAS format (version 1 & 2), use 
```
fastENLOC -go gwas.old_format.vcf.gz -e eqtl.vcf.gz
```


### Running with summary-statistic input

fastENLOC can also run with minimum summary statistics information, namely, the estimated SNP genetic effects and the corresponding standard errors from single-SNP analyses of molecular and complex traits. 
In this case, the input is organized into a single tabular file. The command to start fastENLOC analysis is 
```
fastENLOC -sum summary_stats_file
```


## GTEx v8 multi-tissue eQTL annotations for fastENLOC

We provide pre-formatted GTEx v8 eQTL annotation in fastENLOC VCF format:

+  [Multi-tissue eQTL annotation with hg38 position ID](https://drive.google.com/open?id=1kfH_CffxyCtZcx3z7k63rIARNidLv1_P)
+  [Multi-tissue eQTL annotation with rs ID](https://drive.google.com/open?id=1rSaHenk8xOFtQo7VuDZevRkjUz6iwuj0)


## Citation

1. Wen, X., Pique-Regi, R., Luca, F., 2017. Integrating Molecular QTL Data into Genome-wide Genetic Association Analysis: Probabilistic Assessment of Enrichment and Colocalization. *PLOS Genetics*, 13(3): e1006646.
2. Pividori, M., et al., 2020. PhenomeXcan: Mapping the genome to the phenome through the transcriptome. *Science Advances*, 6(37), p.eaba2083.
3. Hukku, A., Pividori, M., Luca, F., Pique-Regi, R., Im, H.K. and Wen, X., 2021. Probabilistic colocalization of genetic variants from complex and molecular traits: promise and limitations. *The American Journal of Human Genetics*, 108(1), pp.25-35.
4. Hukku, A., Sampson, M.G., Luca, F., Pique-Regi, R. and Wen, X., 2022. Analyzing and Reconciling Colocalization and Transcriptome-wide Association Studies from the Perspective of Inferential Reproducibility.  *The American Journal of Human Genetics*, (in press)