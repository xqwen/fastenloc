# FastENLOC: Fast Enrichment Estimation Aided Colocalization Analysis

**Version 3.1** 


This repository contains the software implementation of FastENLOC, which enables integrative genetic association analysis of molecular QTL data and GWAS data. The statistical model and the key computational procedures are described in \[1\], \[2\], \[3\], and \[4\].

For questions/comments regarding the software package, please contact Xiaoquan (William) Wen (xwen at umich dot edu).


## New Features in Version 3.1

- Unified VCF input format for both complex- and molecular-trait data
- Colocalization computation with only summary-level data (complete compatibility with ``coloc``)
- Improved computational algorithms for Bayesian fine-mapping results using signal clusters/credible sets
- Improved multiple imputation procedures with diagnostic output 
- Utility to convert fine-mapping results from SuSiE and DAP-G
- Complete tutorial and documentation



## License


Software distributed under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [LICENSE](http://www.gnu.org/licenses/gpl-3.0.en.html) for more details.


## Quick Start

Visit the [``tutorial``](https://xqwen.github.io/fastenloc/) for a detailed guide and use case examples.


### Running with Fine-mapping Results

For optimal results, we recommend performing multi-SNP probabilistic fine-mapping prior to running the FastENLOC analysis. The fine-mapping outputs should be organized into two VCF files, one for each trait. (Utilities are available to convert DAP-G and SuSiE fine-mapping results to the FastENLOC-compatible VCF format.)

To run the FastENLOC analysis, use the following command:
```
fastenloc -g gwas.vcf.gz -e eqtl.vcf.gz -tv total_variants
```
If you are using the legacy FastENLOC GWAS format (version 1 or 2) for the GWAS file, use this command instead:
```
fastenloc -go gwas.old_format.vcf.gz -e eqtl.vcf.gz -tv total_variants
```


### Running with Summary Statistic Input

FastENLOC can also operate with minimal summary statistics, specifically the estimated SNP genetic effects and corresponding standard errors from single-SNP analyses of molecular and complex traits.

For this approach, the input should be organized into a single tabular file. To initiate the FastENLOC analysis, use the following command:
```
fastenloc -sum summary_stats_file -tv total_variants
```


## GTEx v8 Multi-tissue eQTL Annotations for FastENLOC
We provide pre-formatted GTEx v8 eQTL annotations (fine-mapped using DAP-G) in the FastENLOC-compatible VCF format.

+  [Multi-tissue eQTL annotation with hg38 position ID](https://drive.google.com/open?id=1kfH_CffxyCtZcx3z7k63rIARNidLv1_P)
+  [Multi-tissue eQTL annotation with rs ID](https://drive.google.com/open?id=1rSaHenk8xOFtQo7VuDZevRkjUz6iwuj0)


## Citation

1. Wen, X., Pique-Regi, R., Luca, F., 2017. Integrating Molecular QTL Data into Genome-wide Genetic Association Analysis: Probabilistic Assessment of Enrichment and Colocalization. *PLOS Genetics*, 13(3): e1006646.
2. Pividori, M., et al., 2020. PhenomeXcan: Mapping the genome to the phenome through the transcriptome. *Science Advances*, 6(37), p.eaba2083.
3. Hukku, A., Pividori, M., Luca, F., Pique-Regi, R., Im, H.K. and Wen, X., 2021. Probabilistic colocalization of genetic variants from complex and molecular traits: promise and limitations. *The American Journal of Human Genetics*, 108(1), pp.25-35.
4. Hukku, A., Sampson, M.G., Luca, F., Pique-Regi, R. and Wen, X., 2022. Analyzing and Reconciling Colocalization and Transcriptome-wide Association Studies from the Perspective of Inferential Reproducibility.  *The American Journal of Human Genetics*, (in press)
