# fastENLOC: fast enrichment estimation aided colocalization analysis

**Version 3** (planned release date: July 2023)


This repository contains the software implementation of fastENLOC, which enables integrative genetic association analysis of molecular QTL data and GWAS data. The statistical model and the key computational procedures are described in \[1\], \[2\], \[3\], and \[4\].

For questions/comments regarding the software package, please contact Xiaoquan (William) Wen (xwen at umich dot edu).


## New features in version 3

1. Unified VCF input format for both complex- and molecular-trait data
2. Improved computational methods for Bayesian fine-mapping results using signal clusters/credible sets
3. Improved multiple imputation procedure
4. Diagnosis output from multiple imputation procedure 
5. Capping option for enrichment parameter ($a_1$) 
6. Updated tutorial and documentation
7. Utility to convert fine-mapping results from other software tools (e.g., SuSiE)

## License

Software distributed under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [LICENSE](http://www.gnu.org/licenses/gpl-3.0.en.html) for more details.


## Tutorial and guideline

A detailed tutorial is provided in [``tutorial``](https://github.com/xqwen/fastenloc/tree/master/tutorial/) directory. Briefly, three main steps are required for a complete analysis

1. Summarize and prepare input from fine-mapping analysis of eQTL data 
2. Summarize and prepare input from fine-mapping analysis of GWAS data
3. Run fastENLOC

We distribute pre-computed eQTL annotations from GTEx (v8) data. In the simplest case, the required GWAS PIPs can be computed from single-SNP association summary-statistics (e.g., z-scores and p-values) using [``torus``](https://github.com/xqwen/torus/). For better accuracy and improved statistical power, both molecular and complex-trait phenotypes should be fine-mapped using individual-level data. 


## GTEx v8 multi-tissue eQTL annotations for fastENLOC

If you prefer to using newly released GTEx v8 eQTL annotation for analysis, please download the following vcf files

+  [Multi-tissue eQTL annotation with hg38 position ID](https://drive.google.com/open?id=1kfH_CffxyCtZcx3z7k63rIARNidLv1_P)
+  [Multi-tissue eQTL annotation with rs ID](https://drive.google.com/open?id=1rSaHenk8xOFtQo7VuDZevRkjUz6iwuj0)


## Citation

1. Wen, X., Pique-Regi, R., Luca, F., 2017. Integrating Molecular QTL Data into Genome-wide Genetic Association Analysis: Probabilistic Assessment of Enrichment and Colocalization. *PLOS Genetics*, 13(3): e1006646.
2. Pividori, M., et al., 2020. PhenomeXcan: Mapping the genome to the phenome through the transcriptome. *Science Advances*, 6(37), p.eaba2083.
3. Hukku, A., Pividori, M., Luca, F., Pique-Regi, R., Im, H.K. and Wen, X., 2021. Probabilistic colocalization of genetic variants from complex and molecular traits: promise and limitations. *The American Journal of Human Genetics*, 108(1), pp.25-35.
4. Hukku, A., Sampson, M.G., Luca, F., Pique-Regi, R. and Wen, X., 2022. Analyzing and Reconciling Colocalization and Transcriptome-wide Association Studies from the Perspective of Inferential Reproducibility.  *The American Journal of Human Genetics*, (in press)