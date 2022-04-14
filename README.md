# fastENLOC: fast enrichment estimation aided colocalization analysis

**Current release: version 2** (April, 2022)


This repository contains the software implementation of fastENLOC, which enables integrative genetic association analysis of molecular QTL data and GWAS data. The statistical model and the key computational procedures are described in \[1\], \[2\], \[3\], and \[4\].

For questions/comments regarding to the software package, please contact Xiaoquan (William) Wen (xwen at umich dot edu).

## New features in version 2

+ Locus-level colcoalization analysis
+ Auto diagnosis of input files
+ Utility scripts for generating input files and computing gene-level colocalization probabilities (GLCP and GRCP). 
+ Code optimization
+ Bug fix

## License

Software distributed under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [LICENSE](http://www.gnu.org/licenses/gpl-3.0.en.html) for more details.


## Tutorial and guideline

A detailed tutorial is provided in [``tutorial``](https://github.com/xqwen/fastenloc/tree/master/tutorial/) directory. Briefly, three main steps are required for a complete analysis

1. Prepare eQTL annotation
2. Prepare GWAS sumary (in term of posterior inclusion probabilities, or PIPs)
3. Run fastenloc

We distribute pre-computed eQTL annotations from GTEx (v8) data. In the simplest case, the required GWAS PIPs can be computed from single-SNP association summary-statistics (e.g., z-scores and p-values) using [``torus``](https://github.com/xqwen/torus/)


## GTEx v8 multi-tissue eQTL annotations for fastENLOC

If you prefer to using newly released GTEx v8 eQTL annotation for analysis, please download the following vcf files

+  [Multi-tissue eQTL annotation with hg38 position ID](https://drive.google.com/open?id=1kfH_CffxyCtZcx3z7k63rIARNidLv1_P)
+  [Multi-tissue eQTL annotation with rs ID](https://drive.google.com/open?id=1rSaHenk8xOFtQo7VuDZevRkjUz6iwuj0)


## Citation

1. Wen, X., Pique-Regi, R., Luca, F., 2017. Integrating Molecular QTL Data into Genome-wide Genetic Association Analysis: Probabilistic Assessment of Enrichment and Colocalization. *PLOS Genetics*, 13(3): e1006646.
2. Pividori, M., et al., 2020. PhenomeXcan: Mapping the genome to the phenome through the transcriptome. *Science Advances*, 6(37), p.eaba2083.
3. Hukku, A., Pividori, M., Luca, F., Pique-Regi, R., Im, H.K. and Wen, X., 2021. Probabilistic colocalization of genetic variants from complex and molecular traits: promise and limitations. *The American Journal of Human Genetics*, 108(1), pp.25-35.
4. Hukku, A., Sampson, M.G., Luca, F., Pique-Regi, R. and Wen, X., 2022. Analyzing and Reconciling Colocalization and Transcriptome-wide Association Studies from the Perspective of Inferential Reproducibility.  *The American Journal of Human Genetics*, (in press)