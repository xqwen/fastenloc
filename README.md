# fastENLOC: fast enrichment estimation aided colocalization analysis


This repository contains the software implementation of fastENLOC, which enables integrative genetic association analysis of molecular QTL data and GWAS data. The statistical model and the key computational procedures are described in \[[1](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006646)\] and \[[2](https://www.biorxiv.org/content/10.1101/833210v1)\]. Compared to the previous implementation of [ENLOC](https://github.com/xqwen/integrative), the new implementation is a standalone C++ program and runs magnitude faster.    

For questions/comments regarding to the software package, please contact Xiaoquan (William) Wen (xwen at umich dot edu).

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

## Building fastenloc
You will need a working C++ compiler and `make` already installed and configured. Installing additional library 
dependencies is described below.
### Debian family (Ubuntu)
```shell
sudo apt install libgsl-dev libboost-iostreams-dev zlib1g-dev
cd src
make
```

### RHEL family (CentOS, Rocky)
```shell
sudo yum install boost-devel gsl-devel zlib-devel
cd src
make
```

### Other Linux, or building libraries from source:
* [GNU Scientific Library](https://www.gnu.org/software/gsl/)
* [ZLib](https://zlib.net/)
* [Boost iostreams](https://www.boost.org/doc/libs/)

### MacOS
You will want [Homebrew](https://brew.sh/) installed to easily install libraries. It's possible to build all of the necessary
library dependencies from scratch if you prefer.

```shell
brew install gsl boost libomp
cd src
make
```


## Citation

1. Wen, X., Pique-Regi, R., Luca, F. Integrating Molecular QTL Data into Genome-wide Genetic Association Analysis: Probabilistic Assessment of Enrichment and Colocalization. *PLOS Genetics*. 2017 Mar 13(3): e1006646.
2. Pividori and Rajagopal et al. PhenomeXcan: Mapping the genome to the phenome through the transcriptome. *BioRxiv* 2019: doi: 10.1101/833210
