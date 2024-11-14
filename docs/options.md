# Command Line Options

Here is the list of all fastENLOC command-line options


## 1. Required

+ ``-total_variant | -tv number_of_total_variants``: required for all running options, except when $p_1$, $p_2$, and $p_{12}$ are specified. 

### 1.1 Probabilistic Fine-mapping Input
+ ``-eqtl | -e eqtl_probabilitic_annotation``
+ ``-gwas | -g gwas_probabilistic_annotation``

### 1.2 Summary Statistics Input

+ ``-summary | -sum combined_summary_input``: combined summary statistics input
+ ``-gs gwas_summary_statistics``: GWAS summary statistics that used in the hybrid option

### 1.3 Legacy GWAS probabilistic input

+ ``-go gwas_probabilistic_legacy_annotation``: for GWAS summary statistics processed by ``TORUS``

### 1.4 Tissue Specification

+ ``-tissue | -t tissue``: required only when multi-tissue QTL annotations are used

<br>
<br>


## 2. Strongly Recommended

+ ``-prefix output_prefix``: specify the output file prefix and directory. Use this option to avoid overwriting the result files with default names. 

<br>
<br>

## 3. Optional

### 3.1 Bypassing Enrichment Analysis

+ ``-p1 p1_value``: frequency of GWAS-only SNPs
+ ``-p2 p2_value``: frequency of QTL-only SNPs
+ ``-p12 p12_value``: frequency of colocalized SNPs

+ ``-coloc_default_prior``: a handy option to specify 
    - $p_1 = 1\times 10^{-4}$
    - $p_2 = 1\times 10^{-4}$
    - $p_{12} = 1\times 10^{-5}$ 

<br>

+ ``-a0 a0_value``
+ ``-a1 a1_value``


### 3.2 Bypassing Colocalization Computation (Enrichment Analysis Only)

+ ``--enrich_only | --enrich``: performing enrichment analysis only without further computing colcoalization probabilities

### 3.3 Multi-thread Processing

+ ``-thread number_of_parallel_threads``: for performance improvement. But for most application cases, single-thread processing is quite efficient

### 3.4 Colocalization Algorithm Options

+ ``--approx | --legacy ``: use legacy computational algorithm implemented in the original ENLOC described in [Wen et al., 2017](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006646).
+ ``--exact``: the default computational algorithm for computing colocalization probabilities, implemented since version 2.
+ ``-shrinkage | -s shrinkage_coefficient``: shrinkage coefficient for enrichment analysis, the default shrinkage coefficient is 1, corresponding to a ${\rm N}(0,1)$ prior on $\alpha_1$. Larger coefficient value results in stronger shrinkage. Set ``-s 0`` for no shrinkage on $\alpha_1$.

+ ``-impute | -imp multiple_imputation_rounds``: specify the number imputations for estimating the enrichment parameters. The default value is 25, which is generally sufficient. 

### 3.5 Output Options

+ ``--output_all | --all``: output colocalization probabilities for all SNPs. By default, the SNP output contains only the subset with SCP (i.e., SNP-level colocalization probability) $ \ge 10^{-4}$.

### 3.6 Miscellaneous Options

+ ``-cap_a1 a1_upper_limit``: overriding the estimated $\hat \alpha_1$ value with a pre-defined upper limit. It is mostly used for sensitivity analysis to examine how $\hat \alpha_1$ affects colocalization computation.
+ ``-seed random_seed``: user defined random seed for multi-imputation in enrichment analysis.
+ ``-sdy_e``: specify a single grid for BF computation in eQTL data
+ ``-sdy_g``: specify a single grid for BF computation in GWAS data
+ ``-conv``: convert summary statistic input to probabilistic fine-mapping format without performing any computational tasks (enrichment analysis or colocalization computation).





