# FastENLOC Output

## Enrichment Analysis Output

### Enrichment Parameter Estimate

The enrichment parameter estimates are saved in a text file with the suffix ``enloc.enrich.out``. An example of such a file is shown below:

```
                Intercept    -7.663           -
Enrichment (no shrinkage)     4.421       0.147
Enrichment (w/ shrinkage)     4.242       0.144


## Alternative (coloc) parameterization: p1 = 4.693e-04, p2 = 7.598e-04, p12 = 2.483e-05
```
where:

+ ``Intercept``: the point estimate of $\alpha_0$
+ ``Enrichment``: the point estimate of $\alpha_1$, followed by the corresponding standard error 

Both the non-shrinkage and shrinkage versions of the $\alpha_1$ estimate are included; however, the shrinkage estimate is used for calculating colocalization probabilities.

The last line provides the equivalent transformation to the $(p_1, p_2, p_{12})$ priors used in ``coloc``.

### Multi-Imputation Details

The text tabular file with the suffix ``enloc.mi.out`` contains the parameter estimates from each multi-imputation run. 
A header is included.
The top lines from an example output are shown below:

```
  a0      a1  p_eqtl  p_gwas
 -7.646   4.416         7.696e-04       5.092e-04
 -7.648   4.347         7.919e-04       5.069e-04
 -7.648   4.332         7.813e-04       5.061e-04
 -7.646   4.354         7.691e-04       5.074e-04
 -7.644   4.296         7.813e-04       5.074e-04
 ```
 

## Colocalization Analysis Output

### SNP-level Colocalization Probability Output

The colocalization probability for each variant within a signal cluster, credible set, or locus is saved in a tab-delimited text file with the suffix ``enloc.snp.out``. A header row is included in the file. Below are the top lines from an example file:
```
Signal  SNP     PIP_qtl PIP_gwas_marginal       PIP_gwas_qtl_prior      SCP
ENSG00000006837:1(@)ENSG00000006837   ENSG00000006837_chr5_133528592_C_A_b38   4.763e-02 3.730e-05    1.491e-04      1.158e-04
ENSG00000006837:1(@)ENSG00000006837   ENSG00000006837_chr5_133528604_C_T_b38   1.391e-01 3.830e-05    3.781e-04      3.472e-04
ENSG00000006837:1(@)ENSG00000006837   ENSG00000006837_chr5_133528881_G_A_b38   1.501e-01 3.757e-05    3.974e-04      3.675e-04
ENSG00000006837:1(@)ENSG00000006837   ENSG00000006837_chr5_133529047_C_T_b38   1.501e-01 3.757e-05    3.974e-04      3.675e-04
ENSG00000006837:1(@)ENSG00000006837   ENSG00000006837_chr5_133529522_C_T_b38   8.186e-02 3.778e-05    2.340e-04      2.015e-04
```

The columns are:

+ ``Signal``: ID of the corresponding signal cluster, credible set, or locus
+ ``SNP``: variant identifier
+ ``PIP_qtl``: the posterior probability of the variant being a QTL, denoted as  $P(\gamma = 1 \mid {\rm  QTL~ data})$
+ ``PIP_gwas_marginal``: the posterior probability of the variant being a GWAS hit, denoted as $P(d = 1 \mid {\rm  GWAS~ data})$
+ ``PIP_gwas_qtl_prior``: the probability of the variant being a GWAS hit, conditional on it being a QTL, denoted as $P(d = 1 \mid {\rm GWAS~ data}, \gamma=1)$   
+ ``SCP``: SNP-level colocalization probability, denoted as $P(d=1, \gamma=1 \mid {\rm GWAS~data,~QTL~data})$

### Signal-level Colocalization Probability Output


The colocalization probability for each signal cluster, credible set, or locus is saved in a tab-delimited text file with the suffix ``enloc.sig.out``. A header row is included in the file. Below are the top lines from an example file:

```
Signal  Num_SNP CPIP_qtl        CPIP_gwas_marginal      CPIP_gwas_qtl_prior     RCP     LCP
ENSG00000006837:1(@)ENSG00000006837     18  9.894e-01 6.649e-04    3.003e-03      2.472e-03     3.053e-03
ENSG00000006837:2(@)ENSG00000006837      1  2.101e-02 1.487e-04    3.346e-04      1.978e-04     1.978e-04
ENSG00000011083:1(@)ENSG00000011083      1  9.953e-01 4.443e-05    2.881e-03      2.881e-03     2.881e-03
ENSG00000011083:2(@)ENSG00000011083     10  9.953e-01 1.865e-03    1.076e-02      9.157e-03     1.075e-02
ENSG00000011083:3(@)ENSG00000011083      4  1.531e-02 7.739e-04    9.195e-04      2.303e-04     2.383e-04
```
The columns are:

+ ``Signal``: ID of the signal cluster, constructed by concatenating the signal cluster IDs from the QTL and GWAS annotations, separated by ``(@)``
+ ``Num_SNP``: the number of the member SNPs in the signal cluster
+ ``PIP_qtl``: the cumulative posterior probability that the signal cluster contains a QTL, denoted as $\sum_i P(\gamma_i = 1 \mid {\rm  QTL~ data})$, where variant $i$ is a member SNP
+ ``CPIP_gwas_marginal``: the cumulative posterior probability that the signal cluster contains a causal GWAS hit, denoted as $\sum_i P(d_i = 1 \mid {\rm  GWAS~ data})$, where variant $i$ is a member SNP 
+ ``CPIP_gwas_qtl_prior``: the cumulative posterior probability that the signal cluster contains a causal GWAS hit, conditional on the hit also being a QTL, denoted as
$\sum_i P(d_i = 1 \mid {\rm  GWAS~ data}, \gamma_i=1)$, where variant $i$ is a member SNP 
+ ``RCP``: Signal/Regional-level colocalization probability, denoted as $\sum_i P(d_i=1, \gamma_i=1 \mid {\rm GWAS~data,~QTL~data})$, where variant $i$ is a member SNP 
+ ``LCP``: Locus-level colcalization probability, denoted as $\sum_{i,j} P(d_i=1, \gamma_j=1 \mid {\rm GWAS~data,~QTL~data})$, where variants $i$ and $j$ are member SNPs. 

Note that LCP is always no smaller than RCP.

### Gene-level Colocalization Probability Output

If the signal clusters in the QTL annotation follows the naming convention ``gene_name:signal_id``, FastENLOC will generate a gene-level colocalization probability output accumulating signal-level colocalization probabilities across the gene with the ID ``gene_name``. 
A header row is included in the file. Below are the top lines from an example file:
```
 Gene            GRCP    GLCP
ENSG00000006837         2.669e-03       3.250e-03
ENSG00000011083         1.224e-02       1.384e-02
ENSG00000013561         3.228e-04       3.228e-04
ENSG00000015479         1.092e-03       1.094e-03
ENSG00000016082         5.946e-03       7.906e-03
ENSG00000019582         1.191e-03       1.191e-03
```
The columns are:

+ ``Gene``: gene ID
+ ``GRCP``: the probability that the gene contains at least one signal cluster with RCP $ \gt 0$ 
+ ``GLCP``: the probability that the gene contains at least one signal cluster with LCP $ \gt 0$

Note that GRCP value quantifies the probability that the gene contains at least one colocalized SNP.