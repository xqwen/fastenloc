# LD mismatching and its impact on colocalization analysis


## Simulation summary

In this set of simulations, we consider M = 4888263 variants across 6977 genes. The genotype data are obtained from the GEUVADIS/1000 Genome project.  

The eQTL data are simulated using the genotype data from FIN. Five GWAS data sets are simulated with the same colocalization configuration but different genotype data from 5 distinct population groups (FIN, GBR, CEU, TSI, and YRI). 

The causal GWAS variants, causal eQTL variants, and the colocalized variants are pre-determined by a sampling process. Particularly, the causal GWAS variants (and the corresponding colocalization sites) remain invariant across all simulated GWAS data sets. In total, there are 2107 genes that truly harbor a colocalized.


## Running analysis

Run ```make``` to start the analysis process. Because of the stochastic nature of the enrichment estimation (by multiple imputation), there may be slight numerical difference from the reported values in the paper.

