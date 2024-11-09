# Enrichment Analysis and Prior Setup

A key feature of FastENLOC is its ability to estimate colocalization priors directly from the provided data. Specifically, FastENLOC approaches the problem of enrichment analysis by considering the following model

$$ \log\left[ \frac{P(d = 1 \mid \gamma)}{P(d = 0 \mid \gamma)} \right] = \alpha_0 + \alpha_1 \gamma, $$

for each genetic variant, where:

+ $d$ is a binary indicator of whether the variant is a GWAS hit,
+ $\gamma$ is a binary indicator of whether the variant is a causal eQTL.
Thus, $\alpha_1$ represents the level of enrichment (log odds ratio) of causal eQTLs among GWAS hits. The input data for FastENLOC provides probabilistic estimates for $d$ and $\gamma$ across all variants, albeit imperfectly.
Treating $d$ and $\gamma$ as missing data, FastENLOC applies a multiple imputation procedure along with an EM algorithm to estimate $(\alpha_0, \alpha_1)$, yielding $(\hat \alpha_0, \hat \alpha_1)$ by pooling information across genome-wide observations. Additionally, two marginal probabilities are estimated as by-products:
	1.	Genome-wide frequency of eQTLs: $~~p_e = P(\gamma = 1)$
	2.	Genome-wide frequency of GWAS hits: $~~p_g = P(d = 1)$

## Alternative Parameterization 

The ``coloc`` method parameterizes the prior for colocalized variant by $p_{12}$. 
It follows that 

$$ p_{12} = P(d=1, \gamma=1) = \frac{\exp(\alpha_0 + \alpha_1)}{1+\exp(\alpha_0 + \alpha_1)}  p_e $$ 

Additionally, it requires

$$ p_1 = P(d=1, \gamma=0) = p_g - p_{12} $$

and 

$$ p_2 = P(d=0, \gamma=1) = p_e - p_{12} $$

Conversely,

$$ p_g = p_{12} + p_1 $$

$$ p_e = p_{12} + p_2 $$

$$ \alpha_0 = \log\left[\frac{p_1}{1-p_1-p_2-p_{12}}\right] $$

and 

$$ \alpha_1 = \log \left[ \frac{(1-p_{12}-p_1 - p_2) p_{12}}{p_1 p_2} \right],$$

indicating that the two different parameterizations are equivalent. 


## Running Enrichment Analysis

The enrichment analysis procedure is integrated into FastENLOC and runs by default. The following command-line options are key for this procedure:


+ ``-total_variants total_variants_number``: specify the total number of genetic variants measured in the GWAS study. The input data, particularly the probabilistic fine-mapping input, typically includes only a small subset of notable variants. However, the genome-wide total variant count is crucial for calibrating $p_g$ and $p_e$. If this option is not specified, FastENLOC will pause execution and prompt the user to confirm the number.
+  ``-impute imputation_runs``: define the number of multiple imputation runs, set to 25 by default, which is generally sufficient based on multiple imputation literature.

+ ``-shrinkage coef``: set the shrinkage coefficient for the $\alpha_1$ estimate. When informative colocalized variants are sparse, the $\alpha_1$ estimate can become unstable, often indicated by a large standard error. In such cases, FastENLOC shrinks $\hat \alpha_1$ toward 0 according to the specified shrinkage coefficient. This coefficient is defined as the inverse of the prior variance on $\alpha_1$ — a larger value increases the shrinkage effect. By default, the shrinkage coefficient is set to 1, which generally performs well across application scenarios.


## Estimating Enrichment Prior Only

To perform only the enrichment analysis without further calculating colocalization probabilities in FastENLOC, specify the ``--enrich_only`` option on the command line.


## Bypassing Enrichment Analysis

Users can bypass the enrichment analysis by specifying required priors directly via command-line options. While not recommended as a standard colocalization approach, this option allows sensitivity analysis of colocalization results relative to prior specification.

To bypass enrichment analysis, use one of the following options:

1.	``-a0 a0_value -a1 a1_value``: specify $\alpha_0$ and $\alpha_1$ values.
2.	``-p1 p1_value -p2 p2_value -p12 p12_value``: specify values for $p_1$, $p_2$, and $p_{12}$.

