# Bayesian False Discovery Rate Control of Colocalization Probabilities


A frequently asked question in probabilistic colocalization analysis is how to control false positives using the posterior probability output. A common practice in bioinformatics is to consider observations with posterior probabilities above a predefined threshold (e.g., 0.8) as noteworthy. This approach effectively controls the local false discovery rate (fdr) at $1 - \text{threshold}$. However, controlling the local fdr is not directly comparable to the false discovery rate (FDR) control commonly used in multiple testing settings.

The process of controlling the False Discovery Rate (FDR) using posterior probabilities is both simple and straightforward to implement. Detailed technical discussions of the methodology can be found in [Efron's book](https://efron.ckirby.su.domains/other/2010LSIexcerpt.pdf), [Newton et al., 2004](https://pubmed.ncbi.nlm.nih.gov/15054023/), and [Stephens, 2017](https://academic.oup.com/biostatistics/article-abstract/18/2/275/2557030?redirectedFrom=fulltext). Here, we focus on the computational procedure using GRCP output from the FastENLOC analysis.

The sample gene-level output file ([download](https://github.com/xqwen/fastenloc/tree/master/sample_data/okamoto_sim.enloc.gene.out)) is from analyzing simulated data in [Okamoto et al. 2023](https://www.cell.com/ajhg/fulltext/S0002-9297(22)00536-5). The file contains GRCP and GLCP values (generated using probabilistic fine-mapping input) for 1,198 genes.


## Bayesian FDR Control using R 

To perform Bayesian FDR control using GRCP values, first load data into R:

```
attach(read.table("okamoto_sim.enloc.gene.out", head=T))
summary(GRCP)
```

The output shows the summary statistics of GRCP values:

```
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
0.0000776 0.0006596 0.0014740 0.0670592 0.0032180 1.0000000
```


The FDR control procedure finds the probability threshold corresponding to the pre-defined FDR control level. In this example we set the control level at 5%.

```
control_level <- 0.05
```
The following R code identifies the corresponding GRCP probability threshold:

```
lfdr <- sort(1-GRCP)
FDR <- cumsum(lfdr)/1:length(lfdr)
thresh <- 1 - lfdr[max(which(FDR<=control_level))]
```

In this example,
```
> thresh
[1] 0.8196  
```

Thus, we reject all the genes with GRCP values $ \ge 0.8196$, i.e.,

```
rej <- Gene[GRCP >= thresh]
length(rej)

> [1] 60
```

Thus, at 5% FDR level, we identify 60 genes harboring a colocalized variant. 


## Conclusions

The Bayesian FDR control procedure described above can be broadly applied to various scenarios involving posterior probabilities, without the need of converting to $p$-values. The core computation requires only four lines of R code.

For FastENLOC, user can also perform FDR control for signal-level colocalization output, such as RCPs and LCPs. 

