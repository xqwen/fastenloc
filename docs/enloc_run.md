# Running Colocalization Analysis



## 1. Running with Probabilistic Fine-mapping Input

The command line syntax to run FastENLOC with the probabilistic fine-mapping input is

```
fastenloc -eqtl eqtl_annotation -gwas gwas_annotation -total_variants total_number_of_variants [-t tissue_name] 
```

Required options:

+ ``-eqtl``: specify the prepared eQTL annotation file
+ ``-gwas``: specify the prepared GWAS PIP file
+ ``-total_variants``: specify the number of genetic variants interrogated genome-wide
+ ``-tissue``: specify the tissue name (only required if multi-tissue eQTL annotation is used)

### 1.1 Running with Legacy GWAS Annotation Files

GWAS annotation files formatted for previous versions of FastENLOC remain compatible with newer versions (3.0 and above). To use these legacy-format files, specify the option ``-go gwas_input_legacy_format`` (instead of ``-gwas``) to indicate the older GWAS annotation format.


## 2. Running with Summary Statistics Input

For combined summary statistics files, run colocalization analysis use the following command:

```
fastenloc -sum combined_summary_input_file -total_variants total_number_of_variants
```

where
+ ``-sum``: specify the combined summary statistics input file


## 3. Running with Hybrid Input

For hybrid input, where the eQTL file is in probabilistic fine-mapping format and the GWAS file is in summary statistics format, use the following command to run FastENLOC:

```
fastenloc -eqtl eqtl_annotation -gs gwas_summary_stats_file -total_variants total_number_variants [-t tissue_name]
```

Note that 

+ ``-gs``: specify the summary statistics file containing only GWAS/complex trait information. 