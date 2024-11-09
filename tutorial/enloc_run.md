# Running Colocalization Analysis



##  Running with Probabilistic Fine-mapping Input

The command line syntax to run fastENLOC with the probabilistic fine-mapping input is

```
fastenloc -eqtl eqtl_annotation -gwas gwas_annotation -total_variants total_number_of_variants [-t tissue_name] 
```

Required options:

+ ``-eqtl``: specify the prepared eQTL annotation file
+ ``-gwas``: specify the prepared GWAS PIP file
+ ``-total_variants``: specify the number of genetic variants interrogated genome-wide
+ ``-tissue``: specify the tissue name (only required if multi-tissue eQTL annotation is used)

### Running with Legacy GWAS Annotation Files

GWAS annotation files formatted for previous versions of fastENLOC remain compatible with newer versions (3.0 and above). To use these legacy-format files, specify the option ``-go gwas_input_legacy_format`` (instead of ``-gwas``) to indicate the older GWAS annotation format.


## Running with Summary Statistics Input

For combined summary statistics files, run colocalization analysis use the following command:

```
fastenloc -sum combined_summary_input_file -total_variants
```

Required option:

+ ``-sum``: specify the combined summary statistics input file

