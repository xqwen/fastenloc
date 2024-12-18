# Input Format


FastENLOC analysis requires summary results from genetic association analyses of a molecular trait and a complex trait to perform colocalization analysis. We refer to these summary results as eQTL and GWAS information for the molecular trait and complex trait analyses, respectively.

The current version of FastENLOC supports two types of summary information:

1. **Single-SNP Summary Statistics Input**: summary statistics from single-SNP association analyses, specifically the estimated effect size, beta, and its standard error, se(beta), for each SNP in each trait.

2. **Probabilistic Fine-mapping Input**: probabilistic association information from multi-SNP fine-mapping analyses of both traits (supported fine-mapping software includes SuSiE and DAP-G).

The probabilistic fine-mapping input is preferred, as it generally provides more accurate results for colocalization analysis, even though the computational cost for obtaining the single-SNP summary statistics input is lower.

All input files can be either compressed (by gzip) or uncompressed, FastENLOC can detect the format automatically.



## 1. Single-SNP Summary Statistics Input

The single-SNP summary statistics input should be provided in a tabular file with the following format:
```
Locus_id   SNP_id   beta_eqtl   se_eqtl    beta_gwas   se_gwas
```
The ``Locus_id`` can be used to represent a gene in an eQTL study. While the signs of ``beta_eqtl`` and ``beta_gwas`` are not critical for the colocalization analysis, it is strongly recommended that both are defined using the same reference allele.

An example input file can be downloaded from [here](https://github.com/xqwen/fastenloc/tree/master/sample_data/coloc_test_data.sum).


## 2. Probabilistic Fine-mapping Input

The input for probabilistic fine-mapping follows the standard VCF format, where association information for colocalization analysis is recorded in the INFO field. Separate files are required for molecular and complex traits. Most importantly, FastENLOC expects molecular and complex trait tiles to have matching variant IDs.
An example input line is shown below:
```
chr1	115746	chr1_115746_C_T_b38	C	T	ENSG00000269981:1@Spleen=2.00812e-01[9.997e-01:4]
```

The first five columns represent the chromosome, position, SNP ID, reference allele, and alternative allele, consistent with a standard VCF file. Note that FastENLOC uses only the SNP ID information and does not verify or utilize the position or allele information.

The INFO field has the following format:
```
locus_id@tissue=pip[cpip:number_of_snps]
```
+ ``Locus_id``: For molecular QTL mapping results, this often represents each (of possibly many) independent association signal  for a gene.
+ ``tissue`` (optional): multi-tissue/cell type QTL annotations can be saved in a single VCF file. (In GWAS file, this entry can be used to record trait name.) For a single tissue QTL study, leave out the tissue entry, i.e., use ``locus_id@=pip[cpip:number_of_snps]``
+ ``pip``: posterior inclusion probability for the variant
+ ``cpip``: cumulative posterior inclusion probability for the variant's signal cluster/credible set
+ ``number_of_snps``: number of total variants representing the corresponding signal cluster/credible set

A single variant may be associated with multiple genes across various tissues. In such cases, concatenate multiple annotations using the delimiter ``|``. For example, 
```
chr1	633264	chr1_633264_T_C_b38	T	C	ENSG00000225972:2@Muscle_Skeletal=5.77565e-01[1.000e+00:7]|ENSG00000229344:1@Liver=1.69146e-01[8.012e-01:9]
```
There is no restriction on the number of annotation entries that can be concatenated.


### 2.1 Use pre-computed GTEx multi-tissue eQTL annotation

Simply download appropriate vcf files below. Need to specify ``-tissue`` command line option in FastENLOC analysis. 

+  [Multi-tissue eQTL annotation with hg38 position ID](https://drive.google.com/open?id=1kfH_CffxyCtZcx3z7k63rIARNidLv1_P)
+  [Multi-tissue eQTL annotation with rs ID](https://drive.google.com/open?id=1rSaHenk8xOFtQo7VuDZevRkjUz6iwuj0)

When making complex trait file, make sure variant IDs match the corresponding QTL file.



### 2.2 Construct fine-mapping input from DAP-G and SuSiE results

The probabilistic fine-mapping input files can be constructed using the utility scripts provided in the [utility](../utility/) directory. 
Currently, ``FastENLOC`` supports ``DAP`` and ``SuSiE`` as both provide signal cluster/credible set information required by the FastENLOC colocalization analysis.  
Both ``DAP`` and ``SuSiE`` analyze a pre-defined genomic region, i.e., a locus, at a time. 

Upon completing fine-mapping analysis, following the procedures described below to construct the probabilistic fine-mapping input for FastENLOC,

1. Name each fine-mapping output file by ``locus_id.postfix``, where ``locus_id`` is the required entry in the FastENLOC input file.  The postfix can be arbitrary but it is required. The utility takes the leading string before the delimiter ``.`` as the locus id.
2. Organize the fine-mapping output files of all loci into a single directory ``fm_rst_dir``. The utility tool assume all files within the directories are fine-mapping output files. Thus avoid place unrelated files into ``fm_rst_dir``.
3. Provide a VCF file to all annotate all variant's position and allele information. No INFO field is expected or used by the utility tool. 



#### 2.2.1 Converting DAP output

The DAP output can be directly converted by the utility script ``dap2enloc`` by running 
```
dap2enloc -dir dap_rst_dir -vcf snp_vcf_file [-tissue tissue_name] | gzip - > fastenloc.dap.annotation.vcf.gz
```

Note that you can (but don't have to) specify the eQTL dataset name through ``-tissue`` option. This option may be useful if you consider merge eQTL annotations from multiple data sets.

See an example with data [here](dap_processing.md).

#### 2.2.2 Converting SuSiE output

For each fine-mapped locus, the SuSiE results should be outputted into a tabular file with the following format
```
CS_ID   SNP_ID   PIP
```
where
+ ``CS_ID``: signal cluster ID
+ ``SNP_ID``: snp ID
+ ``PIP``: SNP pip

The data matrix can contain more column entries, but only the first three are used by the utility script ``susie2enloc``.

The command to run the script is  
```
susie2enloc -dir susie_rst_dir -vcf snp_vcf_file [-tissue tissue_name] | gzip - > fastenloc.susie.annotation.vcf.gz
```
See an example with data [here](susie_processing.md).

## 3. Hybrid Input

In many applications, users may want to use our GTEx probabilistic annotation for eQTLs along with summary statistics for GWAS data. In earlier versions of FastENLOC, this required processing GWAS summary statistics with TORUS. However, starting from version 3.1, this step is no longer necessary, as we have integrated the primary functionality of TORUS directly into FastENLOC.

To use this functionality, users can provide GWAS summary statistics in the following tabular format:

```
Locus_id   SNP_id   beta_gwas   se_gwas
```

This format closely resembles both the combined summary statistics input and the input format for TORUS.