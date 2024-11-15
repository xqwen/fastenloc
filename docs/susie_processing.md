# Processing DAP Fine-mapping Output 

Here, we demonstrate the procedure to convert fine-mapping output from the software ``SuSiE`` using an example dataset.

1. Download the sample data set from [here](https://github.com/xqwen/fastenloc/tree/master/sample_data/susie_example.tgz)

2. Unpack the package in a terminal, using the following command:

```
tar zxf susie_example.tgz
```
You should see a gzipped vcf file ``snp.chr5.vcf.gz`` and a directory ``eqtl_susie_out`` with 1,198 SuSiE fine-mapping output files. Note that each of the fine-mapping file is named as ``gene_name.dap.out``. The processing script will recognize and use the string ``gene_name``. 

Each of the SuSiE fine-mapping file is extracted from the SuSiE R object with the following tabular format with no header needed.

```
Credible_Set SNP PIP
```


3. Run the utility script [``susie2enloc``](https://github.com/xqwen/fastenloc/tree/master/utility/susie2enloc) to obtain the FastENLOC eQTL input:

```
susie2enloc -v snp.chr5.vcf.gz -d eqtl_susie_out/ > eqtl.fastenloc.vcf
```
or obtain a gzipped version:

```
susie2enloc -v snp.chr5.vcf.gz -d eqtl_susie_out/ | gzip - >  eqtl.fastenloc.vcf.gz
```
