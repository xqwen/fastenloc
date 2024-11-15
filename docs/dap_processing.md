# Processing DAP Fine-mapping Output 

Here, we demonstrate the procedure to convert fine-mapping output from the software package DAP using an example dataset.

1. Download the sample data set from [here](https://github.com/xqwen/fastenloc/tree/master/sample_data/dap_example.tgz)

2. Unpack the package in a terminal, using the following command:

```
tar zxf dap_example.tgz
```
You should see a gzipped vcf file ``snp.chr5.vcf.gz`` and a directory ``eqtl_dap_out`` with 1,198 DAP fine-mapping output files. Note that each of the fine-mapping file is named as ``gene_name.dap.out``. The processing script will recognize and use the string ``gene_name``. The files are direct output from ``DAP``, no additional edit/processing is needed. 

3. Run the utility script [``dap2enloc``](https://github.com/xqwen/fastenloc/tree/master/utility/dap2enloc) to obtain the FastENLOC eQTL input:

```
dap2enloc -v snp.chr5.vcf.gz -d eqtl_dap_out/ > eqtl.fastenloc.vcf
```
or obtain a gzipped version:

```
dap2enloc -v snp.chr5.vcf.gz -d eqtl_dap_out/ | gzip - >  eqtl.fastenloc.vcf.gz
```


