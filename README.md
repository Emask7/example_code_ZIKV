# This is an example of how to analyze bulk RNA-seq data, using raw counts from the 1543PC study of ZIKV in baboons
The analysis follows methods detailed in the DESeq2 user manual (http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html), but with a few modifications

Run analysis from the following files in this order:
1. data_import.R
      In this file, you import raw gene counts and sample information from Excel files, then run PCA to identify outliers in the data
2. 
