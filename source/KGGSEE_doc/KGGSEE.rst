.. raw:: html

    <p style="font-size: 28px; line-height: 35px; font-weight: bold" align=center>
        KGGSEE: A biological Knowledge-based mining platform for Genomic and Genetic association Summary statistics using gEne Expression
    </p>

.. centered:: User manual 1.0
.. centered:: Miaoxin Li, Lin Jiang

Introduction
=============

KGGSEE is a standalone Java tool for knowledge-based analyses of genomic and genetic association summary statistics of complex phenotypes by integrating gene expression and related data. It has four major integrative analyses, 1) a multi-strategy conditional gene-based association framework, 2) conditional gene-based association guided by selective expression in tissues or cell types, 3) estimation of phenotype-associated tissues or cell-type based on gene expression in single-cell or bulk cells of different tissues, and 4) causal gene inference for complex diseases and/or traits based-on multiple eQTL. More integrative analysis functions will be added into this analysis platform in the future.

.. image:: ./media/kggsee_pipeline.png
    :align: center

Installation
==============

kggsee.jar
~~~~~~~~~~~~~~~

The main library of KGGSEE, kggsee.jar, does not need an installation procedure as long as its `Java Runtime Environment(JRE) v1.8 <https://www.oracle.com/java/technologies/javase-jre8-downloads.html>`_ (or up) is pre-installed in your machine. The kggsee.jar can be directly executed given a file path

R packages
~~~~~~~~~~~~~

.. note::
    However, you many need to install several R packages which will be called by kggsee.jar. The following are instructions for R package installation under the interactive working interface of R.

- Rserve

    ``install.packages(“Rserve”, dep=TRUE)``
    More installation instruction can be found at https://www.dundas.com/support/support-center/support-articles/installation/install-and-configure-r.

- NNLM

    ``install.packages(“NNLM”, dep=TRUE)``
    More installation instruction can be found at https://github.com/linxihui/NNLM.

- MendelianRandomization

    The first step is to install the PhenoScanner package(and the MendelianRandomization package if you haven't done this previously):

    .. code:: R

        install.packages("devtools")
        library(devtools)
        install_github("phenoscanner/phenoscanner")
        library(phenoscanner)

        install.packages("MendelianRandomization")
        library(MendelianRandomization)

Resource data
~~~~~~~~~~~~~~~~~~

Under the folder of kggsee, there is a folder named **resources**, which contains running resource data, e.g., gene boundary and gene expression. KGGSEE will automatically download required resource files. Users can all manually download the files and put them into the corresponding folders.

Tutorials
=========================

Gene-based association analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Purpose: Detect associated genes of a phenotype by GWAS summary statistics
 - Input data:
     
   1. GWAS summary statistics compressed in a text file(a fabled data set for education purpose): *examples/gwas.sum.stat.gz*
     
   2. Genotypes in VCF format to approximate correction between summary statistics: *resources/hg19/gty/1kg.phase3.v5.shapeit2.eur.hg19.chr1.vcf.gz*

.. code:: shell 

  java -Xmx10g -jar kggsee.jar --gene-assoc --sum-file examples/gwas.sum.stat.gz --vcf-ref resources/hg19/gty/1kg.phase3.v5.shapeit2.eur.hg19.chr1.vcf.gz --keep-ref --nt 10 --out examples/out/geneAssoc


Gene-based causality analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Purpose: Detect causal genes of a phenotype by GWAS summary statistics and eQTL
 - Input data:

   1. GWAS summary statistics compressed in a text file(a fabled data set for education purpose): *examples/gwas.sum.stat.gz*
   
   2. Genotypes in KGGSEE objects(generated last time) to approximate correction between summary statistics: *examples/out/geneAssoc*
   
   3. eQTL summary statistics compressed in a text file: *resources/hg19/eqtl/Brain-FrontalCortex_BA9_.transcript.maf05.p05.gz.eqtl.txt.gz*
     
.. code:: shell  

   java -Xmx10g -jar kggsee.jar --macg --eqtl-file resources/hg19/eqtl/Brain-FrontalCortex_BA9_.transcript.maf05.p05.gz.eqtl.txt.gz --filter-maf-le 0.05 --sum-file examples/gwas.sum.stat.gz --beta-or y --saved-ref  examples/out/geneAssoc --nt 10 --out examples/out/macg --excel


Estimate relevant cell-types of a phenotype
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
Purpose: Estimate relevant cell-types of a phenotype and finely map associated genes according to selective expression
 - Input data:
    
   1. GWAS summary statistics compressed in a text file(a fabled data set for education purpose): *examples/gwas.sum.stat.gz*
     
   2. Genotypes in KGGSEE objects(generated last time) to approximate correction between summary statistics: *examples/out/geneAssoc*
     
   3. Gene expression data compressed in a text file: *resources/hs_scRNA_cluster_mean.tsv.gz*
     
.. code:: shell

     java -Xmx10g -jar kggsee.jar --spa --expression-file resources/hs_scRNA_cluster_mean.tsv.gz --only-hgnc-gene --sum-file examples/gwas.sum.stat.gz --saved-ref  examples/out/geneAssoc --nt 10 --out examples/out/spa --excel
 

eQTL-guided gene-based association analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Purpose: Detect associated genes of a phenotype by GWAS summary statistics guided by eQTLs of a celltype
 - Input data:
     
   1. GWAS summary statistics compressed in a text file(a fabled data set for education purpose): *examples/gwas.sum.stat.gz*
     
   2. Genotypes in KGGSEE objects(generated last time) to approximate correction between summary statistics: *examples/out/geneAssoc*

   3. Gene expression data compressed in a text file: *resources/gtex.v8.transcript.mean.tsv.gz*

   4. eQTL summary statistics compressed in a text file: *resources/hg19/eqtl/Brain-FrontalCortex_BA9_.transcript.maf05.p05.gz.eqtl.txt.gz*

.. code:: shell

  java -Xmx10g -jar kggsee.jar --gene-finemapping --sum-file examples/gwas.sum.stat.gz --saved-ref  examples/out/geneAssoc --expression-file resources/gtex.v8.transcript.mean.tsv.gz --eqtl-file resources/hg19/eqtl/Brain-FrontalCortex_BA9_.transcript.maf05.p05.gz.eqtl.txt.gz --filter-eqtl-p 0.01 --nt 10 --out examples/out/geneAssoceQTL


Functions
=========================

Gene-based association analysis by an effective chi-square statistics(ECS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One can perform gene-based association analysis by an effective chi-square statistics(ECS) with GWAS *p*-values of variants. The *p*-values are converted to chi-square statistics(degree of freedom = 1). The ECS merges all chi-square statistics of a gene after correcting the redundancy of the statistics due to LD. The LD is calculated from genotypes of an ancestrally matched sample in VCF format, e.g. a panel of 1000 Genomes Project. The method of ECS is described in our paper(`Paper Link <http://bing.com>`_).

Required options
----------------------

- ``--gene-assoc``
- ``--sum-file [/path/to/summary/file]``
- ``--vcf-ref [/path/to/vcf/file]``
- ``--keep-ref``
- ``[--saved-ref "previous/output/path/prefix"]``
- ``--out [output/path/prefix]``

**See an analysis example at:** `Gene-based association analysis <#gene-based-association-analysis>`_


Explanations and Optional options
--------------------------------------

- ``--gene-assoc``: The main function option.
- ``--sum-file``: The file containing GWAS summary statistics.

    Three columns of the GWAS summary statistic file, chromosome, physical position and *p*-value are minimal requirement. The default column names are CHR, BP and P respectively. Otherwise, users should specify the name by ``--chrom-col``, ``--pos-col`` and ``--p-col`` respectively.

    .. table::
        :align: center

        === ====== ======
        CHR BP     P
        === ====== ======
        1   751756 0.979957
        1   752566 0.863844
        1   752894 0.55814
        1   753405 0.968401
        1   755890 0.918246
        === ====== ======


- ``--vcf-ref``: The file containing genotypes to calculate genotypic correlations.
- ``--keep-ref``: Save the encoded genotypes in VCF for future usage, which will speed up next analysis.
- ``--saved-ref``: Instead of using ``--vcf-ref``, one can directly specify the path of encoded genotypes last time by specifying last output path.
- ``--filter-maf-le``: Filter out variants with minor allele frequency less or equal than the specified value.
- ``--out``: Specify the path and prefix name of the output files. The main output file of the gene-based analysis is ***.gene.pvalue.txt** or ***.gene.pvalue.xls**. The following

    .. csv-table::
        :file: ./table/demo.gene.pvalue.csv
        :header-rows: 1
        :align: center

    columns in the output file are gene symbol, number of variants in the gene, *p*-values of gene-based association test, and the detailed information of the top variant within the gene(i.e., the variant with smallest *p*-value). These columns include chromosome, physical position, *p*-value, whether the top variant was ignored in the gene-based association analysis, and gene feature annotations according to RefGene and GENCODE.

Infer causal genes based on GWAS summary statistics and eQTLs by Mendelian randomization analysis framework for causal gene estimation(MACG)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One can perform multiple IVs based MR analysis to infer casual gene or transcript by an integrative framework named MACG. MACG adopted two multiple IVs based MR methods for causality test and casual effect estimation of a gene’s expression to a phenotype, median-based MR and ML-based MR. MACG needs two major inputs, GWAS and eQTL summary statistics respectively. The GWAS summary statistics refer to the logarithm of odds ratio or regression coefficients and the corresponding standard errors(SEs) from a large-scale GWAS study, indicating the association between IVs and a phenotype. The eQTL summary statistics are similar to that of the GWAS, indicating association between IVs and expression of genes or transcripts in a tissue or cell type. MACG has integrated the pre-calculated cis-eQTLs in 55 tissues or cell-types with gene-level and transcript-level expression from GTEx(version 8).

Required options
---------------------

- ``--macg``
- ``--eqtl-file [path/to/eQTL/file of genes or transcripts]``
- ``--sum-file [/path/to/summary/file]``
- ``--beta-or [y/n]``
- ``--saved-ref  [previous/output/path]``
- ``--out [output/path/prefix]``

**See an analysis example at:** `Gene-based causality analysis <#gene-based-causality-analysis>`_

Explanations and Optional options
---------------------------------------

- ``--macg``: The main function option.
- ``--eqtl-file``: See above description.
- ``--sum-file``: See above description.
- ``--beta-or``: Indicate whether the coefficients(i.e., betas) in the summary statistics file are conventional odds ratios. If yes, KGGSee will automatically transform the betas and SEs by the natural logarithm. 
- ``--saved-ref``: See above description.
- ``--out``: Specify the path and prefix name of the output files. The main output file is the Mendelian randomization analysis results for causal gene estimation, named ***.mr.gene.txt** or ***. gene.mr.gene.xls**. The following

    .. csv-table::
        :file: ./table/demo.mr.gene.csv
        :header-rows: 1
        :align: center

    columns in the output file are gene symbol, number of variants in the gene, *p*-values of causality tests by Median-based MR, detailed causality estimation by Median-based MR, *p*-values of causality tests by maximal likelihood-based MR, detailed causality estimation by maximal likelihood-based MR, chromosome, top GWAS variant position, *p*-value, beta and SE of the top GWAS variant, *p*-value, beta and SE of the top GWAS variant as an eQTL. When a gene has multiple transcripts, the detailed MR results will show MR analysis of all transcripts. Each MR analysis result has four components, the number IVs for the estimation, the estimated causal effect, the standard error of the estimation, and the *p*-values.

Finely map genes and estimate relevant cell-types of a phenotype by the single-cell(or bulk-cell) type and phenotype cross annotation framework(SPA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One can simultaneously prioritize phenotype associated genes and cell-types with GWAS *p*-values and gene/transcript expression profile. The GWAS *p*-values types and expression were analyzed by an iterative prioritization procedure. In the procedure, phenotype-associated genes were prioritized by a conditional gene-based association(using the ECS again) according to the genes’ selective expression in disease related cell-types while the phenotype related cell-types were prioritized by an enrichment analysis of Wilcoxon rank-sum test for phenotype-associated genes’ selective expression. The phenotype-associated gene list and phenotype related cell-type list were updated by turns until the two list were unchanged. The detailed method is described in our paper(`Paper Link <http://bing.com>`_).

Required options
-------------------

- ``--spa``
- ``--expression-file [path/to/expression/file]``
- ``--only-hgnc-gene``
- ``--sum-file [/path/to/summary/file]``
- ``--saved-ref  [previous/output/path/prefix]``
- ``--filter-maf-le 0.02``
- ``--out [output/path/prefix]``

**See an analysis example at:** `Estimate relevant cell-types of a phenotype <#estimate-relevant-cell-types-of-a-phenotype>`_

Explanations and Optional options
----------------------------------

- ``--spa``: The main function option.
- ``--multiple-testing``: The multiple testing method to select significant genes for the conditional analysis. There are three settings. *bonf*: Standard Bonferroni correction given a family-wise error rate specified by ``--p-value-cutoff``.  *benfdr*: Benjamini-Hochberg method to control the false discovery rate. *fixed*: Filtering by a fixed *p*-value cutoff.
- ``--p-value-cutoff``: The cutoff for the multiple testing.
- ``--only-hgnc-gene``: Only consider genes with hgnc gene symbols.
- ``--expression-file``: The path of gene expression file.

    The expression file contains gene symbols(the first column), expression mean and standard errors of the gene or transcript in a cell types or clusters. One can include the Ensembl transcript ID of a gene in the first column. When a gene has multiple transcripts, each row can only contain the data of transcript. The standard error is not pre-requisite.

    .. csv-table::
        :file: ./table/gene.expression.file.csv
        :header-rows: 1
        :align: center

- ``--sum-file``: See above description. 
- ``--filter-maf-le``: See above description.
- ``--out``: Specify the path and prefix name of the output files. One of main output files is the conditional gene-based analysis results, named ***.finemapping.gene.ecs.txt** or ***. finemapping.gene.ecs.xls**. The following

    .. csv-table::
        :file: ./table/demo.finemapping.gene.ecs.csv
        :header-rows: 1
        :align: center

    columns in the output file are gene symbol, chromosome, transcription start position, transcription end position, number of variants in the gene, the LD group ID of genes, *p*-values of gene-based association test, *p*-values of conditional gene-based association test, and the selective expression score in enriched tissue or cell-types.

    Another main output files is the selective expression enrichment analysis results at different tissues or cell types, named ***.celltype.txt** or ***. celltype.xls**. The following

    .. csv-table::
        :file: ./table/demo.celltype.csv
        :header-rows: 1
        :align: center

    columns in the output file are tissue or cell-type names, the *p*-value of enrichment according to the selective expression derived from the robust regression *z*-score, the logarithm of *p*-value.

eQTL-guided gene/transcript-based association
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One can perform gene/transcript-based association analysis guided by eQTLs. The statistical method is still the effective chi-square statistics(ECS). ECS only combines the GWAS *p*-values of eQTLs of a gene or transcript. The pre-calculated cis-eQTLs of gene-level and transcript-level expression in 55 tissues or cell-types from GTEx(version 8) have been integrated into KGGSEE.

Required options
--------------------

- ``--gene-finemapping``
- ``--eqtl-file [path/to/eQTL/file of genes or transcripts]``
- ``--filter-eqtl-p 0.01``
- ``--expression-file [path/to/expression/file]``
- ``--calcu-selectivity``
- ``--sum-file [/path/to/summary/file]``
- ``--filter-maf-le 0.02``
- ``--saved-ref  [previous/output/path]``
- ``--out [output/path/prefix]``

**See an analysis example at:** `eQTL-Gene-based association analysis <#eqtl-guided-gene-based-association-analysis>`_

Explanations and Optional options
-----------------------------------

- ``--gene-finemapping``: The main function option.
- ``--multiple-testing``: See above description.
- ``--p-value-cutoff``: See above description.
- ``--only-hgnc-gene``: See above description.
- ``--gene-score``: See above description.
- ``--calcu-selectivity``: See above description.
- ``--eqtl-file``: The path of eQTL file based on the gene-level or transcript-level expression.

    The eQTL file has a similar fasta format. The first row is just column names and optional. The eQTL data of a gene or transcript start with the symbol “>”. In the same row, the gene symbol, Ensembl transcript ID and chromosome name are included and delimited by tab characters. The subsequent row contains the summary statistics the eQTL for the gene or transcript. The tab-delimited columns are physical position, reference allele, alternative allele, frequency of alternative allele, estimated effect size, standard error of the estimation, *p*-value, effective sample sizes and determination coefficient in a linear regression respectively. In the regression, the number of alternative alleles is used as an independent variable. On KGGSEE, we have pre-calculated the eQTL data using GTEx data(version 8). Variants within 1MB upstream and downstream of a gene or a transcript are included. 

    .. code::

        #symbol id      chr     pos     ref     alt     altfreq beta    se      p       neff    r2
        >DDX11L1	ENST000456328	1						
        13418	G	A	0.161	-0.03	0.013	0.027	62	0.076
        19391	G	A	0.11	0.065	0.027	0.017	63	0.085
        107970	G	A	0.285	-0.024	0.01	0.018	86	0.063
        >MIR6859	ENST0000612080	1						
        13418	G	A	0.161	-0.03	0.013	0.027	62	0.076
        19391	G	A	0.11	0.065	0.027	0.017	63	0.085
        62578	G	A	0.081	0.062	0.024	7.98E-03	67	0.098
        99334	A	G	0.088	0.071	0.035	0.043	56	0.07
        …	…	…	…	…	…	…	…	…
    
- ``--filter-eqtl-p``: Set the *p*-value cutoff to filter out less significant eQTL for the analysis.
- ``--sum-file``: See above description. 
- ``--filter-maf-le``: See above description.
- ``--out``: Specify the path and prefix name of the output files. There are three main result files. One is the gene-based association result file, named ***.gene.txt** or ***.gene.xls**. The following

    .. csv-table::
        :file: ./table/demo.gene.csv
        :header-rows: 1
        :align: center

    columns in the output file are gene symbol, number of variants in the gene, chromosome, gene start position, gene end position, the position of top variant, the *p*-value, coefficient and standard error of the variant for gene expression as an eQTL. The second is the conditional gene-based analysis results, named ***.finemapping.gene.ecs.txt** or ***. finemapping.gene.ecs.xls**. The third is the selective expression enrichment analysis results at different tissues or cell types, named ***.selectivity.ecs.txt** or ***.selectivity.ecs.xls**. Their file formats are the same as above.

Options Index
===============

Inputs/outputs
~~~~~~~~~~~~~~~~

    .. csv-table::
        :file: ./table/input.output.options.index.csv
        :header-rows: 1
        :align: center

Quality control
~~~~~~~~~~~~~~~~~~~

    .. csv-table::
        :file: ./table/quality.control.options.index.csv
        :header-rows: 1
        :align: center

Functions
~~~~~~~~~~~

    .. csv-table::
        :file: ./table/functions.options.index.csv
        :header-rows: 1
        :align: center

Utilities
~~~~~~~~~~~

    .. csv-table::
        :file: ./table/utilities.options.index.csv
        :header-rows: 1
        :align: center
