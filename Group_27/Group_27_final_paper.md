# **RNA-Seq Data Analysis**
**Differential Expression Analysis	- DESeq2 -**

- Introduction
- Algorithm
- Tutorial
- Analysis
- Application
- Related
- Reference

## **Introduction**

**Differential expression**

Differential expression is an analysis that determines quantitative changes in expression of genes between individuals and groups of people. It uses statistical testing to determine if the difference in gene expression between varying individuals and groups is significant or by random chance.

**DESeq2**

DESeq2 is an analysis tool used to analyze the count data from RNA-Seq. It is available in R and Bioconductor packages. It is used to determine whether the genes are differentially expressed. The expressions of genes are tested using negative binomial generalized linear models (GLM). DESeq2 uses dispersion and logarithmic fold changes based on count data to determine which genes are differentially expressed. It also detects and corrects dispersion estimates that are too low through modeling of the dependence of the dispersion on the average expression strength over all samples.

## **Algorithm**

**Normalization**

Normalization is the first step of the data analysis process. It is also the most important step. With normalization, the expression levels of genes are more comparable between and within samples. After normalization, a GLM is fitted to each gene to determine the overall expression strength of the gene.
The read count is modeled as a negative binomial distribution with mean and dispersion. The mean is calculated using the equation below:

![equation1](./equation_1.png)

To get the mean *u<sub>ij<sub>*, the quantity *q<sub>ij<sub>* is scaled by the gene-specific normalization factor *s<sub>ij<sub>*, which is also known as the size factor. *q<sub>ij<sub>* represents the gene length, which is proportional to the concentration of cDNA fragments from the gene in the sample. Since *s<sub>ij<sub>* is considered as a constant within a sample, it can be replaced by another constant, *s<sub>j<sub>*, which is a constant that can be used for all genes in a sample. *s<sub>j<sub>* accounts for differences in sequencing depth between samples. The constant is estimated using the median-of-ratios method shown below:

![equation2](./equation_2(1).png)

In the method, counts are divided by the sample-specific normalization factors determined by median ratio of gene counts relative to the geometric mean per gene.

The overall expression strength of the gene is determined using GLMs with a logarithmic link:

![equation3](./equation_3(1).png)

The design matrix element *x<sub>jr<sub>* is an element that indicates whether a sample *j* is treated or not. The GLM fit returns the coefficient *𝜷<sub>ir<sub>*, which indicates the overall expression strength of the gene. The GLMs also return the log2 fold change between treatment and control samples.

## **Tutorial** (Coming From FeatureCounts)

**Input Files:**
The input files for Deseq two are a count matrix and column data. The count matrix comes from featureCounts where the ith row and jth column describe the expression level for gene i in sample j. The counts must be raw in order for the statistical model for Deseq2 to hold true. The column data is a table with metadata on the count matrix columns. 

![figure1](./Count_Matrix1.png)

Figure 1: Example of a count matrix where the ith row and jth column describe the expression level for gene i in sample j. 

![figure2](./Column_Data.png)

Figure 2: Example of column data, where the first column is the sample name, the second column is the patient type, the third column is the treatment type, and the fourth column is the amount of time for the treatment. The columns tell the metadata for each sample in a count matrix.

**Creating the DeseqDataSet (DDS) object:**
The DeseqDataSet object is the parameter for Deseq2. Since the count matrix comes from FeatureCounts, the DeseqDataSet object is created using the DeSeqDataSetFromMatrix command. This command takes in three parameters. The first is the count matrix. The second is the column data, and the third is the Design Condition. The design condition explains how the counts for each gene depend on the variables in column data. The command to create the DeseqDataSet object is: Dds <- DESeqDataSetFromMatrix (count matrix, column data, Design Condition).

**Running Deseq2:**
The Deseq2 command takes in a single parameter, which is the DeseqDataSet object. The command to run Deseq2 is: Dds <- DESeq(Dds). 


## **Analysis**

Results Table:
After running Deseq2, running this command: Results <- results(Dds), gives the results table. 

![figure3](./Results_Table1.png)

Figure 3: This is an example results table which shows you base mean (expression), log2FoldChange, lfcSe, Stat, and p-value for each gene.

The three most important results in the results table for differential expression are expression levels, log2Foldchange, and p-Values. The Log2FoldChange shows how much the gene’s expression seems to have changed due to treatment from sample to sample. Log2FoldChange is reported on a logarithmic base 2 scale. This means that a log2FoldChange value of 1.5 indicates that the gene’s expression increased by 2^1.5 = 2.82. The p-values show the confidence of the log2FoldChange value and indicate whether or not there is a differential expression. The p-value represents the probability that the results happened by random chance. 

**MA Plot:**

![figure4](./MA_Plot.png)

Figure 4: MA plot where the x-axis represents average expression and y-axis represents  log2FoldChange

The first common way to visualize the results from the results table is through an MA plot. In an MA plot, each gene is a point, and each point is colored if the p-value for the log2change of that specific gene is less than 0.1. Colored points indicate that there is enough evidence to conclude that the gene is differentially expressed. The x-axis represents the average expression for the gene, while the y-axis represents the Log2FoldChange for each gene. 

**Volcano Plot:**

![figure5](./Volcano_Plot1.png)

Figure 5: Volcano Plot where the x-axis represents log2FoldChange and y-axis represents the significance values


The second common way to visualize the results from the results table is through a volcano plot. In a volcano plot, each point represents a gene, and a point is colored if the log2FoldChange meets a minimum threshold for p-value. A colored point represents a gene that is differentially expressed. In the Volcano plot, the x-axis represents log2Foldchange and the y-axis represents significance values. The right half of the plot represents genes that are upregulated, while the left half represents genes that are downregulated. The top of the plot represents the most statistically significant genes. 


## **Application**

The following are three ways DESeq2 can be applied and an example of such an application in a published experiment.


#### Differentiation


##### Differentiate between two genes or more genes:

Human gut microbes impact host serum metabolome and insulin sensitivity (H. Pedersen, et al.  Nature 2016.)



*   “For associations of microbial functional modules with metabolic syndrome and type 2 diabetes, the ranks were based on Wald statistics for testing differentially abundant [KEGG Orthology gene groups] with a negative binomial test, using the DESeq2 R package” <sup>[3]</sup>


#### Identification


##### Identify sample expressing differently in drug trial or disease state.

Tumor Evolution and Drug Response in Patient-Derived Organoid Models of Bladder Cancer (Suk Hyung Lee, et al. Cell 2018)



*   “Differentially expressed genes in organoid samples were identified compared to their primary tumors using the DESeq2 package” <sup>[4]</sup>


#### Normalization


##### Normalize expression levels for further computation (e.g. clustering).

An Integrated Genome-wide CRISPRa Approach to Functionalize lncRNAs in Drug Resistance (Assaf C. Bester, et al. Cell 2018)



*   “Transcript counts per gene were normalized by DESeq2, and genes with arithmetic means lower than 0.5 were filtered from further analysis” <sup>[5]</sup>
*   “Normalized counts were subjected to the varianceStabilizingTransformation (VST) function from DESeq2 prior to downstream analysis.” <sup>[5]</sup>

  


## **DESeq2 vs. Other Tools**


##### How does DESeq2 compare against other differential expression analysis tools?

In Schurch et al. 2015<sup>[1]</sup>, “How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?” Eleven tools were evaluated on a range of expression data that differ by the number of biological replicates included. Each tool was also evaluated on a dataset of clean replicates to establish a ‘truth value’ relative to each tool. The control and experimental values were compared in the end to derive truth and false positive rates.


Biological replicates are samples produced with the same experimental design but from different specimens. Technical (clean) replicates are repeated measurements of the same specimen. <sup>[2]</sup>

The researchers concluded that for experiments with >20 biological replicates, all eleven tools performed equally well. However, for experiments with &lt;12 biological replicates DESeq2 (and edgeR) outperformed the rest in minimizing false positive rates while still obtaining relatively high true positive rates.


![schurch](./schurch.jpg)

Schurch et al., RNA Society, 2015


## **References**

[1] [https://pubmed.ncbi.nlm.nih.gov/27022035/](https://pubmed.ncbi.nlm.nih.gov/27022035/) 

[2] [www.nigms.nih.gov/training/documents/module3-biological-and-technical-replicates.pdf](https://www.nigms.nih.gov/training/documents/module3-biological-and-technical-replicates.pdf) 

[3] [https://www.nature.com/articles/nature18646](https://www.nature.com/articles/nature18646) 

[4] [https://pubmed.ncbi.nlm.nih.gov/29625057/](https://pubmed.ncbi.nlm.nih.gov/29625057/) 

[5] [https://www.cell.com/cell/pdf/S0092-8674(18)30384-2.pdf](https://www.cell.com/cell/pdf/S0092-8674(18)30384-2.pdf) 

[6] [https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

[7] [http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

[8] [https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf](https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf)

[9] [https://www.genomatix.de/online_help/help_regionminer/DESeq2.pdf](https://www.genomatix.de/online_help/help_regionminer/DESeq2.pdf)

