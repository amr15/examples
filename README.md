SC RNA-seq analysis piepline
============================================================
Time-stamp: "2017-11-02 14:14:48"
# Install 

``` shell
export PATH="$PATH:$PWD/Drop-seq_tools-1.12/"
```

# Overall steps

1. Preprocessing (From Fastq to UMI exprssion matrix)
2. Generate Gene type heatmap

# Detailed steps
## Step1: Preprocessing (From Fastq to UMI exprssion matrix)

* **Code**: `Processing_part1_including_intornic_reads.sh`
* **Input**: two fastq files (one for each R1 and one for R2).
* **The output**: 
  *  Genes by cells matrix with UMIs 
  *  UMIs (transcripts) whose barcodes differ by one basepair are merged together when counting umis/gene 
* *Please See drop-seq tools for explanation of various functions in the actual script*

## Step2: Downstream Analysis
All functions can be found in `post_processing.R`. Require: 
  * [Seurat Package](http://satijalab.org/seurat/) in R 
  * [Seaborn](https://seaborn.pydata.org/) for heatmap visualization in python
### prepare "seurat object" and QC
* Usage: `seurat_object<- process_umi('nameoffile.txt', 'name_you_want_to_use')`
* Basic QC metrics include a plot of the distribution of UMIs in all cells and number of genes in all cells
### Clustering and visualization (Using TSNE and KNN clustering) 

* Usage: `process_umi_step2(seurat_object, k for knn clustering, umi cutoff (to control for dobulets))`
  
* Cell selection thresholds: 
  1. `min_reads`: minimal reads (*10000* Habib et al. 2017) 
  2. `min_genes`: minimal number of genes per nucleus/cell (*200* Habib et al. 2017) 
  3. `max_reads`: Cells that have an abnormally high number of UMI for a given number of genes are considered to be doublets.
  4. The umi cutoff is currently chosen by visualizing the nGene plotted against the nUMI for all cells.

* To estimate k we can run `optimize_k(seurat_object)`  to get out a table.  Ideal cluster choice maximizes Dunn & Shilouette Indicies and minimizes the Davies Bouldin index.  

* If you would like to perform dimensionality reduction using SVD instead of PCA, you can use the function in dimensionality_reduction.R to output an updated seurat object.  Perform this before using `optimize_k(seurat_object)`

* Usage: `seurat_object<-svd(seurat_object)`

### Identify cell type specific marker genes 

* Usage: `identify_clusters(seurat_object, # of clusters, 'name_of_sample')

The name of sample will be used to name the output files, and the number of cells in each cluster will be outputted to the screen.
For each cluster a table with a list of signficantly up or downregulated marker genes as well as a table with the average UMI for each gene will be outputted. 


### Create a heatmap of the cell type specific genes  (python script)
* **Input**: a matrix of the average umi for each gene in each cell type (cluster).
*  **outputs**:
   1. A normalized heat map (z-score)
   2. raw unnormalized heatmap
   3. A csv table showing the order of genes in the normalized and unnormalized heatmap 

### Create a Violinplot displaying raw UMI distribution of marker genes  (python script)
* **Input**: 
   1. A matrix rows(marker genes), columns(cells)  filled with raw UMI counts
   2. A matrix rows(cells) and one column with the cluster identity 
*  **outputs**:
   1. Violinplots of raw UMI distribution for 8 marker genes
   You can change the subplot in the python script to accomodate more or less markers genes.  plt.subplot(241) means create a 2x4 grid,  and place this subplot in the 1st position of that grid







