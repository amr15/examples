SC RNA-seq analysis piepline
============================================================
Time-stamp: "2017-11-01 10:30:16"


# Overall steps

1. Preprocessing (From Fastq to UMI exprssion matrix)
2. Generate Gene type heatmap

# Detailed steps
## Step1: Preprocessing (From Fastq to UMI exprssion matrix)

* **Code**: `Processing_part1_including_intornic_reads.sh`
* **Input**: two fastq files (one for each R1 and one for R2).
* **The output**: 
  *  Genes by cells  matrix of genes (rows) by cells(columns).
  *  UMIs (transcripts) whose barcodes differ by one basepair are merged together when counting umis/gene 
* *Please See drop-seq tools for explanation of various functions in the actual script*

## Step2: Downstream Analysis using Seurat Package in R & Python (Seaborn) for heatmap visualization
All functions can be found in `post_processing.R`

### Creating a "seurat object" for downstream analysis and visualizing basic qc metrics 

* Usage: `seurat_object<- process_umi('nameoffile.txt', 'name_you_want_to_use')`
* Basic QC metrics include a plot of the distribution of UMIs in all cells and number of genes in all cells

### Running clustering and visualization (Using TSNE and KNN clustering) 

* Usage: `process_umi_step2(seurat_object, k for knn clustering, umi cutoff (to control for dobulets))`

The umi cutoff is currently chosen by visualizing the nGene plotted against the nUMI for all cells.  Cells that have an abnormally high number of UMI 
for a given number of genes are considered to be doublets. To estimate k we can use Sqrt(n) where n is the number of cells as a rough estimate.

### Identify cell type specific marker genes 

* Usage: `identify_clusters(seurat_object, # of clusters, 'name_of_sample')`

The name of sample will be used to name the output files, and the number of cells in each cluster will be outputted to the screen.
For each cluster a table with a list of signficantly up or downregulated marker genes as well as a table with the average UMI for each gene will be outputted. 


### d) Create a heatmap of the cell type specific genes  (python script)
* **Input**: a matrix of the average umi for each gene in each cell type (cluster).
*  **outputs**:
   1. A normalized heat map (z-score)
   2. raw unnormalized heatmap
   3. A csv table showing the order of genes in the normalized and unnormalized heatmap 









