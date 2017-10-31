# SC RNA-seq analysis piepline 



## Step1: Preprocessing (From Fastq to UMI exprssion matrix)
1) Run Processing_part1_including_intornic_reads.sh

The input are two fastq files (one for each R1 and one for R2).  The output we care about is a tsv file that is a matrix of genes (rows) by cells(columns).  UMIs (transcripts) whose barcodes differ by one basepair are merged together when counting umis/gene 

Please See drop-seq tools for explanation of various functions in the actual script

## Step2: Downstream Analysis using Seurat Package in R & Python (Seaborn) for heatmap visualization

All functions can be found in post_processing.R 

### a) Creating a "seurat object" for downstream analysis and visualizing basic qc metrics 

Usage: seurat_object<- process_umi('nameoffile.txt', 'name_you_want_to_use')

Basic QC metrics include a plot of the distribution of UMIs in all cells and number of genes in all cells

### b) Running clustering and visualization (Using TSNE and KNN clustering) 

Usage: process_umi_step2(seurat_object, k for knn clustering, umi cutoff (to control for dobulets))

The umi cutoff is currently chosen by visualizing the nGene plotted against the nUMI for all cells.  Cells that have an abnormally high number of UMI 
for a given number of genes are considered to be doublets. To estimate k we can use Sqrt(n) where n is the number of cells as a rough estimate.  

### c) Identify cell type specific marker genes and write two tables (list of marker genes for each cluster & the average umi for each cluster)

Usage: identify_clusters(seurat_object, # of clusters, 'name_of_sample')  

Note:The name of sample will be used to name the output files, and the # of cells in each cluster will be outputted to the screen

### d) Create a heatmap of the cell type specific genes  (python script)










