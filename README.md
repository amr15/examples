
Step1: Preprocessing (From Fastq to UMI exprssion matrix (columns are cells and rows are genes))

Please See drop-seq tools for explanation of various functions

Step2: Downstream Analysis using Seurat Package in R 

a) Creating a "seurat object" for downstream analysis and visualizing basic qc metrics 

Usage: seurat_object<- process_umi('nameoffile.txt', 'name_you_want_to_use')

b) Running clustering and visualization (Using TSNE and KNN clustering) 

Usage: process_umi_step2(seurat_object, k for knn clustering, umi cutoff (to control for dobulets)

note: umi cutoff is currently chosen by visualizing the nGene plotted against the nUMI for all cells.  Cells that have an abnormally high # of UMI 
for a given number of genes are considered to be doublets. 

c) Identify cell type specific marker genes and write two tables (list of marker genes for each cluster & the average umi for each cluster)

Usage: identify_clusters(seurat_object, # of clusters, 'name_of_sample')  

Note:The name of sample will be used to name the output files, and the # of cells in each cluster will be outputted to the screen

d) Create a heatmap of the cell type specific genes 







