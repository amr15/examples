
#### This is for Version 2.0 of Seurat
require(Seurat)
require(tidyverse)
require(data.table)


#1st step in processing the umi file, check out the plots generated here to see how your data look
#function takes in your tsv file, and what you'd like to name your project in quotes ie: "project"
process_umi<-function(file,name){
  tissue.data <- fread(file);  setDF(tissue.data,rownames = tissue.data$V1);  tissue.data$V1 <- NULL
  tissue <- CreateSeuratObject(raw.data = tissue.data,project=name, min.cells=10, min.genes=200, normalization.method = "LogNormalize",do.scale=TRUE)
  mito.genes <- grep("^mt-", rownames(tissue@data), value = T)
  percent.mito <- colSums(as.array(expm1(tissue@data[mito.genes, ]))) / colSums(as.array(expm1(tissue@data)))
  tissue <- AddMetaData(tissue, percent.mito, "percent.mito")
  VlnPlot(tissue, c("nGene", "nUMI", "percent.mito"), nCol = 3)
  GenePlot(tissue, "nUMI", "nGene")
  return (tissue)}

## to generate TSNE Plots (give seurat object, low/high umi estimates for doublet filtering; Density paramter for dbscan clustering)
process_umi_step2 <- function (tissue, density_param,cutoff){
  val<-cutoff  #find_umi_cutoff based on plot of nUMI vs nGenes
  tissue <- SubsetData(tissue, subset.name = "nUMI", accept.high = val)
  tissue <- FindVariableGenes(tissue ,do.plot = FALSE)
  tissue <- RunPCA(tissue,do.print = FALSE, pc.genes = tissue@var.genes)
  tissue=ProjectPCA(tissue,do.print=FALSE)
  tissue <- RunTSNE(tissue, dims.use = 1:11, do.fast = T)
  tissue=FindClusters(tissue,k.param=density_param, dims.use = 1:13)
  TSNEPlot(tissue,do.label = TRUE)
  return (tissue)}

##to create average umi tables with all genes for each clusters
make_umi_table<-function(tissue,cid,outfile)
{
  c<-tissue@raw.data[,rownames(subset(tissue@meta.data, tissue@meta.data$res.0.8==cid))]
  print (length(colnames(c))) # this is the number of cells in the cluster
  write.table(data.frame(rowMeans(c, na.rm=FALSE, dims=1)), sep = "\t",outfile)
}


#to identify the marker genes and write umi tables
identify_clusters <-function(tissue,ncluster,type) {
  tissue@meta.data$orig.ident<-tissue@meta.data$res.0.8
  tissue=SetAllIdent(tissue,id = "orig.ident")
  for (i in (0:ncluster)){
    marker=FindMarkers(tissue,i)
    if (length(colnames(tissue@raw.data[,rownames(subset(tissue@meta.data, tissue@meta.data$res.0.8==i))]))>100){
      marker<-marker[rownames(as.matrix(rowSums(tissue@raw.data[rownames(marker),])) > .05*length(rownames(marker))),]
    }
    else{
      marker<-marker[rownames(as.matrix(rowSums(tissue@raw.data[rownames(marker),])) > 5),]
    }
    write.csv(marker,paste0(type,".markergenes",i,'.csv'))
    make_umi_table(tissue,i,paste0(type,"cluster",i,"averageumi.txt"))}}



# To Plot Genes of your choice in the clusters and clean up your TSNE plot (input is your seurat object and a list of cluster ids that you want to remove from your plot
clean_tsne<-function(tissue, badcluster, genelist){
  TSNEPlot(tissue, do.return = T, no.legend = FALSE,cells.use =rownames(tissue@meta.data[!(tissue@eta.data$res.0.8 %in% badcluster),]))
  FeaturePlot(tissue,genelist,pt.size = 1,col = c('grey','red'))}


###example usage###########
## LI <-process_umi('~/Documents/Ren_lab/sc/LI.counts.tsv','LI') # name of file, what I want to call the seurat object
## LI<-process_umi_step2(LI,5000, 6000, 1) # umi low and high range for doublet limit, density paramater for clustering for TSNE
##identify_clusters<-function(LI,6,"Liver") #of clusters, name for saving
