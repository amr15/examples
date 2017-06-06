install.packages("devtools")
library(devtools)
install_url("https://github.com/satijalab/seurat/releases/download/v1.4.0/Seurat_1.4.0.8.tgz", binary = TRUE)
library(Seurat)
library(dplyr)
library(Matrix)

#1st step in processing the umi file, check out the plots generated here to see how your data look
#function takes in your tsv file, and what you'd like to name your project in quotes ie: "project"
process_umi<- function (file,name){ #taking min cells and max cells to be 1000 and Scaling UMI
 tissue.data = read.table(file, sep = '\t', header=TRUE, row.names=1)
 tissue = new("seurat", raw.data= tissue.data)
 tissue <- Setup(tissue, min.cells = 10, min.genes = 1000, do.scale = T,do.center=F, project = name)
 #qc/sanity plots
 mito.genes <- grep("^mt-", rownames(tissue@data), value = T)
 percent.mito <- colSums(as.array(expm1(tissue@data[mito.genes, ]))) / colSums(as.array(expm1(tissue@data)))
 tissue <- AddMetaData(tissue, percent.mito, "percent.mito")
 VlnPlot(tissue, c("nGene", "nUMI", "percent.mito"), nCol = 3)
 GenePlot(tissue, "nUMI", "nGene")
 return (tissue)}

## To filter out doublets (give seurat object, and low/high umi estimates)
find_umi_cutoff <-function (tissue,low,high){
  expected<- .98* length(tissue@data.info$nGene)
  for (i in low:high){
    if (length((SubsetData(tissue, subset.name = "nUMI", accept.high = i))@data.info$nGene) > expected){
      return(i-1)}else {
        return ("too low")}
  }
}

## to generate TSNE Plots (give seurat object, low/high umi estimates for doublet filtering; Density paramter for dbscan clustering)
process_umi_step2 <- function (tissue, low,high,density_param){
  val<-find_umi_cutoff(tissue, low,high)
  tissue <- SubsetData(tissue, subset.name = "nUMI", accept.high = val) 
  tissue <- MeanVarPlot(tissue ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
  tissue <- PCA(tissue,do.print = TRUE, pcs.print = 5, genes.print = 5)
  tissue=project.pca(tissue,do.print=FALSE)
  tissue <- RunTSNE(tissue, dims.use = 1:11, do.fast = T)
  tissue=DBClustDimension(tissue,1,2,reduction.use = "tsne",G.use =density_param,set.ident = TRUE)
  TSNEPlot(tissue,do.label = TRUE)
  return (tissue)}

##to create average umi tables with all genes for each clusters 
make_umi_table<-function(tissue,cid,outfile)
  {
    c<-tissue@raw.data[,rownames(subset(tissue@data.info, tissue@data.info$DBclust.ident==cid))]
    print (length(colnames(c))) # this is the number of cells in the cluster 
    write.table(data.frame(rowMeans(c, na.rm=FALSE, dims=1)), sep = "\t",outfile) 
  }
  

#to identify the marker genes and write umi tables 
identify_clusters <-function(tissue,ncluster,type) {
  tissue@data.info$orig.ident<-tissue@data.info$DBclust.ident
  tissue=SetAllIdent(tissue,id = "orig.ident")
  for (i in (2:ncluster)){
    marker=FindMarkers(tissue,i)
    #only keep marker genes expressed in 5 % of the cells 
    marker[rownames(as.matrix(rowSums(tissue@raw.data[rownames(marker),])) > .05*length(rownames(marker))),]
    write.csv(marker,paste0(type,".markergenes",i,'.csv'))
    make_umi_table(tissue,i,paste0(type,"cluster",i,"averageumi.txt"))}}
    

###example usage###########
## LI <-process_umi('~/Documents/Ren_lab/sc/LI.counts.tsv','LI') # name of file, what I want to call the seurat object
## LI<-process_umi_step2(LI,5000, 6000, 1) # umi low and high range for doublet limit, density paramater for clustering for TSNE 
##identify_clusters<-function(LI,6,"Liver") #of clusters, name for saving 
