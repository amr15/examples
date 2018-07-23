
#### This is for Version 2.0 of Seurat
require(Seurat)
require(tidyverse)
require(data.table)
require(clusterCrit)

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

##get metrics (median umi, # of genes, # of umi/gene for all cells) and save histogam plots of these distributions
metrics<-function(tissue,name){
  print (median(tissue@meta.data$nGene))
  print (median(tissue@meta.data$nUMI))
  print (median(tissue@meta.data$nUMI/tissue@meta.data$nGene))
  pdf(paste0(name,'nGene'))
  hist(tissue@meta.data$nGene)
  dev.off()
  pdf(paste0(name,'nUMI'))
  hist(tissue@meta.data$nUMI)
  dev.off()
  pdf(paste0(name,'Coverage'))
  hist(tissue@meta.data$nUMI/tissue@meta.data$nGene)
  dev.off()
}


## to generate TSNE Plots (give seurat object, low/high umi estimates for doublet filtering; Density paramter for dbscan clustering)
process_umi_step2 <- function (tissue, density_param,cutoff){
  val<-cutoff  #find_umi_cutoff based on plot of nUMI vs nGenes
  tissue <- SubsetData(tissue, subset.name = "percent.mito", accept.high = 0.05)
  tissue <- SubsetData(tissue, subset.name = "nUMI", accept.high = val) 
  tissue <- SubsetData(tissue, subset.name = "nUMI", accept.low = 200)
  tissue <- FindVariableGenes(tissue ,do.plot = FALSE)
  tissue <- RunPCA(tissue,do.print = FALSE, pc.genes = tissue@var.genes)
  tissue=ProjectPCA(tissue,do.print=FALSE)
  tissue <- RunTSNE(tissue, dims.use = 1:11, do.fast = T)
  tissue=FindClusters(tissue,k.param=density_param, dims.use = 1:11)
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
      marker<-marker[rownames(as.matrix(rowSums(tissue@raw.data[rownames(marker),])) > .05*length(rownames(tissue@meta.data))),]
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

#To calculate clustering metrics to determine the optimal number of clusters, returns a table you can use to optimize choice of K
optimize_k <-function(tissue){
  i = 5
  tissue=FindClusters(tissue,k.param=i, dims.use = 1:10)
  data<-tissue@dr$tsne@cell.embeddings
  partitions<-as.integer(tissue@meta.data[rownames(tissue@dr$tsne@cell.embeddings),'res.0.8'])
  indicies<-intCriteria(data,partitions,c("dav","Dunn","Silhouette"))
  for (i in c(10,15,20,25,30,35,40,45,50,55,60,65,70)){ 
    tissue=FindClusters(tissue,k.param=i, dims.use = 1:10)
    data<-tissue@dr$tsne@cell.embeddings
    partitions<-as.integer(tissue@meta.data[rownames(tissue@dr$tsne@cell.embeddings),'res.0.8'])
    index_new<-intCriteria(data,partitions,c("dav","Dunn","Silhouette"))
    indicies<-rbind(data.frame(index_new), data.frame(indicies))}
  write.csv(indicies, paste0(args[1],"koptimal.csv"))
  return(indicies)
  }

### Combine data sets
combine_data_sets<-function(rep1,rep2,protocol1,protocol2){
  rep1 <- FindVariableGenes(rep1 ,do.plot = FALSE)
  rep2 <- FindVariableGenes(rep2 ,do.plot = FALSE)
  hvg.rep1 <- rownames(head(rep1@hvg.info, 1000))
  hvg.rep2 <- rownames(head(rep2@hvg.info, 1000))
  hvg.union <- union(hvg.rep1, hvg.rep2)
  rep1@meta.data[, "protocol"] <- protocol1
  rep2@meta.data[, "protocol"] <- protocol2
  reps <- RunCCA(rep1, rep2, genes.use = hvg.union,add.cell.id1=protocol1 ,group.by = 'protocol',num.cc = 10,scale.data = TRUE)
  reps <- CalcVarExpRatio(reps, reduction.type = "pca", grouping.var = "protocol")
  reps <- SubsetData(reps, subset.name = "var.ratio.pca", accept.low = 0.5)
  reps <- AlignSubspace(reps, reduction.type = "cca", grouping.var = "protocol", dims.align = 1:10)
  reps <- RunTSNE(reps, reduction.use = "cca.aligned", dims.use = 1:10, do.fast = T)
  reps <- FindClusters(reps, reduction.type = "cca.aligned", dims.use = 1:10, save.SNN = T)
  return (reps)}

###example usage###########
## LI <-process_umi('~/Documents/Ren_lab/sc/LI.counts.tsv','LI') # name of file, what I want to call the seurat object
## LI<-process_umi_step2(LI,5000, 6000, 1) # umi low and high range for doublet limit, density paramater for clustering for TSNE
##identify_clusters<-function(LI,6,"Liver") #of clusters, name for saving
