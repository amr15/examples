
install.packages("Rtsne")
library('Rtsne')
install.packages('ade4')
library('ade4')

#' Turns RNA-seq datda into a binary matrix of 1(expressed) 0 (not expressed) 
#' Input is a seurat object and output is a seurat object with expression values that have been binarized 
#' Retains genes expressed in at least 10 cells and cells that express at least 200 genes 

binarize_matrix <-function(tissue) {
  tissue@raw.data[tissue@raw.data >=1] <-1
  rows<-rownames(as.matrix(which(rowSums(tissue@raw.data)>=10)))
  cols<-rownames(as.matrix(which(colSums(tissue@raw.data)>=200)))
  tissue@scale.data<-tissue@raw.data[rows,cols]
  return(tissue) 
  }

#' Performs TSNE on a scaled distance matrix created from UMI data; goal is to test different distance metrics for visualization
#' Input in a seurat object (tissue) and output is a suerat object with new TSNE coordinates (based on distance metric of choice)
#' n can range from 1 to 10, and corresponds to distance metrics here: https://www.rdocumentation.org/packages/ade4/versions/1.7-6/topics/dist.binary

test_metrics<- function(tissue,n){
  distance<-dist.binary(t(sample@scale.data), method=n, diag = TRUE, upper = FALSE)
  new_viz<-Rtsne(distance)
  sample@dr$tsne@cell.embeddings <- new_viz$Y
  colnames(sample@dr$tsne@cell.embeddings)<-c('tSNE_1', 'tSNE_2')
  rownames(sample@dr$tsne@cell.embeddings)<-colnames(sample@scale.data)
  return(sample)  
}

