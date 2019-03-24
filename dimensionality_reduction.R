install.packages('stylo')
library('stylo')
library(sparsesvd)
library(Rtsne)

#' Performs sparse single value decomposition on the scaled UMI values (dimensions of genes x cells))
#' Default rank of matrix is set to 20 
#' Coordinates from running TSNE on the 'v' matrix are saved back into the original Seurat object
#' Right singular vectors saved as new "PCs"  
#' Returns original seruat object with modified TNSE coordinates and "PCs" 

svd<-function(tissue){
  decomp <-sparsesvd(as(tissue@scale.data,'sparseMatrix'), rank=20L)
  tissue@dr$pca@cell.embeddings <- decomp$v
  projection <-Rtsne(decomp$v)
  tissue@dr$tsne@cell.embeddings <- projection$Y
  colnames(tissue@dr$tsne@cell.embeddings)<-c('tSNE_1', 'tSNE_2')
  rownames(tissue@dr$tsne@cell.embeddings)<-colnames(tissue@scale.data)
  colnames(tissue@dr$pca@cell.embeddings)<-c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14',\
                                             'PC15','PC16','PC17','PC18','PC19','PC20')
  rownames(tissue@dr$pca@cell.embeddings)<-colnames(tissue@scale.data)
  return(tissue) } 
