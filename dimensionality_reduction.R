install.packages('stylo')
lobrary('stylo')
library(sparsesvd)
library(Rtsne)

svd<-function(tissue){
  decomp <-sparsesvd(as(tissue@scale.data,'sparseMatrix'), rank=20L)
  projection <-Rtsne(decomp$v)
  tissue@dr$tsne@cell.embeddings <- projection$Y
  colnames(tissue@dr$tsne@cell.embeddings)<-c('tSNE_1', 'tSNE_2')
  rownames(tissue@dr$tsne@cell.embeddings)<-colnames(tissue@scale.data)
  colnames(tissue@dr$pca@cell.embeddings)<-c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','PC16','PC17','PC18','PC19','PC20')
  rownames(tissue@dr$pca@cell.embeddings)<-colnames(tissue@scale.data)
  return(tissue) } 
