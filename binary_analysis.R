
install.packages("Rtsne")
library('Rtsne')
install.packages('ade4')
library('ade4')

sparse_analysis<-function(tissue) {
  tissue@raw.data[tissue@raw.data >=1] <-1
  rows<-rownames(as.matrix(which(rowSums(tissue@raw.data)>=10)))
  cols<-rownames(as.matrix(which(colSums(tissue@raw.data)>=200)))
  tissue@scale.data<-tissue@raw.data[rows,cols]
  return(tissue) 
  }

test_metrics<- function(tissue,n){
  distance<-dist.binary(t(sample@scale.data), method=n, diag = TRUE, upper = FALSE)
  new_viz<-Rtsne(distance)
  sample@dr$tsne@cell.embeddings <- new_viz$Y
  colnames(sample@dr$tsne@cell.embeddings)<-c('tSNE_1', 'tSNE_2')
  rownames(sample@dr$tsne@cell.embeddings)<-colnames(sample@scale.data)
  return(sample)  
}

