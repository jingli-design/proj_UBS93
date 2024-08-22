

library(GSVA)
library(dplyr)

####function:breast.cancer.predictor####
breast.cancer.predictor <- function(expr.of.sample,expr.of.centroid,marker.gene,HER2.amp.signature.genes,ssgsea.cutoff){
  marker.gene           <- intersect(marker.gene,rownames(expr.of.centroid))
  marker.gene           <- intersect(marker.gene,rownames(expr.of.sample))
  expr.of.sample        <- expr.of.sample[marker.gene,,drop=F]
  expr.of.centroid      <- expr.of.centroid[marker.gene,]
  cor.matrix            <- cor(expr.of.sample,expr.of.centroid,method = "spearman")
  n_samples             <- ncol(expr.of.sample)
  n_subtypes            <- ncol(expr.of.centroid)
  result                <- matrix(NA,nrow = n_samples,ncol = 1)
  for(i in 1:n_samples){
    max_cor                 <- max(cor.matrix[i,])
    max_cor_index           <- which.max(cor.matrix[i,])
    if(is.na(max_cor)){
      result[i,]            <- "cor_is_missing"
    }else{
      result[i,]            <- colnames(cor.matrix)[max_cor_index]
    }}
  result                    <- as.data.frame(result)
  rownames(result)          <- rownames(cor.matrix)
  colnames(result)[1]       <- "subtype"
  HER2.amp.signature.genes               <- HER2.amp.signature.genes
  geneset                                <- list(HER2.amp.signature.genes=HER2.amp.signature.genes)
  sample.res                             <- gsva(expr.of.sample,
                                                 geneset, method="ssgsea",
                                                 verbose=FALSE,ssgsea.norm=FALSE) 
  sample.res                             <- t(sample.res)%>% as.data.frame()
  ssgsea.cutoff                          <- ssgsea.cutoff
  final.subtype                          <- matrix(NA,nrow = n_samples,ncol = 1)
  for (i in 1:n_samples) {
    if(result[i,] == "H_or_L"){
      ssgsea.score                   <- sample.res[i,]
      if(ssgsea.score > ssgsea.cutoff){
        final.subtype[i,]            <- "HER2_amp"
      }else{
        final.subtype[i,]            <- "Luminal"
      }}else{
        final.subtype[i,]            <- result[i,]
      }}
  final.subtype                      <- as.data.frame(final.subtype)
  rownames(final.subtype)            <- rownames(cor.matrix)
  colnames(final.subtype)[1]         <- "subtype"
  predicted.res                      <- list(cor.matrix=cor.matrix,subtype=final.subtype,profiles=expr.of.sample)  
}

