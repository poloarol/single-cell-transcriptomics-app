if (!require('pacman')) install.packages("pacman")

pacman::p_load(DoubletFinder, tidyverse)


pk <- function(data){
  data1 <- paramSweep_v3(data, PCs = 1:10, sct = FALSE)
  data2 <- summarizeSweep(data1, GT = FALSE)
  pk <- find.pK(data2)
  return(pk)
}

pk_ground_truth <- function(data){
  data1 <- paramSweep_v3(data, PCs = 1:10, sct = FALSE)
  ## GT is a vector containing "Singlet" and "Doublet" 
  ## calls recorded using sample multiplexing classification
  ## and/or in silico geneotyping results 
  gt.calls <- data@meta.data[rownames(data1[[1]]), "GT"]
  data2 <- summarizeSweep(data1, GT = TRUE, GT.calls = gt.calls)
  pk <- find.pK(data2)
  return(pk)
}

homotypic <- function(data, annotation, pct = 0.075){
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(pct*nrow(data@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  return(nExp_poi.adj)
}


run_doubletfinder <- function(data, pN, pk, nExp){
  results <- doubletFinder_v3(data, PCs = 1:10, 
                              pN = pN, pK = pK, nExp = nExpi, 
                              reuse.pANN = FALSE, sct = FALSE)
  #seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
}




