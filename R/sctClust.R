##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author dylanmr
##' @export
sctClust <- function(x, dims=NULL, resolution=NULL) {

  future::plan(future::sequential) 
  
  if(is(x, "SingleCellExperiment")) {
    x <- as.Seurat(x)
  } 
  
  if(is.null(resolution)) {
    resolution <- 0.8
  } else {
    resolution <- resolution
  }
  
  x %>% 
    SCTransform(method = "glmGamPoi") %>% 
    RunPCA() -> x
  
  if(is.null(dims)) {
    dims <- round(as.numeric(maxLikGlobalDimEst(data = x@reductions$pca[, 1:50], k = 20)))
  } else {
    dims <- dims
  }
  
  x %>% 
    RunUMAP(dims = seq(dims)) %>%
    FindNeighbors(dims = seq(dims)) %>%
    FindClusters(resolution = resolution) -> x

  return(x)

}
