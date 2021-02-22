##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param sce_object
##' @param dims reduced dimensions to build graph with
##' @param k number of neighbors to build SNN graph with
##' @param weight defaults to jaccard distance (same as seurat); can use rank or numeric as well
##' @param clus_method defaults to louvain (same as seurat clustering); can use infomap, or walktrap
##' @return returns a list with 3 values. this includes assigned clusters, the cell-cell graph, and the bootstrapping results
##' @author dylanmr
##' @export

cluster_across_params <- function(sce_object, dims, k, weight, clus_method) {

  mat <- reducedDim(sce_object, "PCA")[,1:dims]
  orig <- .cluster_fun(mat = mat, k = k, weight = weight, clus_method = clus_method)
  ratios <- bootstrapStability(mat, FUN= .cluster_fun, clusters=orig$clusters, k=k, weight=weight, boot=1, clus_method=clus_method)
  return(list(clusters = orig$clusters, graph = orig$graph, boot_res = ratios))

}


.cluster_fun <- function(mat, dims, k, weight, boot=NULL, clus_method=louvain) {
  
  g <- bluster::makeSNNGraph(mat, k=k, type=weight)
  
  if(clus_method=="louvain") {
    clusts <- igraph::cluster_louvain(g)$membership
  } else if(clus_method == "walktrap") {
    clusts <- igraph::cluster_walktrap(g)$membership
  } else if(clus_method == "info") {
    clusts <- igraph::cluster_infomap(g)$membership
  }
  
  if(is.null(boot)) {
    return(list(clusters = clusts, graph = g))
  } else {
    return(clusts)
  } 

}

