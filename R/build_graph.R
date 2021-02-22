##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param seur
##' @param dims
##' @param k
##' @param reduction
##' @return
##' @author dylanmr
##' @export
build_graph <- function(seur, dims = 50, k = 30, reduction = "pca") {

  df <- Embeddings(seur, reduction = reduction)[,c(1:dims)]
  g <- findKNN(df, k = k, get.distance=T, get.index=T)
  return(g)

}
