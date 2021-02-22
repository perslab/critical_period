##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author dylanmr
##' @export
gen_milo <- function(x, k, d) {
  
  x <- Seurat::as.SingleCellExperiment(x[[1]], assay="SCT")
  milo <- miloR::Milo(x)
  milo <- miloR::buildGraph(milo, k = k, reduced.dim = "PCA", d = d)
  milo <- miloR::makeNhoods(milo, prop = 0.1, k = k, d=d, refined = TRUE, reduced_dims = "PCA")
  milo <- miloR::calcNhoodDistance(milo, d=d)
  milo <- miloR::countCells(milo, meta.data = data.frame(colData(milo)), sample="hash_id")
  milo <- miloR::buildNhoodGraph(milo)
  return(milo)
  
}
