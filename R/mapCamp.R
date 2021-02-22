##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param sc_obj
##' @return
##' @author dylanmr
##' @export
mapCamp <- function(query, dims=NULL, class = NULL) {
  
  if(is.null(dims)) {
    dims <- 40
  } else {
    dims <- dims
  }
  
  camp <- scRNAseq::CampbellBrainData()
  camp <- scater::logNormCounts(camp)
  
  if(is.null(class)) {
    refdata <- camp$clust_all_neurons
  } else if(class == "neur") {
    camp <- camp[,camp$clust_neurons!="miss"]
    refdata <- camp$clust_neurons
  } else {
    refdata <- camp$clust_all_neurons
  }

  camp <- as.Seurat(camp) %>% 
    SCTransform(method="glmGamPoi")
  anch <- FindTransferAnchors(reference = camp, query = query, normalization.method = "SCT", dims = seq(dims))
  predictions <- TransferData(anchorset = anch, refdata = refdata) %>% 
    dplyr::select(predicted.id, prediction.score.max)
  query <- AddMetaData(query, metadata = predictions)
  return(query)

}
