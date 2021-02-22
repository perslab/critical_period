##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param neur_cluster
##' @param huisref
##' @param label
##' @return
##' @author dylanmr
##' @export


integrate_data <- function(query, ref, label) {

  plan(sequential)
  features <- SelectIntegrationFeatures(object.list = list(ref, query), nfeatures = 3000)
  integ.list <- PrepSCTIntegration(object.list = list(ref, query), anchor.features = features, 
                                   verbose = FALSE)
  anchors <- FindIntegrationAnchors(object.list = integ.list, normalization.method = "SCT", 
                                    anchor.features = features)
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", 
                              verbose = FALSE)
  integrated[["proj"]] <- ifelse(is.na(integrated$hash_id), yes = label, no = "cp")
  integrated %>% 
    sctClust(dims=40) -> integrated
  return(integrated)

}
