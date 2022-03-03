##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param sc_obj
##' @return
##' @author dylanmr
##' @export

mapCamp <- function(query, dims=40, class = NULL) {
  
  camp <- scRNAseq::CampbellBrainData()

  neurlabs <- 
    read.table(here::here("data/campbell_meta.txt"), header=T) %>% 
    distinct(NeuronsOnlyClusters, .keep_all=T) %>% 
    filter(NeuronsOnlyClusters != "Non-neuronal") %>% 
    mutate(cluster_labs = stringr::str_split_fixed(NeuronsOnlyClusters, pattern = "[.]", n = 2)[,2],
           numeric_labs = stringr::str_split_fixed(NeuronsOnlyClusters, pattern = "[.]", n = 2)[,1])
  
  glialabs <- 
    read.table(here::here("data/campbell_meta.txt"), header=T)  %>% 
    distinct(AllCellSubclusters, .keep_all=T) %>% 
    filter(NeuronsOnlyClusters == "Non-neuronal") %>% 
    mutate(cluster_labs = stringr::str_split_fixed(AllCellSubclusters, pattern = "[.]", n = 2)[,2],
           numeric_labs = stringr::str_split_fixed(AllCellSubclusters, pattern = "[.]", n = 2)[,1])
  
  if(is.null(class)) {
    refdata <- camp$clust_all_neurons
  } else if(class == "neuron") {
    camp <- camp[,camp$clust_neurons%in%neurlabs$numeric_labs]
    rename <- neurlabs$numeric_labs
    names(rename) <- neurlabs$cluster_labs
    refdata <- fct_recode(camp$clust_neurons, !!!rename)
  } else if(class == "glia"){
    camp <- camp[,camp$clust_all_micro%in%glialabs$numeric_labs]
    rename <- glialabs$numeric_labs
    names(rename) <- glialabs$cluster_labs
    refdata <- fct_recode(camp$clust_all_micro, !!!rename)
  }

  future::plan(future::sequential)
  
  camp <- CreateSeuratObject(counts = counts(camp)) %>% 
    SCTransform(method="qpoisson", vst.flavor = "v2")
  anch <- FindTransferAnchors(reference = camp, query = query, normalization.method = "SCT", recompute.residuals = F)
  predictions <- TransferData(anchorset = anch, refdata = refdata) %>% 
    dplyr::select(predicted.id, prediction.score.max)
  query <- AddMetaData(query, metadata = predictions)
  print("done")
  return(query)

}
