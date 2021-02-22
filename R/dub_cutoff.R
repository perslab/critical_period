##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param sce_objects
##' @return
##' @author dylanmr
##' @export
dub_cutoff <- function(x) {

  # calculate number of neighbors at each proportion that are doublets
  data.frame("prop"=x$proportion_dub_neighbors) %>% 
    group_by(prop) %>% 
    summarize(n=n()) %>% 
    mutate(pct = n/sum(n)) -> data 
  # find point at which we gain very few doublets as proportion increases 
  cut <- data$prop[PCAtools::findElbowPoint(variance = sort(data$n, decreasing = T))+1]
  vec <- if_else(x$proportion_dub_neighbors <= cut, F, T)
  return(vec)
  
}
