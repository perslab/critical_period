##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param sce_object
##' @return
##' @author dylanmr
##' @export
determine_quantile <- function(seur.hto) {
  range <- seq(.8,1,by=0.01)
  # run htodemux on range of values
  dat <- purrr::map(range, function(x) HTODemux(seur.hto, positive.quantile = x, assay = "RNA", verbose = 0)) %>%
    purrr::map_dfr(., function(x) table(x$hash.ID)) %>% mutate(quant=range)
  # find point at which there are more negatives than dubs
  # later convert this to interpolate
  # https://stackoverflow.com/questions/55875829/how-to-interpolate-a-single-point-where-line-crosses-a-baseline-between-two-poin
  dat %>% 
    mutate(diff = Doublet-Negative) %>% 
    pull(diff) -> val
  idx <- min(which(val < 0))
  # two lines wont always intersect at specified values
  # if they dont re-run demux with new parameters
  if(is.infinite(idx)) {
    range <- seq(.99,1,by=0.001)
    dat <- purrr::map(range, function(x) HTODemux(seur.hto, positive.quantile = x, assay = "RNA", verbose = 0)) %>%
      purrr::map_dfr(., function(x) table(x$hash.ID)) %>% mutate(quant=range)
    dat %>% 
      mutate(diff = Doublet-Negative) %>% 
      pull(diff) -> val
    idx <- min(which(val < 0))
    # if they still dont intersect, return .999 as quantile
    if(is.infinite(idx)) {
      quant <- 0.999
    } else {
      quant <-  dat$quant[idx-1]
    }} else {
    quant <- dat$quant[idx-1]
  }
  return(list(dat, quant))

}

