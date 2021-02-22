##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param sce_objects single cell experiment objects with HTO data to demultiplex
##' @param range range of quantiles you want to test
##' @return
##' @author dylanmr
##' @export
autothreshold <- function(x, range = seq(0.3, 1, by = 0.001)) {

  seur.hto <- as.Seurat(altExp(x, "hto"), counts = "counts", data = NULL)  
  plan(sequential)
  seur.hto <- NormalizeData(seur.hto,normalization.method = "CLR", block.size=1)
  
  plan(multicore, workers=30)
  # test different thresholds
  df <- 
    furrr::future_map(range, function(x) HTODemux(seur.hto, positive.quantile = x, verbose = 0, assay="RNA"), .progress=T) %>% 
    setNames(paste0("HTOquant_", range))
  
  dat <-
    df %>% 
    purrr::map_dfr(~.x$RNA_classification.global, .progress=T) %>% 
    mutate(quant = range)%>% 
    pivot_longer(-quant) %>% 
    group_by(quant, value) %>% 
    summarize(value = n()) %>% 
    filter(quant != 1) %>% 
    mutate(class = c("Doublet", "Negative", "Sing")) %>% 
    ungroup() %>% 
    group_by(class) %>% 
    arrange(class, quant) %>% 
    mutate(lag = lag(value),
           diff = value-lag) 
  
  # fit regression line
  loess_reg <-
    dat %>% 
    group_by(class) %>% 
    group_map(~ stats::loess(diff ~ quant, data = .x, span=0.1)$fitted) %>% 
    setNames(c("Doublet", "Negative", "Singlet")) %>% 
    bind_cols() %>% 
    mutate(quant = as.numeric(unique(dat$quant))[-length(unique(dat$quant))])
  
  # determine cquantile as minimal distance between DOublets and negatives after 
  # maximum doublet loss
  lower_bound <- filter(loess_reg, quant > quant[which.min(loess_reg$Doublet)])
  quant <-lower_bound$quant[which.min(abs(lower_bound$Doublet-lower_bound$Negative))]
  print(quant)
  
  # optimal quant
  opt <- df[[grep(paste0(quant,"$"), names(df))]]
  #default quant
  def <- df[["HTOquant_0.99"]]
  
  x$optim_class <- opt$RNA_classification.global
  x$optim_id <- opt$hash.ID
  x$def_class <- def$RNA_classification.global
  x$def_id <- def$hash.ID
  
  return(x)

}
