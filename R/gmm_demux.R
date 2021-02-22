##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param sce_obj
##' @return
##' @author dylanmr
##' @export
gmm_demux <- function(x) {
  
  as.data.frame(x@assays$RNA@data) %>% 
    rownames_to_column("hash") %>% 
    pivot_longer(-hash) -> hash_data
  
  hash_data %>% 
    group_by(hash) %>% 
    group_map(~ Mclust(.x %>% pull(value), G=2, model="V")$classification) %>% 
    bind_cols() %>% 
    magrittr::set_colnames(unique(hash_data$hash)) %>% 
    mutate(row_tot = rowSums(.)) %>%
    transmute(gmm_id = case_when(row_tot > ncol(.) + 1 ~ "Doublet",
                             row_tot == ncol(.) + 1 ~ apply(.[,c(1:length(unique(hash_data$hash)))], 1, 
                                                            function (x) colnames(.)[which.max(x)]),
                             row_tot == ncol(.) ~ "Negative"),
              gmm_global = case_when(gmm_id == "Doublet" ~ "Doublet",
                                     gmm_id == "Negative" ~ "Negative",
                                  T~ "Singlet")
           
    ) -> res
  
  x[["gmm_id"]] <- res$gmm_id
  x[["gmm_global"]] <- res$gmm_global
  
  return(x)
}
