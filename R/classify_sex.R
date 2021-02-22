##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author dylanmr
##' @export
classifySex <- function(x) {

  # score each geneset
  x <- AddModuleScore(x, features = c("Uty","Ddx3y"), name = "sex")
  # update name to access scores 
  acc_name <- paste0("sex",1)
  
  # calculate average score across all groups
  tibble(
    score = unlist(x[[acc_name]]), 
    hash_id = x$hash_id
  ) %>%
    group_by(hash_id) %>% 
    summarise(mean=mean(score))  %>%
    mutate(sex_call = if_else(mclust::Mclust(mean, G=2)$classification == 2, true="M", false="F")) %>% 
    dplyr::select(-mean) -> score
  
  x$sex_call <- score$sex_call[match(x[[]]$hash_id, score$hash_id)]
  x[["sex1"]] <- NULL
  x[["sex2"]] <- NULL
  
  return(x)
}
