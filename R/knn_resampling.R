##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param g igraph object generated using the findKNN function
##' @param k number of nearest neighbors sought in the KNN graph
##' @param meta metadata from the single cell object
##' @param perm number of permutations to generat
##' @return
##' @author dylanmr
##' @export
knn_resampling <- function(g, k, meta, perm = 1000) {

  plan(multicore, workers=10)
  # permute KNN indices
  perm_dat <- furrr::future_map(seq(perm), function(x){
    # set random seed
    set.seed(x)
    # generate permuted vector
    perm_vec <- sample(as.vector(t(g$index)), replace = T)
    # obtain matched data for permuted vector
    data.frame(
      cell = rep(rownames(meta), each = k),
      idx = rep(seq(nrow(meta)), each = k),
      hash_neighb = meta$hash_id[perm_vec],
      age_neighb = meta$age[perm_vec],
      diet_neighb = meta$diet[perm_vec],
      geno_neighb = meta$geno[perm_vec],
      pred_neighb = meta$SCT_snn_res.0.8[perm_vec]
    ) %>%
      mutate(
        hash_idx = meta$hash_id[idx],
        age_idx = meta$age[idx],
        diet_idx = meta$diet[idx],
        geno_idx = meta$geno[idx],
        pred_idx = meta$SCT_snn_res.0.8[idx]
      ) %>% 
      # remove any neighbor which is from the same sample
      #filter(hash_neighb!=hash_idx) %>% 
      group_by(idx, geno_neighb) %>% 
      dplyr::slice(1)
  }, .progress = T) %>% 
    data.table::rbindlist(., idcol = "id")
  
  return(perm_dat)

}
