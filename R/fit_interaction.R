##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @return
##' @author dylanmr
##' @export

# model expression of all detectable genes in agrp neurons over age. 
# Use mgcv to speed up calculations in large data-sets. We are using 
# an ordered factor to calculate "difference smooths". As seen here:
# https://fromthebottomoftheheap.net/2017/12/14/difference-splines-ii/
# include a random intercept for cluster litter, hash, and hashed pool. 
# Sex is included as a by factor - need to check if this is the right way to include this
# Include a random intercept 
# for hash (individual mouse), litter, and hash pool.

.int_model_ordered <- function(data, k, basis) {
  bam(value ~ oGeno +
    s(age, bs = basis, k = k) +
    s(age, by = oGeno, bs = basis, k = k) +
    s(age, sex, bs = "re") +
    s(litter, bs = "re") +
    s(hash, bs = "re") +
    s(hpool, bs="re"),
  data = data, method = "REML", family = nb
  )
}

.int_model_mean <- function(data, k, basis) {
  bam(value ~ geno +
        s(age, m = 2, bs = basis, k=k) +
        s(age, by = oGeno, m = 1, bs = basis, k=k) +
        s(age, sex, bs = "re") +
        s(litter, bs = "re") +
        s(hash, bs = "re") +
        s(hpool, bs="re"),
      data = data, method = "REML", family = nb
  )
}


.get_long_data <- function(seur_obj, min.pct = 0.01) {
  
  # filter genes to those which are detected in at least 1% of cells more than once
  min_cells <- round(min.pct*ncol(seur_obj@assays$SCT@counts))
  keep_genes <- Matrix::rowSums(seur_obj@assays$SCT@counts > 1) > min_cells
  
  
  # filter genes and add metadata to gene expression matrix
  data <-
    seur_obj@assays$SCT@counts[keep_genes,] %>%
    Matrix::t() %>%
    cbind(seur_obj[[]] %>% 
            dplyr::select(
              age, predicted.id, hash_id, sex_call, hash_pool, diet, geno
            ) %>% 
            mutate(
              age = as.numeric(age), 
              sex = factor(sex_call)
            )
    )
  
  # convert matrix to a nested data.table, each gene is a row and all 
  # info is nested in the data column. data.table is significantly faster and
  # we run out of memory using dplyr::nest
  dt.long <- 
    data %>% 
    filter(
      age < 40,
      age > 0 
    ) %>% 
    mutate(
      geno = factor(geno),
      oGeno = factor(geno, levels = c("wt", "mut"), ordered=T),
      group = factor(case_when(
        diet == "HFD" ~ "mHFD",
        geno == "mut" ~ "ob/ob",
        T ~ "gWT")),
      oGroup = factor(group, ordered = T),
      hash = droplevels(hash_id),
      clus = droplevels(factor(predicted.id)),
      litter = droplevels(factor(str_split(hash, "-", simplify = T)[, 1])),
      hpool = droplevels(factor(hash_pool))
    ) %>% 
    dplyr::select(
      -predicted.id, -hash_id, -sex_call, -hash_pool, -diet 
    ) %>% 
    data.table::data.table() %>% 
    tidyfast::dt_pivot_longer(., c(-age, -clus, -hash, -sex, 
                                   -hpool, -litter, -group, -oGroup,
                                   -geno, -oGeno)) %>% 
    tidyfast::dt_nest(., name)
}

.gen_preds <- function(data, model) {
  new_data <- tidyr::expand(data, nesting(hash, age, litter, hpool, oGeno, sex))
  preds <- bind_cols(new_data, as.data.frame(predict(model, newdata = new_data, se.fit = TRUE,
                                                     type="response",  
                                                     exclude = c("s(hash)", "s(litter)",
                                                                 "s(hpool)", "s(age):sexM",
                                                                 "s(age):sexF"))))
  return(preds)
}

# generate model outputs
.get_intmodel_summs <- function(x,k,basis) {
  int_gam <- .int_model_ordered(x, k, basis)
  int_par <- .int_model_mean(x, k, basis)
  summ_gam <- summary(int_gam)
  summ_par <- summary(int_par)
  preds <- .gen_preds(data = x, model = int_gam)
  diff_gam <- tidymv::get_smooths_difference(int_gam, age, list(oGeno = c("mut", "wt")), exclude_random = F)
  return(tibble(summ_gam = list(summ_gam), summ_par = list(summ_par), diff_gam = list(diff_gam), preds = preds))
}


fit_interaction <- function(dat, workers, k, basis) {
  
  # filter genes to those which are detected in at least 1% of cells more than once
  # filter genes and add metadata to gene expression matrix
  dt.long <- .get_long_data(dat, min.pct = 0.01)
  
  print("requesting cores")
  
  if(!is.null(workers)) {
    plan("multisession", workers = workers)
  } else {
    plan("sequential")
  }
  
  print("got cores")
  
  int_age <- furrr::future_map(dt.long$data, ~.get_intmodel_summs(.x, k, basis), .progress = T)

  return(int_age)
}
