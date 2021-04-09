##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param seur_obj seurat object
##' @param group which set of data would you like to model (chow - both wt and mHFD, wt - just chow, mHFD, mut - ob)
##' @param workers number of cores to use for parallelization of model fitting
##' @return
##' @author dylanmr
##' @export

library(mgcv)
library(tidyverse)
library(Seurat)
library(furrr)
library(future)
library(Matrix)
library(progressr)
handlers(global = TRUE)

# build a model to get slope of each gene within each cluster
# by age. Use mgcv to speed up calculations in large data-sets
# include a random intercept for cluster (we dont care that expression 
# may start at different levels across groups). Include a random intercept 
# for hash (individual mouse), litter, and hash pool. Will have to add 
# a random term for sequencing pool. fixed effect to control for the effect 
# of sex

.fit_cluster_slope_re <- function(data) {
  bam(value ~ 0 + age:clus + sex +
    s(clus, bs="re") +
    s(hash, bs = "re") +
    s(litter, bs="re") +
    s(hpool, bs="re"),
  data = data, method = "fREML", family="nb", discrete=T)
}


# build a model to get slope of each gene across all populations by
# by age. Use mgcv to speed up calculations in large data-sets
# include a random intercept for cluster and a random slope that is not correlated 
# to the level of the intercept. Include a random intercept 
# for hash (individual mouse), litter, and hash pool. Will have to add 
# a random term for sequencing pool. fixed effect to control for the effect 
# of sex

.fit_age_re <- function(data) {
  bam(value ~ age + sex +
    s(clus, bs = "re") +
    s(age, clus, bs = "re") +
    s(hash, bs = "re") +
    s(litter, bs="re") +
    s(hpool, bs="re"),
  data = data, method = "fREML", family="nb", discrete=T)
}

# generate model outputs

.gen_preds <- function(data, model) {
  new_data <- tidyr::expand(data, nesting(hash, age, litter, hpool, sex),clus = unique(clus))
  preds <- bind_cols(new_data, as.data.frame(predict(model, newdata = new_data, se.fit = TRUE)))
  return(preds)
}

.get_model_summs <- function(x) {
  p <- progressr::progressor(along = x)
  future.apply::future_lapply(x, function(df, ...) {
    age_shared <- .fit_age_re(df)
    age_spec <- .fit_cluster_slope_re(df)
    summ_shared <- summary(age_shared)
    summ_spec <- summary(age_spec)
    preds <- .gen_preds(data = df, model = age_shared)
    p()
    return(tibble(summ_shared = list(summ_shared), summ_spec = list(summ_spec), preds = list(preds)))
  })
}

fit_maturation_model <- function(seur_obj, group = "wt", workers=40, min.pct = .01) {
  
  if (group == "chow") {
    dat <- subset(seur_obj, subset = diet == "CHOW" & geno == "wt")
  } else if (group == "mut") {
    dat <- subset(seur_obj, subset = geno == "mut")
  } else if (group == "hfd") {
    dat <- subset(seur_obj, subset = diet == "mHFD")
  } else if (group == "wt") {
    dat <- subset(seur_obj, subset = geno == "wt") 
  }
  
  # filter genes to those which are detected in at least 1% of cells more than once
  min_cells <- round(min.pct*ncol(dat@assays$SCT@counts))
  keep_genes <- Matrix::rowSums(dat@assays$SCT@counts > 1) > min_cells
  
  # only model on clusters with > 100 cells in the arcuate
  not_arc <- paste(paste0("n", c("01", "06", 16:18, 29:31, 33:34)), collapse = "|")
  keep_clus <-
    dat[[]] %>%
    dplyr::count(predicted.id) %>%
    filter(n > 100) %>%
    pull(predicted.id)
  keep_clus <- keep_clus[!grepl(not_arc, keep_clus)]
  
  # filter genes and add metadata to gene expression matrix
  data <-
    dat@assays$SCT@counts[keep_genes,] %>%
    Matrix::t() %>%
    cbind(dat[[]] %>% 
            dplyr::select(
              age, predicted.id, hash_id, sex_call, hash_pool
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
      predicted.id %in% keep_clus,
      age < 40,
      age > 0 
      ) %>% 
    mutate(
      hash = droplevels(hash_id),
      clus = droplevels(factor(predicted.id)),
      litter = droplevels(factor(str_split(hash, "-", simplify = T)[, 1])),
      hpool = droplevels(factor(hash_pool))
      ) %>% 
    dplyr::select(
      -predicted.id, -hash_id, -sex_call, -hash_pool 
    ) %>% 
    data.table::data.table() %>% 
    tidyfast::dt_pivot_longer(., c(-age, -clus, -hash, -sex, -hpool, -litter)) %>% 
    tidyfast::dt_nest(., name)
  
  print("requesting cores")
  
  if(!is.null(workers)) {
    plan("multisession", workers = workers)
  } else {
    plan("sequential")
  }
  
  print("got cores")
  
  shared_age <- bind_rows(.get_model_summs(dt.long$data))
  
  shared_age <-  
    shared_age %>% 
    mutate(
      p_age = map_dbl(summ_shared, ~ .x$p.pv[["age"]]),
      p_slope = map(summ_spec, ~ .x$p.pv[2:length(.x$p.pv)]),
      padj_age = p.adjust(p_age),
      age_coef = map_dbl(summ_shared, ~ .x$p.coeff[["age"]]),
      slope_coef = map(summ_spec,  ~ .x$p.coeff[2:length(.x$p.coeff)]),
      gene = dt.long$name
    )
  
  return(shared_age)
}
