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



fit_age_mod <- function(seur_obj, group = "wt", workers=40) {
  
  if (group == "chow") {
    dat <- subset(seur_obj, subset = diet == "CHOW" & geno == "wt")
  } else if (group == "mut") {
    dat <- subset(seur_obj, subset = geno == "mut")
  } else if (group == "hfd") {
    dat <- subset(seur_obj, subset = diet == "mHFD")
  } else if (group == "wt") {
    dat <- subset(seur_obj, subset = geno == "wt") 
  }

  data <-
    dat@assays$SCT@counts[c("Agrp", "Npy", "Lncpint"),] %>%
    Matrix::t() %>%
    cbind(dat[[]] %>% dplyr::select(age, predicted.id, hash_id, sex_call, hash_pool)) %>%
    pivot_longer(c(-age, -predicted.id, -hash_id, -sex_call, -hash_pool)) %>%
    filter(age < 40)

  not_arc <- paste(paste0("n", c("01", "06", 16:18, 29:31, 33:34)), collapse = "|")

  keep <-
    dat[[]] %>%
    count(predicted.id) %>%
    filter(n > 100) %>%
    pull(predicted.id)

  keep <- keep[!grepl(not_arc, keep)]

  nest_data <-
    data %>% 
    filter(predicted.id %in% keep) %>%
    mutate(
      age = as.numeric(age),
      hash = droplevels(hash_id),
      sex = factor(sex_call),
      clus = droplevels(factor(predicted.id)),
      litter = droplevels(factor(str_split(hash, "-", simplify = T)[, 1])),
      hpool = droplevels(factor(hash_pool))
    ) %>%
    nest(data = !name)
  
  if(!is.null(workers)) {
    plan("multisession", workers = workers)
  } else {
    plan("sequential")
  }
  
  shared_age <-
    furrr::future_map(nest_data$data, function(x) {
      age_shared <- .fit_age_re(x)
      age_spec <- .fit_cluster_slope_re(x)
      summ_shared <- summary(age_shared)
      summ_spec <- summary(age_spec)
      return(tibble(summ_shared = list(summ_shared), summ_spec = list(summ_spec)))
    }, .progress = T) %>%
    bind_rows()
  
  shared_age <-  
    shared_age %>% 
    mutate(
      p_age = map_dbl(summ_shared, ~ .x$p.pv[["age"]]),
      p_slope = map(summ_spec, ~ .x$p.pv[2:length(.x$p.pv)]),
      padj_age = p.adjust(p_age),
      age_coef = map_dbl(summ_shared, ~ .x$p.coeff[["age"]]),
      slope_coef = map(summ_spec,  ~ .x$p.coeff[2:length(.x$p.coeff)]),
      gene = nest_data$name
    )
  
  return(shared_age)
}
