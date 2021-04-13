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

.fit_age <- function(data, family, byclus) {
  if (byclus == T) {

    # think about markov random fields
    # build a model to get slope of each gene within each cluster
    # by age. Use mgcv to speed up calculations in large data-sets
    # include a random intercept for cluster (we dont care that expression
    # may start at different levels across groups). Include a random intercept
    # for hash (individual mouse), litter, and hash pool. Will have to add
    # a random term for sequencing pool. fixed effect to control for the effect
    # of sex

    bam(value ~ 0 + age:clus +
      s(clus, sex, bs = "re") +
      s(clus, bs = "re") +
      s(hash, bs = "re") +
      s(litter, bs = "re") +
      s(hpool, bs = "re"),
    data = data, method = "fREML", family = family, discrete = T
    )
  } else {

    # build a model to get slope of each gene across all populations by
    # by age. Use mgcv to speed up calculations in large data-sets
    # include a random intercept for cluster and a random slope that is not correlated
    # to the level of the intercept. Include a random intercept
    # for hash (individual mouse), litter, and hash pool. Will have to add
    # a random term for sequencing pool. fixed effect to control for the effect
    # of sex

    bam(value ~ age +
      s(clus, bs = "re") +
      s(age, clus, bs = "re") +
      s(age, sex, bs = "re") +
      s(hash, bs = "re") +
      s(litter, bs = "re") +
      s(hpool, bs = "re"),
    data = data, method = "fREML", family = family, discrete = T
    )
  }
}

# catch errors
.safe_model <- safely(.fit_age)

# format data for modeling
.get_long_data <- function(seur_obj, nfeat) {

  # filter genes to those which are detected in at least 1% of cells more than once
  var_feats <- .get_var_genes(seur_obj, nfeat)

  # only model on clusters with > 100 cells in the arcuate
  not_arc <- paste(paste0("n", c("01", "06", 16:18, 29:31, 33:34)), collapse = "|")
  keep_clus <-
    dat[[]] %>%
    dplyr::count(predicted.id) %>%
    filter(n > 100) %>%
    pull(predicted.id)
  keep_clus <- keep_clus[!grepl(not_arc, keep_clus)]

  # filter genes and add metadata to gene expression matrix
  # get corrected counts from SCT counts matrix to use for modelling

  # filter genes and add metadata to gene expression matrix
  data <-
    dat@assays$SCT@counts[keep_genes, ] %>%
    Matrix::t() %>%
    cbind(dat[[]] %>%
      dplyr::select(
        age, predicted.id, hash_id, sex_call, hash_pool
      ) %>%
      mutate(
        age = as.numeric(age),
        sex = factor(sex_call)
      ))

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
}


# run models
.run_model <- function(x, family) {
  age_shared <- .fit_age(x, family, byclus == F)
  age_spec <- .fit_age(x, family, byclus == T)
}

# generate model summaries
.gen_summary <- function(models) {
  summ_shared <- summary(models[[1]])
  summ_spec <- summary(models[[2]])
  return(tibble(summ_shared = list(summ_shared), summ_spec = list(summ_spec)))
}


fit_maturation_model <- function(seur_obj, group = "wt", workers = 40, nfeat = 3000, family) {
  if (group == "chow") {
    dat <- subset(seur_obj, subset = diet == "CHOW" & geno == "wt")
  } else if (group == "mut") {
    dat <- subset(seur_obj, subset = geno == "mut")
  } else if (group == "hfd") {
    dat <- subset(seur_obj, subset = diet == "mHFD")
  } else if (group == "wt") {
    dat <- subset(seur_obj, subset = geno == "wt")
  }

  dt.long <- .get_long_data(dat, nfeat)

  print("requesting cores")

  if (!is.null(workers)) {
    plan("multisession", workers = workers)
  } else {
    plan("sequential")
  }

  print("got cores")

  age_models <- furrr::future_map(dt.long$data, ~ .run_model(.x, family), .progress = T)  %>%  
    setNames(dt.long$name) %>% 
    keep(~is.null(.x$error)) %>% 
    map(1)
  
  model_summs <- furrr::future_map(age_models, .gen_summary) %>% 
    bind_rows(.id="gene")
  
  return(shared_age)
}
