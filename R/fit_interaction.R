##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param seur_obj seurat object
##' @param k number of basis functions to fit with GAM
##' @param family distribution on which to model data
##' @param ordered should parametric term in the model be ordered
##' @param basis type of function to model smooth terms
##' @param nfeat number of variable features to find
##' 
##' @return
##' @author dylanmr
##' @export

library(mgcv)
library(tidyverse)
library(Seurat)
library(furrr)
library(future)
library(Matrix)

# model expression of all detectable genes in agrp neurons over age.
# Use mgcv to speed up calculations in large data-sets. We are using
# an ordered factor to calculate "difference smooths". As seen here:
# https://fromthebottomoftheheap.net/2017/12/14/difference-splines-ii/
# include a random intercept for cluster litter, hash, and hashed pool.
# Sex is included as a by factor - need to check if this is the right way to include this
# Include a random intercept
# for hash (individual mouse), litter, and hash pool.

.interaction_model <- function(data, k, basis, family, ordered) {
  if (ordered == T) {
    bam(value ~ oGeno +
      s(age, bs = basis, k = k, m=2) +
      s(age, by = oGeno, bs = basis, k = k, m=1) +
      s(age, sex, bs = "re") +
      s(litter, bs = "re") +
      s(hash, bs = "re") +
      s(hpool, bs = "re"),
    data = data, family = family, method="fREML", discrete=T
    )
  } else {
    bam(value ~ geno +
      s(age, m = 2, bs = basis, k = k) +
      s(age, by = oGeno, m = 1, bs = basis, k = k) +
      s(age, sex, bs = "re") +
      s(litter, bs = "re") +
      s(hash, bs = "re") +
      s(hpool, bs = "re"),
    data = data, method = "fREML", family = family, discrete=T
    )
  }
}

# catch errors
.safe_model <- safely(.interaction_model)

# get variable genes for subsetted object and
# remove non-coding genes
.get_var_genes <- function(seur_obj, nfeat) {
  keep_genes <- !grepl("Rik", rownames(seur_obj))
  seur_obj <- subset(seur_obj, features = rownames(seur_obj)[keep_genes])
  x <- FindVariableFeatures(seur_obj, nfeatures = nfeat)
  x <- x@assays$SCT@var.features
  return(x)
}

# format data for modeling
.get_long_data <- function(seur_obj, nfeat) {

  # filter genes to those which are detected in at least 1% of cells more than once
  var_feats <- .get_var_genes(seur_obj, nfeat)

  # filter genes and add metadata to gene expression matrix
  # get corrected counts from SCT counts matrix to use for modelling
  data <-
    seur_obj@assays$SCT@counts[var_feats, ] %>%
    Matrix::t() %>%
    cbind(seur_obj[[]] %>%
      dplyr::select(
        age, predicted.id, hash_id, sex_call, hash_pool, diet, geno
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
      age < 40,
      age > 0
    ) %>%
    mutate(
      geno = factor(geno),
      oGeno = factor(geno, levels = c("wt", "mut"), ordered = T),
      group = factor(case_when(
        diet == "HFD" ~ "mHFD",
        geno == "mut" ~ "ob/ob",
        T ~ "gWT"
      )),
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
    tidyfast::dt_pivot_longer(., c(
      -age, -clus, -hash, -sex,
      -hpool, -litter, -group, -oGroup,
      -geno, -oGeno
    )) %>%
    tidyfast::dt_nest(., name)
}

# generate predictions from GAM models

.gen_preds <- function(data, model) {
  new_data <- tidyr::expand(data, nesting(hash, age, litter, hpool, oGeno, sex))
  preds <- bind_cols(new_data, as.data.frame(predict(model,
    newdata = new_data, se.fit = TRUE,
    type = "response",
    exclude = c(
      "s(hash)", "s(litter)", "s(hpool)", "s(age,sex)"
    )
  )))
  return(preds)
}

# run models
.run_model <- function(x, k, basis, family) {
  int_gam <- .safe_model(x, k, basis, family, ordered = T) 
  int_par <- .safe_model(x, k, basis, family, ordered = F) 
  return(list(int_gam, int_par))
}

# generate model summaries
.gen_summary <- function(int_ord, int_mean, metadata) {
  summ_gam <- summary(int_ord)
  summ_par <- summary(int_mean)
  preds <- .gen_preds(data = metadata, model = models[[1]])
  diff_gam <- tidymv::get_smooths_difference(models[[1]], age, list(oGeno = c("mut", "wt")), exclude_random = F)
  return(tibble(summ_gam = list(summ_gam), summ_par = list(summ_par), diff_gam = list(diff_gam), preds = preds))
}


fit_interaction <- function(seur_obj, workers, k, basis, family, nfeat = 3000) {

  # filter genes to those which are detected in at least 1% of cells more than once
  # filter genes and add metadata to gene expression matrix
  dt.long <- .get_long_data(seur_obj, nfeat)

  print("requesting cores")

  if (!is.null(workers)) {
    plan("multisession", workers = workers)
  } else {
    plan("sequential")
  }

  print("got cores")

  interaction_ordered <- furrr::future_map(dt.long$data, ~ .safe_model(.x, k, basis, family, ordered = T), .progress = T) %>% 
    setNames(dt.long$name) %>% 
    keep(~is.null(.x$error)) %>% 
    map(1)
  interaction_mean <- furrr::future_map(dt.long$data, ~ .safe_model(.x, k, basis, family, ordered = F), .progress = T) %>% 
    setNames(dt.long$name) %>% 
    keep(~is.null(.x$error)) %>% 
    map(1)
  
  return(interaction_models)
}
