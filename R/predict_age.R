##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param seur_obj seurat object
##' @param features feature on which we are predicting age
##' @return
##' @author dylanmr
##' @export
predict_age <- function(seur_obj, features = genes) {
  
  df <- data.frame(age = seur_obj$age, Matrix::t(seur_obj@assays$SCT@counts[features,]))
  df_sparse <- Matrix::Matrix(data.matrix(df), sparse=T)
  mod <- ranger::ranger(dependent.variable.name="age", data = df_sparse)
  page <- mod$predictions %>% 
    data.frame(page = .) %>% 
    mutate(cells = colnames(seur_obj))
  return(list(page, mod))

}
