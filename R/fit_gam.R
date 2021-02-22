##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data nested tibble with gene names in column1 and long-form
##' data in column2 with age and group
##' @param p to pass progress bar 
##' @param p_grid parameter grid to predict values on
##' @return
##' @author dylanmr
##' @export
fit_gam <- function(data, p, p_grid) {

  mod <- mgcv::gam(value ~ s(age, k=5, m=2, bs="tp") +
                     s(age, by=group, k=5, m=1, bs="tp") +
                     s(group, bs="re"),
                   data=data, method="REML", family="nb")
  # setup prediction data
  pred <- cbind(p_grid, predict(mod, p_grid, se.fit=TRUE, type="response"))
  p()
  return(list(pred, mod))

}


