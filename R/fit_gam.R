##' \description{ This function is used to interpolate gene expression over time
##' With time series data with no branching, this may make more sense
##' To use than to try to construct cell pseudotime, gam model is 
##' inspired from this link: https://peerj.com/articles/6876/?td=tw#supplemental-information
##' here I use a global smoother to share information across groups, but allow for an
##' individual smooth and individual wiggliness across groups} 
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


