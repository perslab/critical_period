##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param rawdat
##' @param preds
##' @return
##' @author dylanmr
##' @export
plot_gam_preds <- function(rawdat, preds) {

  p <- 
    ggplot(data=rawdat, aes(x=age, y=value, group=group)) +
    geom_ribbon(aes(ymin = fit - 2*se.fit, ymax= fit + 2*se.fit, x=age,  fill = group, group = group),
                data=preds, 
                alpha=0.75, 
                inherit.aes=FALSE) +
    geom_line(aes(y=fit), data=preds, color="black", linetype="dashed") +  
    geom_jitter(data = rawdat %>% group_by(group, hash, age) %>% 
                  summarise(mean = mean(value)) %>% filter(age < 40, age >0),
                aes(x=age, y=mean, fill=group), shape=21) +
    theme_minimal()
  
  return(p)
}
