##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param seur
##' @param knn
##' @param col
##' @return
##' @author dylanmr
##' @export
##' 
##' calculate p-value with pbinom
##' inspired by https://github.com/kstreet13/bioc2020trajectories/blob/master/R/imbalance_score.R
##' probably adjust prob for each timepoint?

enrich_score <- function(seur, knn, col, dims) {

  vec <- factor(seur[[col]][,1])
  prop <- as.vector(table(vec) / length(vec))[1]
  group <- levels(factor(vec))[1]
  g <- knn$index
  neighb_mat <- matrix(vec[g], ncol=ncol(g))
  neighb_mat <- neighb_mat == group
  scores <- unlist(lapply(rowSums(neighb_mat), function(x) qnorm(pbinom(x, ncol(g)+1, prop, lower.tail = F))))
  scores <-
    if_else(is.infinite(scores), 
            true = max(scores[is.finite(scores)]) + .5,
            false = scores)
  p_scores <- scales::rescale(scores, c(-1,1))
  # build a data frame to predict labels
  dat <- data.frame(vec = vec, Embeddings(seur, reduction = "pca")[,c(1:dims)])
  # predict labels w/random forest and retain probabilities
  rf_pred <- predictions(ranger::ranger(vec ~ .,data = dat, probability = T))[,1]
  score_name <- paste0(col, "_enrichscore")
  cont_name <- paste0(col, "_rf")
  return(list(binom_scores = p_scores, rf_scores = rf_pred))

}
