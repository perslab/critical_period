##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param seur_obj
##' @param dims
##' @param k
##' @return
##' @author dylanmr
##' @export
detectHotSpot <- function(seur_obj, dims, k, col) {
  
  library(tidygraph)
  print("Building KNN-graph")
  g <- build_graph(seur_obj, dims = dims, k =k, reduction = "pca")
  print("Scoring Cells")
  preds <- enrich_score(seur_obj, knn = g, col = col, dims=dims)
  preds <- map(preds, ~scales::rescale(.x, c(-1,1)))
  print("Converting to igraph object")
  ig <- bluster::makeKNNGraph(x = Embeddings(seur_obj, reduction = "pca")[,c(1:dims)], k=k)
  print("Calculating edge weights")
  E(ig)$weight <- 
    tidygraph::as_tbl_graph(ig) %>% 
    activate(edges) %>% 
    mutate(score_rf = (preds$rf_scores[from]^2 + preds$rf_scores[to]^2)/4,
           score_freq = (preds$binom_scores[from] + (preds$binom_scores[to]))^2) %>% 
    pull(score_rf) %>% 
    as.numeric()
  V(ig)$strength <- preds$rf_scores
  print("Finding communities")
  ig_res <- igraph::cluster_fast_greedy(ig)
  
  df <- data.frame(rf_score = preds$rf_scores, binom_score = preds$binom_scores, hs_clus = ig_res$membership) 
  rownames(df) <- colnames(seur_obj)
  seur_obj <- AddMetaData(seur_obj, df) 
  return(seur_obj)
  
  }
