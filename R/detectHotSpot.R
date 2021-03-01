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
##' 
##' 


library(tidygraph)
library(ranger)

.build_graph <- function(seur_obj, dims, k, reduction = "pca") {
  
  df <- Embeddings(seur_obj, reduction = reduction)[,c(1:dims)]
  g <- BiocNeighbors::findKNN(df, k = k, get.distance=T, get.index=T)
  return(g)
  
}

.enrich_score <- function(seur_obj, knn, col, dims, features = NULL) {
  
  vec <- factor(seur_obj[[col]][,1])
  prop <- as.vector(table(vec) / length(vec))[1]
  group <- levels(factor(vec))[1]
  g <- knn$index
  neighb_mat <- matrix(vec[g], ncol=ncol(g))
  neighb_mat <- neighb_mat == group
  scores <- unlist(lapply(rowSums(neighb_mat), function(x) qnorm(pbinom(x, ncol(g)+1, prop, lower.tail = F))))
  binom_scores <- if_else(is.infinite(scores), 
              true = max(scores[is.finite(scores)]) + .5,
              false = scores)
  # build a data frame to predict labels
  if(is.null(features)) {
    dat <- data.frame(vec = vec, Matrix::t(seur_obj@assays$SCT@scale.data))
  } else if(features=="pca") {
    dat <- data.frame(vec = vec, Embeddings(seur_obj, reduction = "pca")[,c(1:dims)])
  }
  # predict labels w/random forest and retain probabilities
  rf_pred <- predictions(ranger::ranger(vec ~ .,data = dat, probability = T))[,1]
  score_name <- paste0(col, "_enrichscore")
  cont_name <- paste0(col, "_rf")
  return(list(binom_scores = binom_scores, rf_scores = rf_pred))
  
}
  
.calc_edgeweight <- function(preds, g, w = 2, score = "rf") {
  weights <- 
    tidygraph::as_tbl_graph(g) %>% 
    activate(edges) %>% 
    mutate(score_rf = (abs(preds$rf_scores[from]) + abs(preds$rf_scores[to]))^w,
           score_bin = (preds$binom_scores[from]*w + preds$binom_scores[to]*w)^w) %>% 
    as_tibble()
  
  if(score == "rf"){
    return(weights$score_rf)
  }else{
    return(weights$score_bin)
  }
}  
  
.perm_edges <- function(g, nperm) {
  resample_rf <- replicate(nperm, data.frame(rf_scores = sample(V(g)$rf_score), binom_scores = sample(V(g)$binom_score)), simplify = F)
  perm_edges <- lapply(resample_rf, .calc_edgeweight, g=g) %>% 
    bind_cols() %>% 
    mutate(from = tidygraph::as_tbl_graph(g) %>% 
             activate(edges) %>% 
             pull(from),
           clus = V(g)$group[from]) %>% 
    dplyr::select(-from) %>% 
    pivot_longer(-clus) %>% 
    group_by(clus, name) %>% 
    summarise(sum = sum(value)) 
  return(perm_edges)
}

.calc_significance <- function(g, perm_edges) {
  true_val <-
    as_tbl_graph(g) %>% 
    activate(edges) %>% 
    as_tibble() %>% 
    mutate(clus = V(g)$group[from]) %>% 
    group_by(clus) %>% 
    summarise(true_mod = sum(weight))
  dat <- 
    perm_edges %>% 
    nest(-clus) %>% 
    mutate(ecdf = map(data, ~ecdf(.x$sum)),
           sumstats = map(data, function(x) data.frame(mean = mean(x$sum, na.rm=T), sd = sd(x$sum, na.rm=T)))) %>% 
    left_join(true_val) %>% 
    mutate(auc = unlist(map2(ecdf, true_mod, ~.x(.y))),
           pval = unlist(map2(sumstats, true_mod, ~pnorm(.y, .x$mean, .x$sd, lower.tail = F))),
           padj = p.adjust(pval)) %>% 
    dplyr::select(-data)
}

detectHotSpot <- function(seur_obj, dims, k, col, return_extra, nperm, ...) {
  
  print("Building KNN-graph")
  g <- .build_graph(seur_obj, dims = dims, k = k, reduction = "pca")
  print("Scoring Cells")
  preds <- .enrich_score(seur_obj, knn = g, col = col, dims=dims)
  preds <- purrr::map(preds, ~scales::rescale(.x, c(-1,1)))
  print("Converting to igraph object")
  ig <- bluster::makeKNNGraph(x = Embeddings(seur_obj, reduction = "pca")[,c(1:dims)], k=k)
  print("Calculating edge weights")
  E(ig)$weight <- .calc_edgeweight(preds, ig, ...)
  print("Finding communities")
  ig_res <- igraph::cluster_fast_greedy(ig)
  # assign values to nodes for downstream analysis
  V(ig)$group <- ig_res$membership
  V(ig)$rf_score <- preds$rf_score
  V(ig)$binom_score <- preds$binom_score
  # permute edges for null distribution
  print("Permuting edges and calculating significance")
  perm_edges <- .perm_edges(ig, nperm)
  sig_cluster <- .calc_significance(ig, perm_edges)
  # arrange output data
  df <- data.frame(rf_score = preds$rf_scores, binom_score = preds$binom_scores, hs_clus = ig_res$membership) %>% 
    left_join(sig_cluster %>% dplyr::select(clus, padj), by = c("hs_clus" = "clus")) %>% 
    mutate(padj = if_else(padj > 0.05, NA_character_, paste0(hs_clus, "_sig")))
  rownames(df) <- colnames(seur_obj)
  seur_obj <- AddMetaData(seur_obj, df) 
  # decide what to return
  if(return_extra) {
    return(list(g = ig, perm = sig_cluster, seur = seur_obj)) 
  } else {
    return(seur_obj)
  }
  
}

