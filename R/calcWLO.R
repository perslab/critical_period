##' .. content for \description{
##' goal of this function is to find genes that are specifically expressed in each cluster, it is also
##' significantly faster than the FindMarkers function from Seurat
##' } (no empty lines) ..
##'
##' .. content for \details{
##' pseudobulk calculation inspired by: https://jef.works/blog/2020/04/06/quickly-creating-pseudobulks/
##' weighted log odds from: https://github.com/juliasilge/tidylo
##' }
##'
##' @title
##' @param srt seurat object
##' @param group define groups which you want to split; default is srt ident
##' @param compare vector of groups you want to compare #add in
##' @return
##' @author dylanmr
##' @export

require(tidylo)

.countsperclus <- function(seur_obj, group=group, min.cell=100) {
  if(is.null(group)) {
    vec <- factor(as.character(Idents(seur_obj)))
  } else {
    vec <- factor(seur_obj@meta.data[[group]])
  }
  mat.sparse <- seur_obj@assays$SCT@counts
  mm <- model.matrix(~ 0 + vec)
  colnames(mm) <- paste0("clus_", levels(vec))
  mm <- mm[,colSums(mm)>min.cell]
  mat.sum <- mat.sparse %*% mm
  keep <-  Matrix::rowSums(mat.sum > 0) >= ncol(mat.sum)/3
  mat.sum <- mat.sum[keep, ]
  return(mat.sum)
}


calcWLO <- function(seur_obj, prior, group= NULL, id1 = NULL, id2 = NULL, ...) {

  cpg <- .countsperclus(seur_obj, group = group, ...)
  
  if(grepl("un", prior)) {
    cpg <-
      cpg %>% 
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene) %>%
      dplyr::rename(group = name) %>%
      mutate(group = as.factor(group))
    
    if(is.null(id1)) {
      cpg %>% 
        bind_log_odds(set = group, feature = gene, n = value, 
                      uninformative = TRUE, unweighted = TRUE) -> dat
    } else if(!is.null(id1)) {
      if(is.null(id2)) {
        id1 <- paste0("clus_", id1)
        cpg %>% 
          mutate(group = if_else(group == id1, "group_1", "group_2")) %>% 
          group_by(gene, group) %>% 
          summarise(value = sum(value)) %>% 
          bind_log_odds(set = group, feature = gene, n = value, 
                        uninformative = TRUE, unweighted = TRUE) -> dat
      } else {
        id1 <- paste0("clus_", id1)
        id2 <- paste0("clus_", id2)
        cpg %>% 
          mutate(group = case_when(group %in% id1 ~ "group_1",
                                   group %in% id2 ~ "group_2",
                                   T ~ "remove")) %>%
          group_by(gene, group) %>% 
          summarise(value = sum(value)) %>% 
          filter(group!="remove") %>% 
          bind_log_odds(set = group, feature = gene, n = value, 
                        uninformative = TRUE, unweighted = TRUE) -> dat
      }
    }
  } else {
    cpg %>% 
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene) %>%
      dplyr::rename(group = name) %>%
      mutate(group = as.factor(group)) %>%
      bind_log_odds(set = group, feature = gene, n = value) -> dat
  }
  
  return(dat)
}
