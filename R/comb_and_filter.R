##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param sce_objects
##' @return
##' @author dylanmr
##' @export


comb_and_filter <- function(sce_objects, meta, cellnames) {
  
  # remove alternative experiment data for each object
  sce <- 
    map(sce_objects, removeAltExps) %>% 
    map2(cellnames, function(x,y) {
      y <- str_match(y, "/data/sc-10x/data-kallisto/SCOP_37/kallisto/(.*)/.*")[2]
      colnames(x) <- paste0(colnames(x), ".", y)
      return(x)
    })
  
  # extract an overlapping set of rownames for each object
  overlap <- 
    map(sce, rownames) %>% 
    Reduce(intersect, .)
  
  # filter all objects to have same gene names and bind to form final object
  sce_comb <- 
    map(sce, function(x) x[overlap,]) %>% 
    do.call(cbind, .) 
  
  # save ensembl names in rowdata 
  rowData(sce_comb)$ensembl <- rownames(sce_comb)
  
  # replace rownames with gene names
  rownames(sce_comb) <- rowData(sce_comb)$gene_name
  sce_comb <- sce_comb[!duplicated(rownames(sce_comb)),]
  
  # remove all cells labelled as discard by QC
  sce_comb <- sce_comb[,sce_comb$discard==F]
  
  # retain only cells classified as singlets
  sce_comb <- sce_comb[,sce_comb$hto_global=="Singlet"]
  
  # remove manually called doublets
  sce_comb <- sce_comb[,sce_comb$predicted_dub_std == F]
  print(table(sce_comb$predicted_dub_cut))
  sce_comb <- sce_comb[,sce_comb$predicted_dub_cut == F]

  # some cell names are duplicated, make unique
  colnames(sce_comb) <- make.unique(colnames(sce_comb), sep="_")

  # add metadata to object
  colData(sce_comb) <- cbind(colData(sce_comb), meta[match(colData(sce_comb)$hash_id, meta$hash),])
  colData(sce_comb)$hash_id <- droplevels(colData(sce_comb)$hash_id)
  
  # create column to identify if cell is from ob or wt
  reducedDim(sce_comb, "PCA") <- NULL
  return(sce_comb)

}
