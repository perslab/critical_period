##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param sce_objects
##' @return
##' @author dylanmr
##' @export
gen_full_object <- function(sce_objects) {
  # remove alternative experiment data for each object
  sce <- sce_objects %>% 
    purrr::map(removeAltExps)
  # extract an overlapping set of rownames for each object
  overlap <- purrr::map(sce, rownames) %>% 
    Reduce(intersect, .)
  # filter all objects to have same gene names and bind to form final object
  sce_comb <- purrr::map(sce, function(x) x[overlap,]) %>% 
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
  # remove UNK hash
  sce_comb <- sce_comb[,sce_comb$hash_id!="UNK"]
  sce_comb <- sce_comb[,sce_comb$hash_id!="ObA-1"]
  # remove manually called doublets
  sce_comb <- sce_comb[,sce_comb$predicted_dub_std == F]
  print(table(sce_comb$predicted_dub_cut))
  sce_comb <- sce_comb[,sce_comb$predicted_dub_cut == F]
  # some cell names are duplicated, make unique
  colnames(sce_comb) <- make.unique(colnames(sce_comb), sep="_")
  # pull in metadata
  dat <- sync_metadata()
  # add metadata to object
  colData(sce_comb) <- cbind(colData(sce_comb), dat[match(colData(sce_comb)$hash_id, dat$hash_id),])
  # create column to identify if cell is from ob or wt
  colData(sce_comb)$geno <- ifelse(grepl("ob", sce_comb$hash_id, ignore.case = T), yes = "mut", no = "wt")
  # run SCTransform - this function creates a matrix of corrected UMI, and log1p of corrected UMI
  # this takes top variable genes and extracts pearson residuals, stored as the used assay for 
  # downstream analysis
  seur_comb <- as.Seurat(sce_comb)
  seur_comb %>% 
    SCTransform(method = "qpoisson") %>% 
    RunPCA() %>% 
    RunUMAP(dims = seq(40)) %>%
    FindNeighbors(dims = seq(40)) %>%
    FindClusters(resolution = 0.8) -> seur
  return(seur)
}
