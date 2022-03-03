##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param exp_files batched list of files pointing to location of gene expression data
##' @param hto_files batched list of files pointing to location of hashtag data
##' @return sce object which has been demultiplexed, but retained all called cells, qc metrics are in metadata
##' @author dylanmr
##' @export
##' 


.dub_cutoff <- function(x) {
  
  # calculate number of neighbors at each proportion that are doublets
  data.frame("prop"=x$proportion_dub_neighbors) %>% 
    group_by(prop) %>% 
    summarize(n=n()) %>% 
    mutate(pct = n/sum(n)) -> data 
  # find point at which we gain very few doublets as proportion increases 
  cut <- data$prop[PCAtools::findElbowPoint(variance = sort(data$n, decreasing = T))+1]
  vec <- if_else(x$proportion_dub_neighbors <= cut, F, T)
  return(vec)
  
}

build_sce_w_hto <- function(exp_files, hto_files, seed=426) {
  
  set.seed(seed)
  
  run <- unlist(str_match_all(exp_files, ".*kallisto/(.*)\\/counts.*"))[2]
  
  # Read in kallisto counts
  ct_mat <- BUSpaRse::read_velocity_output(
    spliced_dir = exp_files,
    spliced_name = "spliced",
    unspliced_dir = exp_files,
    unspliced_name = "unspliced"
  )
  
  rownames(ct_mat$unspliced) <-  rownames(ct_mat$unspliced) %>% str_split_fixed(., "[.]", n=2) %>% data.frame() %>% pull(1)
  rownames(ct_mat$spliced) <-  rownames(ct_mat$spliced) %>% str_split_fixed(., "[.]", n=2) %>% data.frame() %>% pull(1)
  
  # Read HTO mat
  hto <- BUSpaRse::read_count_output(
    dir = hto_files, 
    name = "cells_x_features", 
    tcc = FALSE
  )
  
  # find barcodes that overlap between spliced mat, unspliced mat, and hto mat
  bc_overlap <- intersect(colnames(ct_mat$spliced), 
                          colnames(ct_mat$unspliced)
                          ) %>% 
    intersect(., colnames(hto))
  
  # only keep genes which can map to busparse mapping
  # code to update transcript and genes -- BUSpaRse::tr2g_ensembl("Mus musculus")
  mm_gene <- readr::read_tsv("/nfsdata/projects/dylan/cp_targets/tr2g.tsv", col_names = F) %>% 
    set_colnames(c("transcript", "gene", "gene_name", "type")) %>% 
    separate(gene, into = c("gene", NA)) %>% distinct(gene, .keep_all=T)
  tot_genes <- Matrix::rowSums(ct_mat$unspliced)
  genes_use <- rownames(ct_mat$unspliced)[tot_genes > 1] %>% 
    intersect(mm_gene$gene)
  
  # Calculate barcode rank and classify cells
  tot_count <- Matrix::colSums(ct_mat$unspliced[,bc_overlap])
  bc_uns <- DropletUtils::barcodeRanks(ct_mat$unspliced[,bc_overlap])
  bcs_use <- names(tot_count)[tot_count > metadata(bc_uns)$inflection]
  if(length(bcs_use) > 2e4) {
    bcs_use <- sort(tot_count, decreasing=T)[1:8e3] %>% names()
  }
  print(paste0(length(bcs_use), " cells have been detected"))

  # build SCE objects for all assays
  sce.spliced <- SingleCellExperiment(list(counts = ct_mat$spliced[genes_use,bcs_use]))
  sce.unspliced <- SingleCellExperiment(list(counts = ct_mat$unspliced[genes_use,bcs_use]))
  sce.hash <- SingleCellExperiment(list(counts = hto[,bcs_use]))
  sce.full <- SingleCellExperiment(list(counts = ct_mat$spliced[genes_use,bcs_use]+ct_mat$unspliced[genes_use,bcs_use]))
  
  # pull down matching gene and transcript names
  rowData(sce.full)$gene_name <-  rownames(sce.full) %>% 
    data.frame("ens"=.) %>% 
    inner_join(mm_gene %>% 
                 distinct(gene, .keep_all=T),
               by = c("ens"="gene")) %>% 
    pull(gene_name)
  
  # add assays as alternative experiments to SCE object
  altExp(sce.full, "hto") <- sce.hash
  altExp(sce.full, "spliced") <- sce.spliced
  altExp(sce.full, "unspliced") <- sce.unspliced
  
  # calculate quantile for use in HTODemux
  quant <- autothreshold(sce.full, run=run)
  
  # convert hash to seurat for demultiplexing of HTO using k-medoid clustering
  seur.hto <- CreateSeuratObject(counts = counts(altExp(sce.full, "hto")))
  seur.hto <- subset(seur.hto, features = rownames(seur.hto)[rowMeans(seur.hto)>1])
  seur.hto <- NormalizeData(seur.hto,normalization.method = "CLR")
  seur.hto <- HTODemux(seur.hto, positive.quantile = quant, verbose = 0, assay="RNA")
  print(table(seur.hto$RNA_classification.global))
  
  # add doublet information to SCE object
  sce.full$hto_global <- seur.hto$RNA_classification.global
  sce.full$hash_max <- seur.hto$RNA_maxID
  sce.full$hash_id <- seur.hto$hash.ID
  sce.full$doublet <- if_else(sce.full$hto_global == "Doublet", true = T, false = F)
  
  # identify intra-hash doublets
  library(scran)
  sce.full <- scuttle::logNormCounts(sce.full)
  dec.hash <- scran::modelGeneVar(sce.full)
  top.hash <- scran::getTopHVGs(dec.hash, n=1000)
  sce.full <- scater::runPCA(sce.full, subset_row=top.hash, ncomponents=20)
  
  # Recovering the intra-sample doublets:
  hashed.doublets <- scDblFinder::recoverDoublets(sce.full, 
                                                  use.dimred="PCA",
                                                  doublets=sce.full$doublet, 
                                                  samples=table(sce.full$hash_max)
                                                  )
  
  sce.full$proportion_dub_neighbors <- hashed.doublets$proportion
  sce.full$predicted_dub_std <- hashed.doublets$predicted
  sce.full$predicted_dub_cut <- .dub_cutoff(sce.full)
  
  # qc filtering
  sce.full <- scuttle::addPerCellQC(sce.full)
  reason <- data.frame(
    "qc.umi.low" = scuttle::isOutlier(sce.full$sum, log=TRUE, type="lower"),
    "qc.umi.high" = scuttle::isOutlier(sce.full$sum, type="higher"),
    "qc.nexprs.low" = scuttle::isOutlier(sce.full$detected, log=TRUE, type="lower"),
    "qc.nexprs.high" = scuttle::isOutlier(sce.full$detected, type="higher")
  ) %>% 
    mutate(reason = pmap_int(.,~any(.==T)),
           reason = if_else(reason == 1, true=T, false=F)) %>% 
    dplyr::pull(reason)
  
  # annotate object to discard
  # will discard in later steps, build rmarkdown file to visualize all QC steps
  sce.full$discard <- reason
  
  # object to return
  return(sce.full)

}
