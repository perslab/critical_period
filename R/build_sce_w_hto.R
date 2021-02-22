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

build_sce_w_hto <- function(exp_files, hto_files) {

  #get paths to directories
  exp_directory <- unlist(str_match_all(exp_files[[1]], "(.*/counts_unfiltered/)\\.*"))[2]
  hto_directory <- unlist(str_match_all(hto_files[[1]], "(.*/counts_unfiltered/)\\.*"))[2]
  
  # Read in kallisto counts
  ct_mat <- BUSpaRse::read_velocity_output(
    spliced_dir = exp_directory,
    spliced_name = "spliced",
    unspliced_dir = exp_directory,
    unspliced_name = "unspliced"
  )
  
  # Read HTO mat
  hto <- BUSpaRse::read_count_output(dir = hto_directory, name = "cells_x_features", tcc = FALSE)
  
  # find barcodes that overlap between spliced mat, unspliced mat, and hto mat
  bc_overlap <- intersect(colnames(ct_mat$spliced), colnames(ct_mat$unspliced)) %>% 
    intersect(., colnames(hto))
  
  # only keep genes which can map to busparse mapping
  # code to update transcript and genes -- BUSpaRse::tr2g_ensembl("Mus musculus")
  mm_gene <- readr::read_tsv("/nfsdata/projects/dylan/cp_targets/tr2g.tsv", col_names = F)
  colnames(mm_gene) <- c("transcript", "gene", "gene_name", "type")
  tot_genes <- Matrix::rowSums(ct_mat$unspliced)
  genes_use <- rownames(ct_mat$unspliced)[tot_genes > 1] %>% 
    intersect(mm_gene$gene)
  
  # Calculate barcode rank and classify cells
  tot_count <- Matrix::colSums(ct_mat$unspliced[,bc_overlap])
  bc_uns <- DropletUtils::barcodeRanks(ct_mat$unspliced[,bc_overlap])
  bcs_use <- names(tot_count)[tot_count > metadata(bc_uns)$inflection]
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
  
  # convert hash to seurat for demultiplexing of HTO using k-medoid clustering
  seur.hto <- as.Seurat(altExp(sce.full, "hto"), counts = "counts", data = NULL)  
  seur.hto <- NormalizeData(seur.hto,normalization.method = "CLR")
  # calculate quantile for use in Demux
  quant <- determine_quantile(seur.hto)
  print(quant[[2]])
  seur.hto <- HTODemux(seur.hto, positive.quantile = quant[[2]], assay = "RNA")
  print(table(seur.hto$RNA_classification.global))
  
  # add doublet information to SCE object
  sce.full$hto_global <- seur.hto$RNA_classification.global
  sce.full$hash_max <- seur.hto$RNA_maxID
  sce.full$hash_id <- seur.hto$hash.ID
  sce.full$doublet <- if_else(sce.full$hto_global == "Doublet", true = T, false = F)
  
  # identify intra-hash doublets
  library(scran)
  sce.full <- logNormCounts(sce.full)
  dec.hash <- modelGeneVar(sce.full)
  top.hash <- getTopHVGs(dec.hash, n=1000)
  set.seed(1011110)
  sce.full <- runPCA(sce.full, subset_row=top.hash, ncomponents=20)
  
  # Recovering the intra-sample doublets:
  hashed.doublets <- scDblFinder::recoverDoublets(sce.full, use.dimred="PCA",
                                                  doublets=sce.full$doublet, samples=table(sce.full$hash_max))
  
  sce.full$proportion_dub_neighbors <- hashed.doublets$proportion
  sce.full$predicted_dub_std <- hashed.doublets$predicted
  sce.full$predicted_dub_cut <- dub_cutoff(sce.full)
  
  # qc filtering
  sce.full <- addPerCellQC(sce.full)
  reason <- data.frame(
    "qc.umi.low" = isOutlier(sce.full$sum, log=TRUE, type="lower"),
    "qc.umi.high" = isOutlier(sce.full$sum, type="higher"),
    "qc.nexprs.low" = isOutlier(sce.full$detected, log=TRUE, type="lower"),
    "qc.nexprs.high" = isOutlier(sce.full$detected, type="higher")
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
