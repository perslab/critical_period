library(targets)
library(tarchetypes)
library(tidyverse)

options(tidyverse.quiet = TRUE)
future::plan(future::multisession)

tar_option_set(packages = c(
  "mclust", "Seurat", "tidyverse", "tidymodels", "yardstick", "R2HTML",
  "DropletUtils", "Matrix", "rjson", "googlesheets4", "stringr", "future",
  "org.Mm.eg.db", "scDblFinder", "AnnotationDbi", "SingleCellExperiment",
  "scran", "scater", "glmGamPoi", "miloR", "BiocNeighbors", "data.table")
)

# define parameters for static branching for clustering

clustering_params <- tidyr::expand_grid(
  dims = seq(20, 50, by = 10),
  k = seq(5, 30, by = 5),
  weight = c("jaccard", "rank"),
  clus_method = c("louvain", "walktrap")
) %>%
  unite("name", clus_method, weight, dims, k, remove = F)

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

# ensure metadata is up to date
tar_pipeline(
  tar_target(exp_paths, path_to_data(), iteration = "list"),
  tar_target(exp_files, exp_paths, format = "file", pattern = map(exp_paths)),
  tar_target(hto_paths, path_to_data(type = "hto"), iteration = "list"),
  tar_target(hto_files, hto_paths, format = "file", pattern = map(hto_paths)),
  tar_target(web_summaries, gen_summary(files = exp_files), format = "file", pattern = map(exp_files)),
  tar_target(sce_objects, build_sce_w_hto(exp_files, hto_files), pattern = map(exp_files, hto_files), iteration = "list"),
  tar_render_rep(
    report,
    "qc_report.Rmd",
    params = tibble(
      par = c(1:29),
      output_file = paste0("sc_qc_", 1:29, ".html"))
  ),
  tar_target(
    sc_comb,
    comb_and_filter(sce_objects) %>% 
      sctClust(dims=40) %>% 
      mapCamp() %>% 
      classifySex()
  ),
  tar_target(
    split_seur,
    list(subset(sc_comb, subset=predicted.id=="neuron") %>% 
           sctClust(dims=40) %>% 
           mapCamp(class="neur"), 
         subset(sc_comb, subset=predicted.id!="neuron") %>% 
           sctClust(dims=40) %>% 
           mapCamp())
  ),
  tar_target(
    velo_obj,
    comb_and_filter(sce_objects, keep_alt=T)
  ),
  # tar_target(
  #   milo, 
  #   gen_milo(x = split_seur, k=30, d=30),
  #   pattern = map(split_seur), 
  #   iteration="list"
  # ),
  tar_target(
    pseudobulk,
    Seurat::as.SingleCellExperiment(split_seur[[1]], assay="SCT") %>%
      scater::aggregateAcrossCells(., use.assay.type="counts", id=colData(.)[,c("predicted.id", "hash_id")]),
    pattern = map(split_seur), 
    iteration="list"
  ),
  tar_target(
    knn_graph,
    build_graph(split_seur[[1]]),
    pattern = map(split_seur), 
    iteration="list"
  ),
  tar_target(
    resampled_graph,
    knn_resampling(knn_graph, k = 30, meta = split_seur[[1]][[]]),
    pattern = map(knn_graph, split_seur),
    iteration = "list"
  )
  # tar_target(
  #   integ_with_embr,
  #   build_ds("/data/pub-others/huisman-natcom-2019/GSE126480_Raw_counts.txt.gz") %>% 
  #     integrate_data(sc_comb, ., label = "huisman")
  # )
  #),
  # tar_map(
  #   values = clustering_params,
  #   names = "name", # Select columns from `values` for target names.
  #   tar_target(
  #     clusters,
  #     cluster_across_params(sce_final, dims = dims, k = k, weight = weight, clus_method = clus_method)
  #   )
  #)
)
