##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param files file list of expression data in default kallisto format
##' @param output_folder name of output folder
##' @return web summary file with QC info about the quality of the run
##' @author dylanmr
##' @import DropletUtils R2HTML rjson ggplot2 Matrix
##' @export

source("R/print_HTML.R")

gen_summary <- function(files, output_folder = "output") {
  
  # create output folder if does not exist, will spit out warning but continue if folder already exists
  dir.create(file.path(here::here(output_folder), "web_summary"))
  # get hash-pool name from expression file name
  run <- unlist(str_match_all(files, ".*kallisto/(.*)\\/counts.*"))[2]
  print(run)
  # load raw mtx
  raw_mtx <- Matrix::readMM(paste0(files, "/unspliced.mtx"))
  # load genes
  genes <- data.table::fread(paste0(files, "/unspliced.genes.txt"), header = F)[, 1]
  colnames(raw_mtx) <- genes[["V1"]]
  # attach barcodes
  rownames(raw_mtx) <- data.table::fread(paste0(files, "/unspliced.barcodes.txt"), header = F)[["V1"]]
  # transpose matrix for downstream steps
  raw_mtx <- Matrix::t(raw_mtx)
  # get probability that each barcode is a cell from the dropletUtils package (could probably be updated)
  tot_count <- Matrix::colSums(raw_mtx)
  bc_uns <- DropletUtils::barcodeRanks(raw_mtx)
  bcs_use <- names(tot_count)[tot_count > bc_uns@metadata$inflection]
  # subset raw mtx to remove empty drops
  filt_mtx <- raw_mtx[, bcs_use]
  # load run info
  kb_stats <- c(
    rjson::fromJSON(file = paste0(unlist(str_match_all(files[[1]], "(.*)\\/counts_unfiltered.*"))[2], "/inspect.json")),
    rjson::fromJSON(file = paste0(unlist(str_match_all(files[[1]], "(.*)\\/counts_unfiltered.*"))[2], "/run_info.json"))
  )
  # determine chemistry version
  tech <- grep("10x(.*)", strsplit(kb_stats$call, "\\s")[[1]][8], value = T)
  # get sequencing/alignment stats
  seq_stats <- data.frame(
    stat = c(
      "Sequencing technology", "Number of reads processed", "% reads pseudoaligned",
      "% reads on whitelist"
    ),
    value = prettyNum(c(
      tech, kb_stats$n_processed, kb_stats$p_pseudoaligned,
      round(kb_stats$percentageReadsOnWhitelist, 2)
    ), big.mark = ",")
  )
  # calculate cell stats and save to df
  p_cnts_in_cells <- round((sum(filt_mtx) / sum(raw_mtx)) * 100, 2)
  med_cnts_cell <- median(Matrix::colSums(filt_mtx))
  med_genes_cell <- median(apply(filt_mtx, 2, function(x) sum(x >= 1)))
  tot_genes_detected <- sum(Matrix::rowSums(filt_mtx) >= 1)
  cell_stats <- data.frame(
    stat = c(
      "Estimated number of cells", "% counts in cells",
      "Median counts per cell", "Median genes per cell", "Total genes detected"
    ),
    value = prettyNum(c(
      ncol(filt_mtx), p_cnts_in_cells, med_cnts_cell,
      med_genes_cell, tot_genes_detected
    ), big.mark = ",")
  )

  # get rank stats
  stats <- DropletUtils::barcodeRanks(raw_mtx)
  # create barcode rank plot png
  plain <- function(x, ...) {
    format(x, ..., scientific = FALSE, drop0trailing = TRUE, big.mark = ",")
  }

  cells <- colnames(raw_mtx) %in% colnames(filt_mtx)
  keep <- !duplicated(stats$total)
  plot_df <- data.frame(Rx = stats$rank, Tx = stats$total, cell = cells)
  plot_df <- plot_df[keep, ]
  bc_plot <- ggplot(subset(plot_df, plot_df$Tx > 0), aes(x = Rx, y = Tx, col = cell, alpha = cell)) +
    geom_point(size = 3) +
    geom_hline(yintercept = stats@metadata$knee, lty = 2, col = "#0972D5", size = 1.5) +
    annotate("text", x = max(plot_df$Tx), y = stats@metadata$knee + 10000, label = "Knee", color = "#0972D5", size = 5.5) +
    geom_hline(yintercept = stats@metadata$inflection, lty = 2, col = "#09CFD5", size = 1.5) +
    annotate("text", x = max(plot_df$Tx), y = stats@metadata$inflection + 500, label = "Inflection", color = "#09CFD5", size = 5.5) +
    scale_x_log10(labels = plain, breaks = scales::trans_breaks("log10", function(x) round(10^x, 0))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) floor(10^x)), labels = plain) +
    scale_color_manual(values = c("#8595A8", "#6406B6"), name = NULL, labels = c("Background", "Cells")) +
    scale_alpha_manual(values = c(0.5, 1)) +
    labs(x = "Barcodes", y = "UMI counts", title = "Barcode Rank Plot") +
    guides(alpha = FALSE, colour = guide_legend(reverse = TRUE, override.aes = list(size = 5))) +
    theme_linedraw() +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 15),
      legend.text = element_text(size = 19),
      legend.background = element_rect(fill = "transparent"),
      legend.position = c(0.15, 0.15)
    )

  # save raw plot to output folder
  ggsave(here::here(output_folder, "web_summary", paste0(run, "_barcode_rank.png")), bc_plot, width = 7, height = 4)

  # output a HTML summary of the run
  ## called from defined function in R folder
  print_HTML(seq_stats = seq_stats, cell_stats = cell_stats, dir = here::here(output_folder, "web_summary/"), sample_id = run)

  # return file path so that it can be tracked
  return(here::here(output_folder, "web_summary", paste0(run, "_summary.html")))
}
