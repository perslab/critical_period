##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param seq_stats
##' @param cell_stats
##' @param dir
##' @param sample_id
##' @return
##' @author sarah145 https://github.com/Sarah145/scRNA_pre_process/blob/master/scripts/functions.R
##' @export

print_HTML <- function(seq_stats, cell_stats, dir , sample_id) {
    
    library(R2HTML)
  
    system(paste0('base64 ', dir, sample_id, '_barcode_rank.png > ', dir, sample_id, '_barcode_rank.txt'))
    b64_bc <- readChar(paste0(dir, sample_id, '_barcode_rank.txt'), file.info(paste0(dir, sample_id, '_barcode_rank.txt'))$size)
    target <- R2HTML::HTMLInitFile(dir, filename=paste0(sample_id, '_summary'))
    HTML('<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Roboto">', file=target)
    HTML("<div class='title'>", file=target)
    HTML.title(' Pre-Processing Summary', HR=1, file = target)
    HTML("</div>", file = target)
    HTML.title(sample_id, HR=2, file = target)
    HTML("<div id='wrapper'>", file=target)
    HTML("<div class='boxed' id='left' align='center'>", file=target)
    HTML.title('Sequencing/Alignment Stats', HR=3, file=target)
    HTML('<table style="width:100%">', file=target)
    HTML(paste('<tr> <td>', seq_stats$stat, '</td> <td align="right">', seq_stats$value, '</td> </tr>'), file=target)
    HTML('</table> <hr>', file=target)
    HTML.title('Cell Stats', HR=3, file=target)
    HTML('<table style="width:100%">', file=target)
    HTML(paste('<tr> <td>', cell_stats$stat, '</td> <td align="right">', cell_stats$value, '</td> </tr>'), file=target)
    HTML('</table>', file=target)
    HTML("</div>", file = target)
    HTML("<div class='boxed' id='right' align='center'>", file=target)
    HTML(paste0("<img src='data:image/png;base64,", b64_bc, "' width=90%>"), file=target)
    HTML("</div>", file = target)
    HTML("</div>", file = target)
    HTML('<style type="text/css">
         .title {
         background-color: #0972D5;
         padding: 8px;
         color: white;
         position: fixed;
         top: 0;
         left: 0;
         z-index: 999;
         width: 100%;
         }
         .boxed {
         border: 1px solid #868D96;
         padding: 10px;
         margin: 20px;
         }
         h1 {
         font-family: "Roboto";
         font-size: 33px;
         }
         h2 {
         font-family: "Roboto";
         font-size: 26px;
         }
         h3 {
         font-family: "Roboto";
         font-size: 18px;
         }
         #wrapper {
         display: flex;
  }
         #left {
         width: 50%;
         }
         #right {
         width: 50%;
         }
         table tr:nth-child(even) {
         background-color: #eee;
         }
         table tr:nth-child(odd) {
         background-color: #fff;
         }
         table {
         font-family: "Roboto";
         font-size: 20px;
         border: 1px solid #868D96;
         }
         #mathplayer{
         height: 80px;
         }
         </style> </head>', file=target)
	HTMLEndFile()
}


