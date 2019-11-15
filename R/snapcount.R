#' snapcount: an R package for interfacing with Snaptron
#'
#' snapcount is a client interface to the Snaptron webservice which supports
#' querying by gene name or genomic region.
#'
#' Results include raw expression counts derived from alignment of RNA-seq
#' samples and/or various summarized measures of expression across one or more
#' regions/genes per-sample (e.g. percent spliced in).
#'
#' To learn more about snapcount, check out the vignette:
#' \code{browseVignettes(package = "snapcount")}
#'
#' @section Package options:
#' \describe{
#' \item{\code{snapcount.host}}{Change the host that snapcount uses when
#'   connecting to Snaptron. Default: \code{snaptron.cs.jhu.edu}}
#' \item{\code{snapcount.port}}{Change the port that snapcount uses when
#'   connecting to Snaptron. Default: \code{80}}
#' }
#'
#' @import data.table
#' @importFrom assertthat assert_that
#' @importFrom rlang !!
#' @importFrom magrittr %>%
#' @importFrom stats median
#' @importFrom R6 R6Class
#' @importFrom GenomicRanges GRanges
#' @importFrom SummarizedExperiment SummarizedExperiment
"_PACKAGE"
