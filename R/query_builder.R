#' A Reference Class for building Snaptron queries
#'
#' @field compilation A single string containing the name of the Snaptron datasource
#'
#' @field regions Either a list of 1 or more `HUGO` gene names `(e.g. "BRCA1")` or a
#'   GRanges-class object containing one or more genomic intervals `(e.g. "chr1:1-1000")`.
#'   Strand information is ignored.
#'
#' @field row_filters A list of strings defining range-related contraints
#'
#' @field column_filters A list of strings defining sample-related contraints
#'
#' @field sids A list of rail_ids (integer sample IDs) to filter results on. Only
#'   records which have been found in at least one of these samples will be returned.
#'
#' @field coordinate_modifier Snaptron coordinate modifier enum. Variants include:
#'
#'   Coordinate$Exact - Consrains the results so that the start/end coordinates
#'              match the start/end of the specified range.
#'
#'   Coordinate$Within - Constrains the results so that the coordinates are
#'              within (inclusive) the specified range.
#'
#'   Coordinate$StartIsExactOrWithin - Constrains the results so that the start
#'              coorindate matches or is within the boundaries of the specified range.
#'
#'   Coorindate$EndIsExactOrWithin - Constrains the results so that the end
#'              coordinate matches or is within the boundaries of the specified range.
#'
#'   Coordinate$Overlaps - Constrains the results so that the coordinates overlap the
#'              specified range.
#'
#' @section Methods:
#' \itemize{
#'   \item{\code{compilation} Get or set query compilation. See
#'   \code{\link{Compilation}} for details.}
#'
#'   \item{\code{regions} Get or set query region which can be either a chromosome
#'   range or HUGO gene name.}
#'
#'   \item{\code{row_filters} Get or set query range filters. Range filters can be
#'   set to a list or vector of either strings or R boolean expressions. If
#'   filters are R boolean expressions then the right-hand-side of the R
#'   expression should evaluate to either a string, float or whole number.
#'   The expression will then be converted to a string before being sent out
#'   to Snaptron.}
#'
#'   \item{\code{column_filters} Get or set query sample filters. When setting sample
#'   filters the argument is expected to follow the sample rules as described
#'   for \code{row_filters}.}
#'
#'   \item{\code{sids} Get or set query sample ids.}
#'
#'   \item{\code{coordinate_modifier} Get or set coordinate modifiers for query. See
#'   \code{\link{Coordinates}} for details.}
#'
#'   \item{\code{query_jx} call \code{\link{query_jx}} function}
#'
#'   \item{\code{query_gene} call \code{\link{query_gene}} function}
#'
#'   \item{\code{query_exon} call \code{\link{query_exon}} function}
#'
#'   \item{\code{print} print query builder object}
#'
#'   \item{\code{from_url} use URL to instantiate QueryBuilder object }
#' }
#'
#' @examples
#' sb <- SnaptronQueryBuilder$new()
#' sb$compilation("gtex")$regions("CD99")$query_jx()
#' sb$from_url("http://snaptron.cs.jhu.edu/gtex/snaptron?regions=chr1:1-100000")
#' sb$column_filters(SMTS == "Brain")
#'
#' # The RHS of an R boolean expression can itself be an expression,
#' # but the expression has to evaluate to either a String, Float or Integer.
#' a <- 10
#' b <- 10
#' sb$row_filters(samples_count >= (a + b))
#' sb$query_jx(return_rse = FALSE)
#' @export
SnaptronQueryBuilder <- R6Class("SnaptronQueryBuilder",
    public = list(
        initialize = function(...) {
            private$query <- list(...)
        },
        compilation = function(compilation = NULL) {
            if (!missing(compilation)) {
                assert_that(compilation %in% names(Compilation),
                            msg = paste0(compilation, ": is not a valid compilation"))
                private$query$compilation <- compilation
                invisible(self)
            } else {
                private$query$compilation
            }
        },
        regions = function(regions = NULL) {
            if (!missing(regions)) {
                if (is(regions, "GRanges")) {
                    private$query$regions <- regions
                } else if (is_genes_or_intervals(regions)) {
                    private$query$regions <- regions
                } else {
                    message <-
                        paste("regions must contain strings representing",
                              "HUGO genes or chromosome intervals, of the form",
                              "chr:start-end or chr:start-end:strand,",
                              "or GRanges object.")
                    private$query$regions <-
                        tryCatch(GenomicRanges::GRanges(regions),
                                 error = function(e) {
                                     stop(message)
                                 })
                }
                invisible(self)
            } else {
                private$query$regions
            }
        },
        row_filters = function(...) {
            if (!missing(...)) {
                private$query$row_filters <-
                    bool_expressions_to_strings(rlang::enquos(...))
                invisible(self)
            } else {
                private$query$row_filters
            }
        },
        column_filters = function(...) {
            if (!missing(...)) {
                private$query$column_filters <-
                    bool_expressions_to_strings(rlang::enquos(...))
                invisible(self)
            } else {
                private$query$column_filters
            }
        },
        sids = function(sids = NULL) {
            if (!missing(sids)) {
                assert_that(is.wholenumber(sids),
                    msg = "sids should be whole numbers")
                private$query$sids <- sids
                invisible(self)
            } else {
                private$query$sids
            }
        },
        coordinate_modifier = function(coordinate_modifier = NULL) {
            if (!missing(coordinate_modifier)) {
                private$query$coordinate_modifier <- coordinate_modifier
                invisible(self)
            } else {
                private$query$coordinate_modifier
            }
        },
        query_jx = function(return_rse = TRUE) {
            private$call("query_jx", list(sb = self, return_rse = return_rse))
        },
        query_exon = function(return_rse = TRUE) {
            private$call("query_exon", list(sb = self, return_rse = return_rse))
        },
        query_gene = function(return_rse = TRUE) {
            private$call("query_gene", list(sb = self, return_rse = return_rse))
        },
        from_url = function(url) {
            url <- httr::parse_url(url)
            if (url$hostname != "snaptron.cs.jhu.edu") {
                stop("URL does not point to Snaptron server", stop. = FALSE)
            }
            resp <- httr::HEAD(url)
            if (resp$status_code != 200 ||
                httr::http_type(resp) != "text/plain") {
                stop(sprintf("%s: is not a valid URL", url), call. = FALSE)
            }
            query <- list()
            for (i in seq_along(url$query)) {
                name <- switch(n <- names(url$query[i]),
                    rfilter = "row_filters",
                    sfilter = "column_filters",
                    regions = "regions",
                    n)

                if (name == "sids") {
                    query[[name]] <-
                        scan(textConnection(url$query[[i]]), sep = ",")
                } else if (name == "contains") {
                    if (url$query[[i]] == "1") {
                        query[["coordinate_modifier"]] <- Coordinates$Within
                    }
                } else if (name == "exact") {
                    if (url$query[[i]] == "1") {
                        query[["coordinate_modifier"]] <- Coordinates$Exact
                    }
                } else if (name == "either") {
                    if (url$query[[i]] == "1") {
                        query[["coordinate_modifier"]] <-
                            Coordinates$StartIsExactOrWithin
                    } else if (url$query[[i]] == "2") {
                        query[["coordinate_modifier"]] <-
                            Coordinates$EndIsExactOrWithin
                    }
                } else {
                    query[[name]] <- c(query[[name]], url$query[[i]])
                }
            }

            query$compilation <- strsplit(url$path, "/")[[1]][1]
            private$query <- query

            invisible(self)
        },
        print = function() {
            cat("<SnaptronQueryBuilder>\n")
            for (param in names(private$query)) {
                if (is.null(private$query[[param]])) {
                    next
                } else if (param == "coordinate_modifier") {
                    desc <- switch(private$query[[param]],
                                   Exact = "exact",
                                   Within = "contains",
                                   StartIsExactOrWithin = "either=1",
                                   EndIsExactOrWithin = "either=2",
                                   "overlaps")
                } else {
                    desc <- paste(private$query[[param]], collapse = ",")
                }
                cat("   ", param, ": ", desc, "\n", sep = "")
            }
        }
    ),
    private = list(
        query = list(),
        call = function(fn_name, args) {
            if (is.null(self$compilation())) {
                stop(
                    paste("Please set a compilation before running", fn_name),
                    call. = FALSE
                )
            }
            if (is.null(self$regions())) {
                stop(
                    paste("Please specify query regions before running", fn_name),
                    call. = FALSE
                )
            }
            fn <- get(fn_name, parent.frame())
            do.call(fn, args)
        }
    )
)
