# nocov start
#' Construct a QueryBuilder object given a compilation and one or regions.
#'
#' @param compilation A single string containing the name of the
#'   Snaptron data source. Any variant of the `Compilation` enum can also be
#'   passed an argument.
#' @param regions Either a list of 1 more `HUGO` gene names as strings
#'   `e.g. "BRCA1"` or a Granges class object containing one or more geneomic
#'   intervals `(e.g. "chr1:1-1000")`.
#'
#' @return
#' A QueryBuilder object.
#' @examples
#' # contruct a query builder for GTEX data source and BRAC1 gene
#' qb <- QueryBuilder(compilation = Compilation$gtex, regions = "BRCA1")
#'
#' # contruct a query builder for TCGA data source and chromosome region
#' qb <- QueryBuilder(compilation = "tcga", regions = "chr1:1-1000")
#'
#' # construct a query builder for TCGA data source using GRanges object
#' library(GenomicRanges)
#' qb <- QueryBuilder(compilation = "tcga", regions = GRanges("chr1", "1-1000"))
#' @export
QueryBuilder <- function(compilation, regions) {
    SnaptronQueryBuilder$new(compilation = compilation, regions = regions)
}

#' Get or set query compilation
#'
#' @param qb A QueryBuilder object constructed using the
#'   \code{\link{QueryBuilder}} function.
#' @param compilation A single string containing the name of the Snaptron data
#'   source. Any variant of the `Compilation` enum can also be passed as an
#'   argument.
#'
#' @return
#' \code{get_compilation} returns the current compilation as string.
#' \code{set_compilation} returns a new \code{QueryBuilder} object with
#'   the compilation set to the value of \code{compilation}.
#'
#' @examples
#' qb <- QueryBuilder(compilation = "gtex", regions = "CD99")
#' get_compilation(qb)
#' qb <- set_compilation(qb, Compilation$tcga)
#' get_compilation(qb)
#' @export
get_compilation <- function(qb) {
    qb$compilation()
}

#' Get or set query regions
#'
#' @param qb A QueryBuilder object constructed using the
#'   \code{\link{QueryBuilder}} function.
#' @param regions Either a list of 1 more `HUGO` gene names as strings
#'   `e.g. "BRCA1"` or a Granges class object containing one or more geneomic
#'   intervals `(e.g. "chr1:1-1000")`.
#'
#' @return
#' \code{get_regions} returns the current regions as a list of strings.
#' \code{set_regions} returns a new \code{QueryBuilder} object with the
#'   regions set to the value of \code{regions}.
#'
#' @examples
#' qb <- QueryBuilder(compilation = "gtex", regions = "CD99")
#' get_regions(qb)
#' qb <- set_regions(qb, "chr1:1-1000")
#' get_regions(qb)
#' qb <- set_regions(qb, GenomicRanges::GRanges("chr1", "1-1000"))
#' get_regions(qb)
#' @export
get_regions <- function(qb) {
    qb$regions()
}

#' Get or set sample-related contraints for query
#'
#' @param qb a QueryBuilder object constructed using the
#'   \code{\link{QueryBuilder}} function
#' @param ... one or more boolean predicates as either strings or unevaluated
#'   expressions
#'
#' @return
#' \code{get_column_filters} returns the current filters as a list of strings.
#' \code{set_column_filters} returns a new \code{QueryBuilder} object
#'   with the column filters set to the value of \code{column_filters}.
#'
#' @examples
#' qb <- QueryBuilder(compilation = "gtex", regions = "CD99")
#' # column filters set using a string
#' qb <- set_column_filters(qb, "SMTS == Brain")
#' get_column_filters(qb)
#' # column filters set using unevaluated expression
#' qb <- set_column_filters(qb, SMTS == "Spleen")
#' get_column_filters(qb)
#' @export
get_column_filters <- function(qb) {
    qb$column_filters()
}

#' Get or set range-related contraints for query
#'
#' @param qb a QueryBuilder object constructed using the
#'   \code{\link{QueryBuilder}} function.
#' @param ... one or more boolean predicates as either strings or unevaluated
#'   expressions.
#'
#' @return
#' \code{get_row_filters} returns the current row filters as list of strings.
#' \code{set_row_filters} returns a new \code{QueryBuilder} object with
#'   the row filters set to the value of \code{row_filters}.
#'
#' @examples
#' qb <- QueryBuilder(compilation = "gtex", regions = "CD99")
#' # row filters set as a string
#' qb <- set_row_filters(qb, "strand == +")
#' get_row_filters(qb)
#' # row filters set using unevaluated expression
#' qb <- set_row_filters(qb, strand == "+")
#' get_row_filters(qb)
#' @export
get_row_filters <- function(qb) {
    qb$row_filters()
}

#' Get or set query sample ids
#'
#' @param qb a QueryBuilder object constructed using the
#'   \code{\link{QueryBuilder}} function.
#' @param sids a vector or 1 or more whole numbers to filter results on.
#'
#' @return
#' \code{get_sids} returns the current sample ids as a vector of integers.
#' \code{set_sids} returns a new \code{QueryBuilder} object with the
#'   sample ids set to the value of \code{sids}.
#'
#' @examples
#' qb <- QueryBuilder(compilation = "gtex", regions = "CD99")
#' qb <- set_sids(qb, c(1, 2, 3))
#' get_sids(qb)
#' @export
get_sids <- function(qb) {
    qb$sids()
}

#' Get or set coordinate modifiers for the query.
#'
#' @param qb a QueryBuilder object constructed using the
#'   \code{\link{QueryBuilder}} function.
#' @param coordinate_modifier any of the variants of the \link{Coordinates}
#'   enum.
#'
#' @return
#' \code{get_coordinate_modifier} returns the current coodinate modifier
#'    as a string.
#' \code{set_coordinate_modifier} returns a new \code{QueryBuilder}
#'   object with the coordinate modifier set to the value of
#'   \code{coordinate_modifier}.
#'
#' @examples
#' qb <- QueryBuilder(compilation = "gtex", regions = "CD99")
#' qb <- set_coordinate_modifier(qb, Coordinates$Within)
#' get_coordinate_modifier(qb)
#' @export
get_coordinate_modifier <- function(qb) {
    qb$coordinate_modifier()
}

#' @rdname get_compilation
#' @export
set_compilation <- function(qb, compilation) {
    cloned_qb <- qb$clone(deep = TRUE)$compilation(compilation)
    cloned_qb
}

#' @rdname get_regions
#' @export
set_regions <- function(qb, regions) {
    cloned_qb <- qb$clone(deep = TRUE)$regions(regions)
    cloned_qb
}

#' @rdname get_column_filters
#' @export
set_column_filters <- function(qb, ...) {
    quosures <- rlang::enquos(...)
    cloned_qb <- qb$clone(deep = TRUE)$column_filters(quosures)
    cloned_qb
}

#' @rdname get_row_filters
#' @export
set_row_filters <- function(qb, ...) {
    quosures <- rlang::enquos(...)
    cloned_qb <- qb$clone(deep = TRUE)$row_filters(quosures)
    cloned_qb
}

#' @rdname get_sids
#' @export
set_sids <- function(qb, sids) {
    cloned_qb <- qb$clone(deep = TRUE)$sids(sids)
    cloned_qb
}

#' @rdname get_coordinate_modifier
#' @export
set_coordinate_modifier <- function(qb, coordinate_modifier) {
    cloned_qb <- qb$clone(deep = TRUE)$coordinate_modifier(coordinate_modifier)
    cloned_qb
}

#' Constructs a \code{QueryBuilder} object from the given url
#'
#' @param url a well-formed url preferably obtained from a call to the
#'   \code{\link{uri_of_last_successful_request}} function
#'
#' @return
#' Returns a \code{QueryBuilder} object with attributes set from the
#'   parsed url.
#'
#' @examples
#' sb <- from_url("http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99")
#' get_regions(sb)
#' get_compilation(sb)
#' @export
from_url <- function(url) {
    SnaptronQueryBuilder$new()$from_url(url)
}

#nocov end
