# nocov start
enum <- function(...) {
    variants <- rlang::ensyms(...)
    values <- purrr::map(variants, rlang::as_string)

    res <- stats::setNames(values, values)
    res <- as.environment(as.list(res))

    lockEnvironment(res, bindings = TRUE)

    res
}
#' Enum for Snaptron Coordinate modifiers
#'
#' @field Exact Return junctions whose start and end
#'   coordinates match the boundaries of the
#'   region requested.
#' @field Within Return junctions whose start and end
#'   coordinates are within the boundaries of the
#'   region requested.
#' @field StartIsExactOrWithin Return junctions whose
#'   start coordinate matches, or is within, the
#'   boundaries of the region requested.
#' @field EndIsExactOrWithin Return junctions whose
#'   end coordinate matches, or is within, the
#'   boundaries of the region requested.
#'
#' @examples
#' qb <- QueryBuilder(compilation = "gtex", regions = "CD99")
#' qb <- set_coordinate_modifier(qb, Coordinates$Exact)
#' qb
#' @export
Coordinates <- enum(Exact, Within, StartIsExactOrWithin, EndIsExactOrWithin)

#' Enum for Snaptron compilations
#'
#' The variants for this enum will be populated dynamically after the
#' package has been loaded. If the package cannot connect to the internet
#' the variants will default to:
#'
#' \itemize{
#'   \item{gtex}
#'   \item{tcga}
#'   \item{srav2}
#'   \item{sra}
#' }
#'
#' @seealso \url{http://snaptron.cs.jhu.edu/data.html} for more information
#'   about Snaptron compilations.
#'
#' @examples
#' qb <- QueryBuilder(compilation = Compilation$gtex, regions = "KCNIP4")
#' query_jx(qb)
#' @export
Compilation <- enum(gtex, tcga, srav2, sra)

#nocov end
