#' Imperative wrappers around SnaptronQueryBuilder class methods
#'
#' This wrappers are intended for users not familiar with the R6
#' Object Oriented Programming  system.
#' @export
QueryBuilder <- function(compilation, regions) {
    SnaptronQueryBuilder$new(compilation = compilation, regions = regions)
}

#' @export
get_compilation <- function(qb) {
    qb$compilation()
}

#' @export
get_regions <- function(qb) {
    qb$regions()
}

#' @export
get_column_filters <- function(qb) {
    qb$column_filters()
}

#' @export
get_row_filters <- function(qb) {
    qb$row_filters()
}

#' @export
get_sids <- function(qb) {
    qb$sids()
}

#' @export
get_coordinate_modifier <- function(qb) {
    qb$coordinate_modifier()
}

#' @export
set_compilation <- function(qb, compilation) {
    cloned_qb <- qb$clone(deep = TRUE)$compilation(compilation)
    cloned_qb
}

#' @export
set_regions <- function(qb, regions) {
    cloned_qb <- qb$clone(deep = TRUE)$regions(regions)
    cloned_qb
}

#' @export
set_column_filters <- function(qb, ...) {
    cloned_qb <- qb$clone(deep = TRUE)$column_filters(...)
    cloned_qb
}

#' @export
set_row_filters <- function(qb, ...) {
    cloned_qb <- qb$clone(deep = TRUE)$row_filters(...)
    cloned_qb
}

#' @export
set_sids <- function(qb, sids) {
    cloned_qb <- qb$clone(deep = TRUE)$sids(sids)
    cloned_qb
}

#' @export
set_coordinate_modifier <- function(coordinate_modifier) {
    cloned_qb <- qb$clone(deep = TRUE)$coodinate_modifier(coordinate_modifier)
    cloned_qb
}

#' @export
from_url <- function(url) {
    SnaptronQueryBuilder$new()$from_url(url)
}
