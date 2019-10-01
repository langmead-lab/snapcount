is_hugo_gene <- function(str) {
    assert_that(is.character(str))
    stringr::str_detect(str, "^[A-Za-z-0-9_]+@?$")
}

is_chromosome_interval <- function(str) {
    assert_that(is.character(str))
    stringr::str_detect(str, "chr[1-22XYM]:\\d+-\\d+")
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    assert_that(is.numeric(x))
    all(abs(x - round(x)) < tol)
}

is_query_builder <- function(object) {
    all(c("SnaptronQueryBuilder", "R6") == class(object))
}

is_list_of_query_builders <- function(list) {
    assert_that(is.list(list))
    purrr::reduce(lapply(list, is_query_builder), all)
}

is_list_of_query_builder_groups <- function(list_of_groups) {
    assert_that(is.list(list_of_groups))
    purrr::reduce(lapply(list_of_groups, is_list_of_query_builders), all)
}