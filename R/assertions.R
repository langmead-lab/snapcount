is_hugo_gene <- function(str) {
    assert_that(is.character(str))
    stringr::str_detect(str, "^[A-Za-z0-9_]+@?$")
}

is_genes_or_intervals <- function(strings) {
    bools <-
        Map(function(s) is_hugo_gene(s) || is_chromosome_interval(s), strings)
    Reduce(`&&`, bools, TRUE)
}

is_chromosome_interval <- function(str) {
    assert_that(is.character(str))
    stringr::str_detect(str, "chr([1-9]|1[0-9]|2[0-2]|[XYM]):\\d+-\\d+")
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    assert_that(is.numeric(x))
    all(x >= 0 & (abs(x - round(x)) < tol))
}

is_query_builder <- function(object) {
    is(object, "SnaptronQueryBuilder") && is(object, "R6")
}

is_list_of_query_builders <- function(list) {
    assert_that(is.list(list))
    purrr::reduce(lapply(list, is_query_builder), all)
}

is_list_of_query_builder_groups <- function(list_of_groups) {
    assert_that(is.list(list_of_groups))
    purrr::reduce(lapply(list_of_groups, is_list_of_query_builders), all)
}
