context("QueryBuilderWrappers")

global_qb <- QueryBuilder(compilation = "gtex", regions = "CD99")

test_that("test initializer", {
    qb <- QueryBuilder(compilation = "gtex", regions = "CD99")

    expect_equal(qb$compilation(), "gtex")
    expect_equal(qb$regions(), "CD99")
})

test_that("test get/set compilation", {
    qb <- set_compilation(global_qb, "tcga")

    expect_equal(get_compilation(qb), "tcga")
})

## test_that("test get/set regions", {
##     qb <- set_regions(global_qb, "chr4:1-1000")

##     expect_equal(get_regions(qb), "ch4:1-1000")
## })

test_that("test get/set column filters as strings", {
    qb <- set_column_filters(global_qb, "foo == bar", "a < b")

    filters <- get_column_filters(qb)
    expect_length(filters, 2)
    expect_equal(filters[[1]], "foo == bar")
    expect_equal(filters[[2]], "a < b")
})

test_that("test get/set column filters as expressions", {
    qb <- set_column_filters(global_qb, foo == "bar", a < (1 + 2))

    filters <- get_column_filters(qb)
    expect_length(filters, 2)
    expect_equal(filters[[1]], "foo == bar")
    expect_equal(filters[[2]], "a < 3")

})

test_that("test get/set row filters as strings", {
    qb <- set_row_filters(global_qb, "foo == bar", "a < b")

    filters <- get_row_filters(qb)
    expect_length(filters, 2)
    expect_equal(filters[[1]], "foo == bar")
    expect_equal(filters[[2]], "a < b")
})

test_that("test get/set row filters as expressions", {
    qb <- set_row_filters(global_qb, foo == "bar", a < (1 + 2))

    filters <- get_row_filters(qb)
    expect_length(filters, 2)
    expect_equal(filters[[1]], "foo == bar")
    expect_equal(filters[[2]], "a < 3")

})

test_that("get get/set sids", {
    qb <- set_sids(global_qb, 1:10)

    expect_equal(get_sids(qb), 1:10)
})

test_that("get/set coordinate_modifiers", {
    qb <- set_coordinate_modifier(global_qb, Coordinates$Exact)

    expect_equal(get_coordinate_modifier(qb), "Exact")
})

test_that("query builder from url", {
    qb <- from_url("http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99")

    expect_equal(get_compilation(qb), "gtex")
    expect_equal(get_regions(qb), "CD99")
})
