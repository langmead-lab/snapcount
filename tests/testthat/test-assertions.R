context("Test assertions")

test_that("Empty compilation", {
    sb <- SnaptronQueryBuilder$new()
    expect_error(sb$query_jx(), "Please set a compilation")
})

test_that("Empty region", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    expect_error(sb$query_jx(), "Please specify query regions")
})

test_that("Invalid chromosome range", {
    expect_error(query_jx(compilation = "tcga", regions = "ch0:100-200"))
    expect_error(query_jx(compilation = "tcga", regions = "chr23:100-200"))
    expect_error(query_jx(compilation = "tcga", regions = "chr22X:100-200"))
})

test_that("Invalid sample filters", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    sb <- set_column_filters(sb, SNTS == "Brain")
    expect_error(query_gene(sb),
                 "(is not a valid sample filter)|(Invalid sample filter)")
})

test_that("Invalid compilation", {
    sb <- SnaptronQueryBuilder$new()
    expect_error(sb$compilation("abcde"), "abcde: is not a valid compilation")
})

test_that("Invalid coordinate modifier", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    sb$coordinate_modifier("none")
    expect_error(query_exon(sb), "Invalid coordinate modifier")
})

test_that("Invalid SIDs", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    expect_error(sb$sids(c("1", "2")), "sids should be whole numbers")
    expect_error(sb$sids(c(1.2, 3.4)), "sids should be whole numbers")
})
