context("Test assertions")

test_that("Invalid chromosome range", {
    expect_error(query_jx(compilation = "tcga", regions = "ch0:100-200"))
    expect_error(query_jx(compilation = "tcga", regions = "chr23:100-200"))
    expect_error(query_jx(compilation = "tcga", regions = "chr22X:100-200"))
})

test_that("Invalid sample filters", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    sb$sample_filters(SNTS == "Brain")
    expect_error(query_gene(sb),
                 "(is not a valid sample filter)|(Invalid sample filter)")
})

test_that("Invalid compilation", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("abcde")
    sb$regions("CD99")
    expect_error(query_gene(sb), "Invalid compilation")
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
    sb$sids(c("1", "2"))
    expect_error(query_exon(sb), "sids should be whole numbers")
    sb$sids(c(1.2, 3.4))
    expect_error(query_exon(sb), "sids should be whole numbers")
})
