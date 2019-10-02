context("Test assertions")

test_that("Invalid chromosome range", {
    expect_error(query_jx(compilation = "tcga", regions = "ch0:100-200"))
    expect_error(query_jx(compilation = "tcga", regions = "chr23:100-200"))
    expect_error(query_jx(compilation = "tcga", regions = "chr22X:100-200"))
})

test_that("Invalid sample filters", {
    expect_error(query_gene(compilation = "gtex", regions = "CD99",
                            sample_filters = SNTS == "Brain"),
                 "(is not a valid sample filter)|(Invalid sample filter)")
})

test_that("Invalid compilation", {
    expect_error(query_gene(compilation = "abcde", regions = "CD99"),
                 "Invalid compilation")
})

test_that("Invalid coordinate modifier", {
    expect_error(query_exon(compilation = "gtex", regions = "CD99",
                            coordinate_modifier = "none"),
                 "Invalid coordinate modifier")
})

test_that("Invalid SIDs", {
    expect_error(query_exon(compilation = "gtex", regions = "CD99",
                            sids = c("1", "2")),
                 "sids should be whole numbers")
    expect_error(query_exon(compilation = "gtex", regions = "CD99",
                            sids = c(1.2, 3.4)),
                 "sids should be whole numbers")
})
