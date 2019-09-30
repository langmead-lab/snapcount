context("SnaptronQueryBuilder")

test_that("empty query builder object", {
    sb <- SnaptronQueryBuilder$new()
    expect_output(sb$print(), "<SnaptronQueryBuilder>")
})

test_that("print method with set attributes", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("tcga")
    sb$regions("CD99")
    sb$sids(1:3)

    expect_output(sb$print(), "<SnaptronQueryBuilder>\\n   compilation: tcga\\n   regions: CD99\\n   sids: 1,2,3")
})

test_that("get/set methods", {
    sb <- SnaptronQueryBuilder$new()

    sb$compilation("gtex")
    expect_equal(sb$compilation(), "gtex")

    sb$regions("CD99")
    expect_equal(sb$regions(), "CD99")

    sb$range_filters(samples_count <= 10)
    expect_equal(as.character(sb$range_filters()), "samples_count <= 10")

    sb$sample_filters(SMTS == "Brain")
    expect_equal(as.character(sb$sample_filters()), "SMTS == Brain")

    sb$sids(1:10)
    expect_equal(sb$sids(), 1:10)
})

test_that("query methods", {
    options(test_ctx = TRUE)
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")$regions("CD99")$sample_filters(SMTS == "Brain")$range_filters(samples_count <= 10)$sids(1:5)

    sb$query_jx()
    expect_equal(uri_of_last_successful_request(), "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&rfilter=samples_count<:10&sfilter=SMTS:Brain&sids=1,2,3,4,5")

    sb$query_gene()
    expect_equal(uri_of_last_successful_request(), "http://snaptron.cs.jhu.edu/gtex/genes?regions=CD99&rfilter=samples_count<:10&sfilter=SMTS:Brain&sids=1,2,3,4,5")

    sb$query_exon()
    expect_equal(uri_of_last_successful_request(), "http://snaptron.cs.jhu.edu/gtex/exons?regions=CD99&rfilter=samples_count<:10&sfilter=SMTS:Brain&sids=1,2,3,4,5")

    sb$query_coverage()
    expect_equal(uri_of_last_successful_request(), "http://snaptron.cs.jhu.edu/gtex/bases?regions=CD99&sids=1,2,3,4,5")

    options(test_ctx = NULL)
})

test_that("test building from invalid url", {
    sb <- SnaptronQueryBuilder$new()
    expect_error(sb$from_url("http://snap.cs.jhu.edu/gtex/snaptron?regions=CD99&rfilter=samples_count<:10"))
})

test_that("test building from url", {
    options(test_ctx = TRUE)

    sb <- SnaptronQueryBuilder$new()
    sb$from_url("http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&rfilter=samples_count<:10&sfilter=SMTS:Brain&sids=1,2,3,4,5")

    expect_equal(sb$compilation(), "gtex")
    expect_equal(sb$regions(), "CD99")
    expect_equal(sb$range_filters(), "samples_count<:10")
    expect_equal(sb$sample_filters(), "SMTS:Brain")
    expect_equal(sb$sids(), 1:5)

    options(test_ctx = NULL)
})
