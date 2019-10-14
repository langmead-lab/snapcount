context("SnaptronQueryBuilder")

orig.options <- options()

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
    options(test_context = TRUE)
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")$regions("CD99")$sample_filters(SMTS == "Brain")$range_filters(samples_count <= 10)$sids(c(50099,50102))

    sb$query_jx()
    expect_equal(uri_of_last_successful_request(), "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&rfilter=samples_count<:10&sfilter=SMTS:Brain&sids=50099,50102")

    sb$query_gene()
    expect_equal(uri_of_last_successful_request(), "http://snaptron.cs.jhu.edu/gtex/genes?regions=CD99&rfilter=samples_count<:10&sfilter=SMTS:Brain&sids=50099,50102")

    sb$query_exon()
    expect_equal(uri_of_last_successful_request(), "http://snaptron.cs.jhu.edu/gtex/exons?regions=CD99&rfilter=samples_count<:10&sfilter=SMTS:Brain&sids=50099,50102")

    sb$query_coverage()
    expect_equal(uri_of_last_successful_request(), "http://snaptron.cs.jhu.edu/gtex/bases?regions=CD99&sids=50099,50102")
    options(orig.options)
})

test_that("test building from invalid url", {
    sb <- SnaptronQueryBuilder$new()
    expect_error(sb$from_url("http://snap.cs.jhu.edu/gtex/snaptron?regions=CD99&rfilter=samples_count<:10"))
})

test_that("test building from url", {
    sb <- SnaptronQueryBuilder$new()
    sb$from_url("http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&rfilter=samples_count<:10&sfilter=SMTS:Brain&sids=50099,50102")

    expect_equal(sb$compilation(), "gtex")
    expect_equal(sb$regions(), "CD99")
    expect_equal(sb$range_filters(), "samples_count<:10")
    expect_equal(sb$sample_filters(), "SMTS:Brain")
    expect_equal(sb$sids(), c(50099,50102))
})

test_that("test building from url with coordinate modifiers", {
    url <-
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99"
    sb <- SnaptronQueryBuilder$new()

    sb$from_url(paste(url, "either=1", sep = "&"))
    expect_equal(sb$coordinate_modifier(), "StartIsExactOrWithin")
    expect_output(sb$print(), "<SnaptronQueryBuilder>\\n   regions: CD99\\n   coordinate_modifier: either=1\\n   compilation: gtex")

    sb$from_url(paste(url, "either=2", sep = "&"))
    expect_equal(sb$coordinate_modifier(), "EndIsExactOrWithin")
    expect_output(sb$print(), "<SnaptronQueryBuilder>\\n   regions: CD99\\n   coordinate_modifier: either=2\\n   compilation: gtex")

    sb$from_url(paste(url, "exact=1", sep = "&"))
    expect_equal(sb$coordinate_modifier(), "Exact")
    expect_output(sb$print(), "<SnaptronQueryBuilder>\\n   regions: CD99\\n   coordinate_modifier: exact\\n   compilation: gtex")

    sb$from_url(paste(url, "contains=1", sep = "&"))
    expect_equal(sb$coordinate_modifier(), "Within")
    expect_output(sb$print(), "<SnaptronQueryBuilder>\\n   regions: CD99\\n   coordinate_modifier: contains\\n   compilation: gtex")
})

## do this again in case there are any test failures and
## we are not able to restore the original env before
## failure
options(orig.options)
