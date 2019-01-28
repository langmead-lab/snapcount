context("Test query functions")

options(test_ctx = TRUE)

test_that("genes_or_intervals and compilations are mandatory args", {
    expect_error(query_jx(compilation = "tcga"), "please specify either a gene or an interval")
    expect_error(query_jx(genes_or_intervals = "CD99"), "compilation is a required argument")
})

test_that("simple junction query", {
    query_jx(compilation = "srav2", genes_or_intervals = "CD99")
    expect_equal(uri_of_last_successful_request(), "http://snaptron.cs.jhu.edu/srav2/snaptron?regions=CD99")
})

test_that("junction query with one sample filter", {
    query_jx(compilation = "gtex", genes_or_intervals = "CD99", sample_filters = exprs(SMTS == "Brain"))
    expect_equal(uri_of_last_successful_request(),
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&sfilter=SMTS:Brain")
})

test_that("junction query with mutiple sample filters", {
    query_jx(compilation = "gtex", genes_or_intervals = "CD99", sample_filters = exprs(description == "Cortex", library_strategy == "RNA-Seq"))
    expect_equal(uri_of_last_successful_request(),
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&sfilter=description:Cortex&sfilter=library_strategy:RNA-Seq")
})

test_that("junction query with one range filter", {
    query_jx(compilation = "gtex", genes_or_intervals = "CD99", range_filters = exprs(samples_count >= 5))
    expect_equal(uri_of_last_successful_request(),
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&rfilter=samples_count>:5")
})

test_that("junction query with mutiple range filters", {
    query_jx(compilation = "gtex", genes_or_intervals = "CD99", range_filters = exprs(samples_count <= 10, coverage_sum < 3))
    expect_equal(uri_of_last_successful_request(),
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&rfilter=samples_count<:10&rfilter=coverage_sum<3")
})

test_that("junction query with sids", {
    query_jx(compilation = "gtex", genes_or_intervals = "CD99", sids = 1:3)
    expect_equal(uri_of_last_successful_request(),
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&sids=1,2,3")
})

test_that("sids must be integers", {
    expect_error(query_jx(compilation = "tcga", genes_or_intervals = "CD99", sids = c("1", "2", "3")))
})
