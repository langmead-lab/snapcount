context("Test query functions")

orig.options <- options(test_context = TRUE)

## test_that("regions and compilations are mandatory args", {
##     expect_error(query_jx(compilation = "tcga"), "argument \"regions\" is missing, with no default")
##     expect_error(query_jx(regions = "CD99"), "argument \"compilation\" is missing, with no default")
## })

test_that("simple junction query", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("srav2")
    sb$regions("CD99")

    query_jx(sb)
    expect_equal(uri_of_last_successful_request(), "http://snaptron.cs.jhu.edu/srav2/snaptron?regions=CD99")
})

test_that("using genomic ranges", {
    x1 <- "chr2:100-200:-"
    g_range <- GenomicRanges::GRanges(x1)
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions(g_range)
    sb$column_filters(SMTS == "Brain")

    query_jx(sb)
    expect_equal(uri_of_last_successful_request(),
                 "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=chr2:100-200&rfilter=strand:-&sfilter=SMTS:Brain")
})

test_that("junction query with NSE sample filter", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    sb$column_filters(SMTS == "Brain")

    query_jx(sb)
    expect_equal(uri_of_last_successful_request(),
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&sfilter=SMTS:Brain")
})

test_that("junction query with mutiple NSE sample filters", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("srav2")
    sb$regions("CD99")
    sb$column_filters(library_name == "HG00115.6", study_accession == "ERP001942")

    query_jx(sb)
    expect_equal(uri_of_last_successful_request(),
        "http://snaptron.cs.jhu.edu/srav2/snaptron?regions=CD99&sfilter=library_name:HG00115.6&sfilter=study_accession:ERP001942")
})

test_that("junction query with one NSE range filter", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    sb$row_filters(samples_count >= 5)

    query_jx(sb)
    expect_equal(uri_of_last_successful_request(),
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&rfilter=samples_count>:5")
})

test_that("junction query with mutiple NSE range filters", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    sb$row_filters(samples_count <= 10, coverage_sum < 3)

    query_jx(sb)
    expect_equal(uri_of_last_successful_request(),
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&rfilter=samples_count<:10&rfilter=coverage_sum<3")
})

test_that("invalid sample filter name", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    sb$column_filters(SNTS == "Brain")
    expect_error(
        query_jx(sb),
        "`SNTS' is not a valid sample filter")
})

test_that("invalid sample filter value", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    sb$column_filters(SMTS == 2)
    expect_error(
        query_jx(sb),
        "`SMTS' filter expects value of type String, but got Integer")
})

test_that("junction query with sids", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    sb$sids(1:3)

    query_jx(sb)
    expect_equal(uri_of_last_successful_request(),
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&sids=1,2,3")
})

test_that("query with non-numeric sids", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("tcga")
    sb$regions("CD99")
    sb$sids(c("1", "2", "3"))

    expect_error(query_jx(sb))
})

test_that("test coordinate Coordinates$Exact", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    sb$coordinate_modifier(Coordinates$Exact)

    query_jx(sb)
    expect_equal(uri_of_last_successful_request(),
                 "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&exact=1")
})

test_that("test coordinate Coordinate$Within", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    sb$coordinate_modifier(Coordinates$Within)

    query_jx(sb)
    expect_equal(uri_of_last_successful_request(),
                 "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&contains=1")
})

test_that("test coordinate Coordinates$StartIsExactorWithin", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    sb$coordinate_modifier(Coordinates$StartIsExactOrWithin)

    query_jx(sb)
    expect_equal(uri_of_last_successful_request(),
                 "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&either=1")
})

test_that("test coordinate Coordinates$EndIsExactOrWithin", {
    sb <- SnaptronQueryBuilder$new()
    sb$compilation("gtex")
    sb$regions("CD99")
    sb$coordinate_modifier(Coordinates$EndIsExactOrWithin)

    query_jx(sb)
    expect_equal(uri_of_last_successful_request(),
                 "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&either=2")
})

options(orig.options)
