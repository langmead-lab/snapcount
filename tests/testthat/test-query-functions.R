## context("Test query functions")

## orig.options <- options(test_context = TRUE)

## test_that("regions and compilations are mandatory args", {
##     expect_error(query_jx(compilation = "tcga"), "argument \"regions\" is missing, with no default")
##     expect_error(query_jx(regions = "CD99"), "argument \"compilation\" is missing, with no default")
## })

## test_that("simple junction query", {
##     query_jx(compilation = "srav2", regions = "CD99")
##     expect_equal(uri_of_last_successful_request(), "http://snaptron.cs.jhu.edu/srav2/snaptron?regions=CD99")
## })

## test_that("using genomic ranges", {
##     library(GenomicRanges)
##     x1 <- "chr2:100-200:-"
##     g_range <- as(x1, "GRanges")
##     query_jx(compilation = "gtex", regions = g_range, sample_filters = SMTS == "Brain")
##     expect_equal(uri_of_last_successful_request(),
##                  "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=chr2:100-200&rfilter=strand:-&sfilter=SMTS:Brain")
## })

## test_that("junction query with NSE sample filter", {
##     query_jx(compilation = "gtex", regions = "CD99", sample_filters = SMTS == "Brain")
##     expect_equal(uri_of_last_successful_request(),
##         "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&sfilter=SMTS:Brain")
## })

## test_that("junction query with mutiple NSE sample filters", {
##     query_jx(compilation = "srav2", regions = "CD99", sample_filters = c(library_name == "HG00115.6", study_accession == "ERP001942"))
##     expect_equal(uri_of_last_successful_request(),
##         "http://snaptron.cs.jhu.edu/srav2/snaptron?regions=CD99&sfilter=library_name:HG00115.6&sfilter=study_accession:ERP001942")
## })

## test_that("junction query with one NSE range filter", {
##     query_jx(compilation = "gtex", regions = "CD99", range_filters = samples_count >= 5)
##     expect_equal(uri_of_last_successful_request(),
##         "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&rfilter=samples_count>:5")
## })

## test_that("junction query with mutiple NSE range filters", {
##     query_jx(compilation = "gtex", regions = "CD99", range_filters = list(samples_count <= 10, coverage_sum < 3))
##     expect_equal(uri_of_last_successful_request(),
##         "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&rfilter=samples_count<:10&rfilter=coverage_sum<3")
## })

## test_that("invalid sample filter name", {
##     expect_error(
##         query_jx(compilation = "gtex", regions = "CD99", sample_filters = SNTS == "Brain"),
##         "`SNTS' is not a valid sample filter")
## })

## test_that("invalid sample filter value", {
##     expect_error(
##         query_jx(compilation = "gtex", regions = "CD99", sample_filters = SMTS == 2),
##         "`SMTS' filter expects value of type String, but got Integer")
## })

## test_that("NSE rhs does not evaluate to basic type", {
##     expect_error(
##         query_jx(compilation = "gtex", regions = "CD99", sample_filters = SMTS == c(1,2,3)),
##         "does not evaluate to a basic type")
## })

## test_that("junction query with sids", {
##     query_jx(compilation = "gtex", regions = "CD99", sids = 1:3)
##     expect_equal(uri_of_last_successful_request(),
##         "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&sids=1,2,3")
## })

## test_that("query with non-numeric sids", {
##     expect_error(query_jx(compilation = "tcga", regions = "CD99", sids = c("1", "2", "3")))
## })

## test_that("test coordinate Coordinates$Exact", {
##     query_jx(compilation = "gtex", regions = "CD99", coordinate_modifier = Coordinates$Exact)
##     expect_equal(uri_of_last_successful_request(),
##                  "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&exact=1")
## })

## test_that("test coordinate Coordinate$Within", {
##     query_jx(compilation = "gtex", regions = "CD99", coordinate_modifier = Coordinates$Within)
##     expect_equal(uri_of_last_successful_request(),
##                  "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&contains=1")
## })

## test_that("test coordinate Coordinates$StartIsExactorWithin", {
##     query_jx(compilation = "gtex", regions = "CD99", coordinate_modifier = Coordinates$StartIsExactOrWithin)
##     expect_equal(uri_of_last_successful_request(),
##                  "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&either=1")
## })

## test_that("test coordinate Coordinates$EndIsExactOrWithin", {
##     query_jx(compilation = "gtex", regions = "CD99", coordinate_modifier = Coordinates$EndIsExactOrWithin)
##     expect_equal(uri_of_last_successful_request(),
##                  "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=CD99&either=2")
## })

## options(orig.options)
