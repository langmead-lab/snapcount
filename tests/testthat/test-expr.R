context("exprs")

test_that("convert simple expression to string", {
    expect_match(exprs(samples_count < 10), c("samples_count < 10"))
})

test_that("convert two simple expressions to strings", {
    expect_setequal(exprs(samples_count > 10, samples_count <= 100), c("samples_count > 10", "samples_count <= 100"))
})

test_that("convert expression contain variable to string", {
    a <- 10
    expect_match(exprs(samples_count >= a), c("samples_count >= 10"))
})

test_that("convert expression containing arithmetic rhs to string", {
    a <- 10
    b <- 11
    expect_match(exprs(samples_count >= a + b), c("samples_count >= 21"))
})

test_that("assert lhs does not get evaluated on conversion", {
    samples_count <- 10
    expect_match(exprs(samples_count >= samples_count), c("samples_count >= 10"))
})

test_that("fail when converting an expression with non-logical operator", {
    expect_error(exprs(samples_count / 10))
})
