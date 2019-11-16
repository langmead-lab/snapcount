context("Test high level functions")

get_full_path_name <- function(file_name) {
    full_path_name <-
        system.file("testdata", file_name, package = "snapcount")
}

test_that("junction union", {
    sb1 <- SnaptronQueryBuilder$new()
    sb1$compilation("encode1159")
    sb1$regions("chr1:1879786-1879786")
    sb1$coordinate_modifier(Coordinates$EndIsExactOrWithin)
    sb1$range_filters(strand == "-")

    sb2 <- sb1$clone(deep = TRUE)
    sb2$compilation("rpc")

    expected <- readRDS(file = get_full_path_name("test_junction_union_output.rds"))
    result <- junction_union(sb1, sb2)

    expect_true(all(expected == result))
})

test_that("shared sample count", {
    ## group 1
    sb1 <- SnaptronQueryBuilder$new()
    sb1$compilation("gtex")
    sb1$regions("chr1:1879786-1879786")
    sb1$coordinate_modifier(Coordinates$EndIsExactOrWithin)
    sb1$range_filters(strand == "-")

    sb2 <- sb1$clone(deep = TRUE)
    sb2$regions("chr1:1879903-1879903")
    sb2$coordinate_modifier(Coordinates$StartIsExactOrWithin)
    group1 <- list(sb1, sb2)

    ## group 2
    sb3 <- sb1$clone(deep = TRUE)
    sb3$regions("chr1:9664595-9664595")
    sb3$range_filters(strand == "+")

    sb4 <- sb2$clone(deep = TRUE)
    sb4$regions("chr1:9664759-9664759")
    sb4$range_filters(strand == "+")
    group2 <- list(sb3, sb4)

    ## group 3
    sb5 <- sb1$clone(deep = TRUE)
    sb5$regions("chr6:32831148-32831148")

    sb6 <- sb2$clone(deep = TRUE)
    sb6$regions("chr6:32831182-32831182")
    group3 <- list(sb5, sb6)

    result <- shared_sample_counts(group1, group2, group3)
    expected <- readRDS(file = get_full_path_name("test_ssc_output.rds"))

    expect_equal(expected, result)

    sb0 <- sb5$clone(deep = TRUE)$coordinate_modifier(Coordinates$Within)
    expect_error(shared_sample_counts(group1, group2, list(sb0, sb6)),
                 "group1 returned no results")

    expect_error(shared_sample_counts(group1, group2, list(sb5, sb0)),
                 "group2 returned no results")
})

test_that("junction inclusion ratio", {
    sb1 <- SnaptronQueryBuilder$new()
    sb1$compilation("srav2")
    sb1$regions("chr2:29446395-30142858")
    sb1$coordinate_modifier(Coordinates$Within)
    sb1$range_filters(strand == "-")

    sb2 <- sb1$clone(deep = TRUE)
    sb2$regions("chr2:29416789-29446394")

    result <- junction_inclusion_ratio(list(sb1), list(sb2))
    expected <- readRDS(file = get_full_path_name("test_jir_output.rds"))

    expect_equal(expected, result)

    sb0 <- sb1$clone(deep = TRUE)$regions("CD99")$coordinate_modifier(Coordinates$Exact)
    expect_error(junction_inclusion_ratio(list(sb0), list(sb2)),
                 "group1 returned no results")

    expect_error(junction_inclusion_ratio(list(sb1), list(sb0)),
                 "group2 returned no results")
})

test_that("percent spliced in", {
    ## inclusion group 1
    sb1 <- SnaptronQueryBuilder$new()
    sb1$compilation("srav2")
    sb1$regions("chr1:94468008-94472172")
    sb1$coordinate_modifier(Coordinates$Exact)
    sb1$range_filters(strand == "+")

    ## inclusion group 2
    sb2 <- sb1$clone(deep = TRUE)
    sb2$regions("chr1:94472243-94475142")

    ## exclusion group
    sb3 <- sb1$clone(deep = TRUE)
    sb3$regions("chr1:94468008-94475142")

    result <- percent_spliced_in(list(sb1), list(sb2), list(sb3), min_count = 1)
    expected <- readRDS(file = get_full_path_name("test_psi_output.rds"))

    expect_equal(expected, result)

    sb0 <- sb1$clone(deep = TRUE)$regions("CD99")
    expect_error(percent_spliced_in(list(sb0), list(sb2), list(sb3)),
                 "inclusion_group1 returned no results")

    expect_error(percent_spliced_in(list(sb1), list(sb0), list(sb3)),
                 "inclusion_group2 returned no results")

    expect_error(percent_spliced_in(list(sb1), list(sb2), list(sb0)),
                 "exclusion_group returned no results")
})

test_that("tissue specificity", {
    sb1 <- SnaptronQueryBuilder$new()
    sb1$compilation("gtex")
    sb1$regions("chr4:20763023-20763023")
    sb1$coordinate_modifier(Coordinates$EndIsExactOrWithin)
    sb1$range_filters(strand == "-")

    sb2 <- sb1$clone(deep = TRUE)
    sb2$regions("chr4:20763098-20763098")
    sb2$coordinate_modifier(Coordinates$StartIsExactOrWithin)

    result <- tissue_specificity(list(sb1, sb2))
    expected <- readRDS(file = get_full_path_name("test_ts_output.rds"))

    expect_equal(expected, result)

    sb0 <- sb1$clone(deep = TRUE)$coordinate_modifier(Coordinates$Within)
    expect_error(tissue_specificity(list(sb0, sb2)),
                 "group1 returned no results")

    expect_error(tissue_specificity(list(sb1, sb0)),
                 "group2 returned no results")
})
