#' Get the union of junctions from 2 or more compilations
#' which are on the same reference
#'
#' This function queries 2 or more compilations which are on the same
#' reference version (e.g. hg38) and merges the resulting junctions
#' into a single output table, unioning the sample coverage columns
#' and the snaptron_id (jx ID) columns (the latter delimiter will
#' be ":").  All sample IDs will be disjoint between compilations.
#'
#' Union is based on the following fields (combined into a comparison key):
#' \itemize{
#'   \item group
#'   \item chromosome
#'   \item start
#'   \item end
#'   \item strand
#' }
#'
#' The goal is to have a single list of junctions where every junction
#' occurs in at least one compilation and if a junction occurs in > 1
#' compilations it still only has a single row representing all the
#' samples across compilations that it appears in.
#' Sample aggregate statistics will be recalculated for junctions which
#' are merged across *all* samples from all compilations:
#'
#'
#' \itemize{
#'   \item sample_count
#'   \item coverage_sum
#'   \item coverage_avg
#'   \item coverage_median
#' }
#'
#'
#' @param ... One or more SnaptronQueryBuilder objects
#' @return A RangedSummarizedExperiment of junctions appearing in at least
#' one compilation
#' @seealso [junction_intersection()]
#' @export
#' @examples
#' sb1 <- SnaptronQueryBuilder$new()
#' sb1$compilation(Compilation$gtex)
#' sb1$regions("chr1:1879786-1879786")
#' sb1$coordinate_modifier(Coordinates$EndIsExactOrWithin)
#' sb1$range_filters(strand == "-")
#'
#' sb2 <- SnaptronQueryBuilder$new()
#' sb2$compilation("tcga")
#' sb2$regions("chr1:1879786-1879786")
#' sb2$coordinate_modifier(Coordinates$EndIsExactOrWithin)
#' sb2$range_filters(strand == "-")
#'
#' junction_union(sb1, sb2)
junction_union <- function(...) {
    assert_that(is_list_of_query_builders(list(...)),
                msg = "junction_union expects 1 or more SnaptronQueryBuilder objects")
    merge_compilations(..., all = TRUE)
}

#' Get the intersection of junctions from 2 or more compilations
#' which are on the same reference
#'
#' This function operates similar to the `junction_union()` function, i.e
#' it’s cross compilation and merging of the same junction from multiple
#' compilations will be handled exactly the same way. But instead
#' of every junction which appears in at least one compilation, only
#' the junctions which appear in *every* compilation will be returned.
#'
#' @param ... One or more SnaptronQueryBuilder objects
#' @return A RangedSummarizedExperiment of junctions common across compilations
#' @seealso [junction_union()]
#' @export
#' @examples
#' sb1 <- SnaptronQueryBuilder$new()
#' sb1$compilation("gtex")
#' sb1$regions("chr1:1879786-1879786")
#' sb1$coordinate_modifier(Coordinates$EndIsExactOrWithin)
#' sb1$range_filters(strand == "-")
#'
#' sb2 <- SnaptronQueryBuilder$new()
#' sb2$compilation("tcga")
#' sb2$regions("chr1:1879786-1879786")
#' sb2$coordinate_modifier(Coordinates$EndIsExactOrWithin)
#' sb2$range_filters(strand == "-")
#'
#' junction_intersection(sb1, sb2)
junction_intersection <- function(...) {
    assert_that(is_list_of_query_builders(list(...)),
                msg = "junction_intersection expects 1 or more SnaptronQueryBuilder objects")
    merge_compilations(..., all = FALSE)
}

merge_compilations <- function(..., all) {
    metadata_list <- list()
    datasets <- lapply(list(...), function(sb) {
        df <- sb$query_jx(return_rse = FALSE)
        compilation <- sb$compilation()
        compilation_metadata <- get_compilation_metadata(sb$compilation())

        if (is.null(metadata_list[[compilation]])) {
            metadata_list[[compilation]] <<- compilation_metadata
        }
        return(df)
    })

    if (length(datasets) == 1) {
        return(rse(datasets[[1]], metadata_list[[1]]))
    }

    for (i in seq_along(datasets)) {
        if (i == 1) {
            res <- datasets[[1]]
        } else {
            res <- merge(res, datasets[[i]], all = all,
                         by = c("chromosome", "start", "end", "strand")) %>%
                finalize_merge(col_names = names(datasets[[1]]))
        }
    }

    for (i in seq_along(metadata_list)) {
        if (i == 1) {
            metadata <- metadata_list[[1]]
        } else {
            metadata <- merge(metadata, metadata_list[[i]], by = "rail_id", all = TRUE)
        }
    }

    rse(res, metadata)
}

finalize_merge <- function(dt, col_names) {
    dt$samples <- str_cat(dt$samples.x, dt$samples.y)
    dt$snaptron_id <- str_cat(dt$snaptron_id.x, dt$snaptron_id.y, sep = ",")

    dt$`DataSource:Type` <- str_cat(dt$`DataSource:Type.x`,
                                    dt$`DataSource:Type.y`, sep = ",")
    dt$source_dataset_id <- str_cat(dt$source_dataset_id.x,
                                    dt$source_dataset_id.y, sep = ",")

    dt$coverage_sum <- purrr::map2_int(dt$coverage_sum.x,
                                   dt$coverage_sum.y, sum, na.rm = TRUE)
    dt$samples_count <- purrr::map2_int(dt$samples_count.x,
                                    dt$samples_count.y, sum, na.rm = TRUE)
    dt$coverage_avg <- dt$coverage_sum / dt$samples_count
    dt$coverage_median <-
        stringr::str_extract_all(dt$samples, "\\b\\d+\\b") %>%
        lapply(calculate_coverage_median) %>% unlist()

    dt$length <- choose_non_na(dt$length.x, dt$length.y)
    dt$left_motif <- choose_non_na(dt$left_motif.x, dt$left_motif.y)
    dt$right_motif <- choose_non_na(dt$right_motif.x, dt$right_motif.y)
    dt$annotated <- choose_non_na(dt$annotated.x, dt$annotated.y)
    dt$left_annotated <- choose_non_na(dt$left_annotated.x, dt$left_annotated.y)
    dt$right_annotated <- choose_non_na(dt$right_annotated.x, dt$right_annotated.y)

    dt[, col_names, with = FALSE]
}

str_cat <- function(..., sep = "") {
    strings <- lapply(list(...), stringr::str_replace_na, replacement = "")
    paste(strings[[1]], strings[[2]], sep = sep)
}

choose_non_na <- function(x, y) {
    ifelse(is.na(x), y, x)
}

calculate_coverage_median <- function(samples) {
    samples[c(FALSE, TRUE)] %>% as.numeric() %>% median()
}

#' Relative measure of splice variant usage similar to PSI that allows
#' for 2 arbitrarily defined groups of junctions (not limited to
#' cassette exons).
#'
#' Calculates a coverage summary statistic per sample of the normalized
#' coverage difference between two sets of separate junctions defined
#' by at least two basic queries and organized into two groups.
#'
#' The summary statistic is as follows:
#' If the coverage of the first group is "A" and the second is "B":
#'
#' `JIR(A,B)=(A - B) / (A+B+1)`
#'
#' This is calculated for every sample that occurs in one or the other
#' (or both) groups’ results.
#'
#' @param group1,group2 Each group is a list of 1 or more SnaptronQueryBuilder objects
#' @param group_names Optional vector of strings representing the group names
#' @return A DataFrame of samples, with their JIR score and metadata, which had > 0 coverage
#'  in at least one resulting row in at least one of the groups
#'
#' @examples
#' sb1 <- SnaptronQueryBuilder$new()
#' sb1$compilation("srav2")
#' sb1$regions("chr2:29446395-30142858")
#' sb1$coordinate_modifier(Coordinates$Within)
#' sb1$range_filters(strand == "-")
#'
#' sb2 <- SnaptronQueryBuilder$new()
#' sb2$compilation("srav2")
#' sb2$regions("chr2:29416789-29446394")
#' sb2$coordinate_modifier(Coordinates$Within)
#' sb2$range_filters(strand == "-")
#'
#' junction_inclusion_ratio(list(sb1), list(sb2))
#' @export
junction_inclusion_ratio <- function(group1, group2, group_names = NULL) {
    assert_that(is_list_of_query_builders(group1),
                is_list_of_query_builders(group2))

    query_results <- run_queries(group1, group2)
    if (is.null(query_results[[1]])) {
        stop("Unable to calculate JIR: group1 returned no results")
    }
    if (is.null(query_results[[2]])) {
        stop("Unable to calculate JIR: group2 returned no results")
    }

    jir <- merge(query_results[[1]], query_results[[2]],
                 by = "sample_id", all = TRUE) %>% replace_na(0)

    ## jir[, jir := calc_jir(coverage.y, coverage.x)]
    jir$jir <- calc_jir(jir$coverage.y, jir$coverage.x)

    if (is.null(group_names)) {
        group_names <- c("group1_coverage", "group2_coverage")
    } else {
        group_names <- paste(group_names, "coverage", sep = "_")
    }

    data.table::setnames(jir, old = "coverage.x", new = group_names[[1]])
    data.table::setnames(jir, old = "coverage.y", new = group_names[[2]])

    jir[order(jir, decreasing = TRUE)]
}

calc_jir <- function(a, b) {
    (a - b)/(a + b + 1)
}

#' Relative measure of splice variant usage, limited currently to
#' cassette exon splice variants
#'
#' Similar to the JIR, this calculates Percent Spliced In (PSI) statistics
#' for the definition of 2 different groups: inclusion and exclusion.
#' Currently this function only supports the cassette exon use case.
#'
#' Inclusion typically defines 2 basic queries, one for the junction
#' preceding the cassette exon, and the second for the junction following
#' the cassette exon.  The exclusion group contains one basic query
#' which defines the junction which skips the cassette exon.
#'
#' The PSI itself is implemented as:
#'
#'
#' \code{PSI(inclusion1, inclusion2, exclusion) =
#'   mean(inclusion1, inclusion2) / (mean(inclusion1, inclusion2) + exclusion)}
#'
#' where each term denotes the coverage of junctions that resulted
#' from the basic queries in that group in the current sample.
#'
#' @param inclusion_group1,inclusion_group2,exclusion_group Where each is a list of 1 or
#'   more SnaptronQueryBuilder objects
#' @param min_count minimum total count (denominator) required to not be assigned -1
#' @param group_names Optional vector of strings representing the group names
#' @return A DataFrame of samples, with their PSI score and metadata, which had > 0 coverage
#'  in at least one resulting row in at least one of the groups
#'
#' @examples
#' inclusion_group1 <- SnaptronQueryBuilder$new()
#' inclusion_group1$compilation("srav2")
#' inclusion_group1$regions("chr1:94468008-94472172")
#' inclusion_group1$coordinate_modifier(Coordinates$Exact)
#' inclusion_group1$range_filters(strand == "+")
#'
#' inclusion_group2 <- SnaptronQueryBuilder$new()
#' inclusion_group2$compilation("srav2")
#' inclusion_group2$regions("chr1:94468008-94472172")
#' inclusion_group2$coordinate_modifier(Coordinates$Exact)
#' inclusion_group2$range_filters(strand == "+")
#'
#' exclusion_group <- SnaptronQueryBuilder$new()
#' exclusion_group$compilation("srav2")
#' exclusion_group$regions("chr1:94468008-94475142")
#' exclusion_group$coordinate_modifier(Coordinates$Exact)
#' exclusion_group$range_filters(strand == "+")
#'
#' percent_spliced_in(list(inclusion_group1), list(inclusion_group2), list(exclusion_group))
#' @export
percent_spliced_in <- function(inclusion_group1, inclusion_group2,
                               exclusion_group, min_count = 20,
                               group_names = NULL)
{
    query_results <-
        run_queries(inclusion_group1, inclusion_group2, exclusion_group)
    if (is.null(query_results[[1]])) {
        stop("Unable to calculate PSI: inclusion_group1 returned no results")
    }
    if (is.null(query_results[[2]])) {
        stop("Unable to calculate PSI: inclusion_group2 returned no results")
    }
    if (is.null(query_results[[3]])) {
        stop("Unable to calculate PSI: exclusion_group returned no results")
    }

    psi <-
        merge(query_results[[1]], query_results[[2]],
                 by = "sample_id", all = TRUE) %>%
        merge(query_results[[3]], by = "sample_id", all = TRUE) %>%
        replace_na(0)

    ## psi[, psi := calc_psi(coverage.x, coverage.y, coverage, min_count)]
    psi$psi <- calc_psi(psi$coverage.x, psi$coverage.y, psi$coverage, min_count)

    if (is.null(group_names)) {
        group_names <- c("inclusion_group1_coverage",
                         "inclusion_group2_coverage",
                         "exclusion_group_coverage")
    } else {
        group_names <- paste(group_names, "coverage", sep = "_")
    }

    data.table::setnames(psi, old = "coverage.x", new = group_names[[1]])
    data.table::setnames(psi, old = "coverage.y", new = group_names[[2]])
    data.table::setnames(psi, old = "coverage", new = group_names[[3]])

    psi
}

calc_psi <- function(inclusion1, inclusion2, exclusion, min_count) {
    mean_inclusion <- (inclusion1 + inclusion2) / 2.
    total <- mean_inclusion + exclusion
    psi <- mean_inclusion / total

    ifelse(inclusion1 == 0 | inclusion2 == 0 | total < min_count, -1, psi)
}

#' Tissue Specificity (TS): produces a list of samples with their tissues
#' marked which either contain queried junctions (1) or not (0); can be used
#' as input to significance testing methods such as Kruskal-Wallis to look for
#' tissue enrichment (currently only works for the GTEx compilation).
#'
#' Lists the number of samples labeled with a specific tissue type.

#' Samples are filtered for ones which have junctions across all the
#' user-specified groups. That is, if a sample only appears in the results of
#' some of the groups (from their basic queries) it will be assigned a 0,
#' otherwise if it’s in all of the groups’ results it will be assigned a 1.
#' This is similar to the SSC high level query type, but doesn’t sum the
#' coverage.
#'
#' The samples are then grouped by their tissue type (e.g. Brain).
#' This is useful for determining if there’s an enrichment for a specific
#' tissue in the set of junctions queried.  Results from this can be fed to a
#' statistical test, such as the Kruskal-wallis non-parametric rank test.

#' This query is limited to GTEx only due to the fact that GTEx is one of the
#' few compilations that has consistent and complete tissue metadata.
#'
#' @param ... One or more SnaptronQueryBuilder objects
#' @param group_names Optional vector of strings representing the group names
#' @return A DataFrame of all samples in the compilation with either a 0 or 1
#'  indicating their occurrence and shared status (if > 1 group passed in).
#'  Occurrence here is if the sample has at least one result with > 0 coverage,
#'  and further, if > 1 group passed in, if it occurs in the results of all groups.
#'  Also includes the sample tissue type and sample_id.
#'
#' @examples
#' inclusion_group1 <- SnaptronQueryBuilder$new()
#' inclusion_group1$compilation("gtex")
#' inclusion_group1$regions("chr4:20763023-20763023")
#' inclusion_group1$coordinate_modifier(Coordinates$EndIsExactOrWithin)
#' inclusion_group1$range_filters(strand == "-")
#'
#' inclusion_group2 <- SnaptronQueryBuilder$new()
#' inclusion_group2$compilation("gtex")
#' inclusion_group2$regions("chr4:20763098-20763098")
#' inclusion_group2$coordinate_modifier(Coordinates$StartIsExactOrWithin)
#' inclusion_group2$range_filters(strand == "-")
#'
#' tissue_specificity(list(inclusion_group1, inclusion_group2))
#' @export
tissue_specificity <- function(..., group_names = NULL) {
    list_of_groups <- list(...)
    assert_that(is_list_of_query_builder_groups(list_of_groups),
                msg = paste("tissue_specificity expects 1 or",
                            "more list of SnaptronQueryBuilder objects"))
    num_groups <- length(list_of_groups)

    if (is.null(group_names)) {
        group_names <- paste0("g", seq_along(num_groups))
    }

    dfs <- lapply(seq_along(num_groups), function(i) {
        g <- list_of_groups[[i]]
        name <- group_names[[i]]
        if (length(g) == 1) {
            tissue_specificity_per_group(g[[1]], g[[1]], name)
        } else{
            tissue_specificity_per_group(g[[1]], g[[2]], name)
        }
    })

    do.call(rbind, dfs)
}

tissue_specificity_per_group <- function(group1, group2, group_name) {
    if (is_query_builder(group1)) {
        group1 <- list(group1)
    }

    if (is_query_builder(group2)) {
        group2 <- list(group2)
    }

    stopifnot(is.list(group1), is.list(group2))
    query_data <- run_queries(group1, group2, summarize = FALSE)
    if (is.null(query_data[[1]])) {
        stop("Unable to calculate TS: group1 returned no results")
    }
    if (is.null(query_data[[2]])) {
        stop("Unable to calculate TS: group2 returned no results")
    }

    merged_data <-
        merge(query_data[[1]], query_data[[2]],
              by = "sample_id", all = TRUE) %>% replace_na(0)
    merged_data$shared <- shared(merged_data$coverage.x, merged_data$coverage.y)
    merged_data <- merged_data[, !c("coverage.x", "coverage.y")]

    metadata <- get_compilation_metadata(group1[[1]]$compilation())
    metadata <- metadata[, c("rail_id", "SMTS")]

    ts <- merge(merged_data, metadata, by.x = "sample_id",
                by.y = "rail_id", all = TRUE)
    ts <- replace_na(ts, 0, colnames = setdiff(names(ts), "SMTS"))

    ts$group <- rep(group_name, nrow(ts))
    data.table::setnames(ts, old = "SMTS", new = "tissue")

    unique(ts)
}

shared <- function(cov1, cov2) {
    cov1 != 0 & cov2 != 0
}

#' Shared Sample Count (SSC): counts total number of samples in which 2
#' different junctions both occur in.
#'
#' This produces a list of user-specified groups and the read coverage of the
#' junctions in all the samples which were shared across all the basic queries
#' occurring in each group.
#'
#' Example: User defines a single group of junctions "GroupA" made up of 2
#' separate regions (two basic queries).
#'
#' An SSC query will return a single line for GroupA which will have the total
#' number of samples which had at least one junction which was returned from both basic
#' queries. It will also report a summary statistic of the total number of
#' groups which had one or more samples that were shared across the basic
#' queries, in this case it would be 1.  Also, it will report the number of
#' groups which had at least one shared sample and which had matching
#' junctions (from the query) which were fully annotated.

#' This function can be used to determine how much cross-sample support there
#' is for a particular junction configuration (typically a cassette exon).
#'
#' @param ... One or more lists of SnaptronQueryBuilder objects
#' @param group_names Optional vector of character strings representing group names
#' @return A DataFrame of results based on the list of groups passed in via "group_names".
#'  Each group is reported with the # of unique samples which occurred in all of its defined set
#'  of related basic queries (e.g. two inclusion basic queries in a cassette exon scenario).
#'
#' @examples
#' group1 <- SnaptronQueryBuilder$new()
#' group1$compilation("gtex")
#' group1$regions("chr1:1879786-1879786")
#' group1$coordinate_modifier(Coordinates$EndIsExactOrWithin)
#' group1$range_filters(strand == "-")
#'
#' group2 <- SnaptronQueryBuilder$new()
#' group2$compilation("gtex")
#' group2$regions("chr1:1879903-1879903")
#' group2$coordinate_modifier(Coordinates$StartIsExactOrWithin)
#' group2$range_filters(strand == "-")
#'
#' ssc<-shared_sample_counts(list(group1, group2))
#' @export
shared_sample_counts <- function(..., group_names = NULL) {
    list_of_groups <- list(...)
    assert_that(is_list_of_query_builder_groups(list_of_groups),
                msg = paste("shared_sample_counts expects 1 or",
                            "more lists of SnaptronQueryBuilder objects"))

    counts <- lapply(list_of_groups, function(g) {
        shared_sample_count(g[[1]], g[[2]])
    }) %>% unlist()

    if (is.null(group_names)) {
        group_names <- paste0("g", seq_along(list_of_groups))
    }

    data.table(group = group_names, counts = counts)
}

shared_sample_count <- function(group1, group2) {
    if (is_query_builder(group1)) {
        group1 <- list(group1)
    }

    if (is_query_builder(group2)) {
        group2 <- list(group2)
    }

    stopifnot(is.list(group1), is.list(group2))

    query_results <- run_queries(group1, group2)
    if (is.null(query_results[[1]])) {
        stop("Unable to calculate SSC: group1 returned no results")
    }
    if (is.null(query_results[[2]])) {
        stop("Unable to calculate SSC: group2 returned no results")
    }
    intersect(query_results[[1]]$sample_id, query_results[[2]]$sample_id) %>%
        length()
}

run_queries <- function(..., summarize = TRUE) {
    lapply(list(...), count_samples, summarize = summarize)
}

count_samples <- function(group, summarize = TRUE) {
    dfs <- lapply(group, function(e) {
        q <- e$query_jx(return_rse = FALSE)
        if (is.null(q)) {
            return()
        }

        data.table::data.table(sample = extract_samples(q)) %>%
            tidyr::separate("sample", into = c("sample_id", "coverage"), convert = TRUE)
    })

    res <- do.call(rbind, dfs)

    if (summarize && !is.null(res)) {
        ## res[, .(coverage = sum(coverage)), by = .(sample_id)]
        res[, lapply(.SD, sum), by = c("sample_id"), .SDcols = c("coverage")]
    } else {
        res
    }
}

replace_na <- function(dt, replacement, colnames = NULL) {
    if (!is.null(colnames)) {
        for (name in colnames) {
            set(dt, which(is.na(dt[[name]])), name, replacement)
        }
    } else {
        for (i in seq_along(dt)) {
            set(dt, which(is.na(dt[[i]])), i, replacement)
        }
    }

    dt
}
