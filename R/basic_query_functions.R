#' A Reference Class for building Snaptron queries
#'
#' @section Methods:
#' \code{compilation} Get or set query compilation. See
#'   \code{\link{Compilation}} for details.
#'
#' \code{regions} Get or set query region which can be either a chromosome range
#'   or HUGO gene name.
#'
#' \code{range_filters} Get or set query range filters. Range filters can be set
#'   to a list or vector of either strings or R boolean expressions. If filters
#'   are R boolean expressions then the right-hand-side of the R expression
#'   should evaluate to either a string, float or whole number.  The expression
#'   will then be converted to a string before being sent out to Snaptron.
#'
#' \code{sample_filters} Get or set query sample filters. When setting sample
#'   filters the argument is expected to follow the sample rules as described
#'   for \code{range_filters}.
#'
#' \code{sids} Get or set query sample ids.
#'
#' \code{coordinate_modifier} Get or set coordinate modifiers for query. See
#'   \code{\link{Coordinates}} for details.
#'
#' \code{query_jx} call \code{\link{query_jx}} function
#'
#' \code{query_gene} call \code{\link{query_gene}} function
#'
#' \code{query_exon} call \code{\link{query_exon}} function
#'
#' \code{print} print query builder object
#'
#' \code{from_url} use URL to instantiate SnaptronQueryBuilder object
#'
#' @export
#' @examples
#' sb <- SnaptronQueryBuilder$new()
#' sb$compilation("gtex")$regions("CD99")$query_jx()
#' sb$from_url("http://snaptron.cs.jhu.edu/gtex/snaptron?regions=chr1:1-100000")
#' sb$sample_filters(SMTS == "Brain")
#'
#' # The RHS of an R boolean expression can itself be an expression,
#' # but the expression has to evaluate to either a String, Float or Integer.
#' a <- 10
#' b <- 10
#' sb$range_filters(samples_count >= (a + b))
#' sb$query_jx(return_rse = FALSE)
SnaptronQueryBuilder <-
    R6Class("SnaptronQueryBuilder",
            public = list(
                initialize = function(...) {
                    private$query <- list(...)
                },
                compilation = function(compilation = NULL) {
                    if (!missing(compilation)) {
                        private$query$compilation <- compilation
                        invisible(self)
                    } else {
                        private$query$compilation
                    }
                },
                regions = function(regions = NULL) {
                    if (!missing(regions)) {
                        private$query$regions <- regions
                        invisible(self)
                    } else {
                        private$query$regions
                    }
                },
                range_filters = function(range_filters = NULL) {
                    if (!missing(range_filters)) {
                        private$query$range_filters <-
                            bool_expressions_to_strings(
                                rlang::enexpr(range_filters)
                            )
                        invisible(self)
                    } else {
                        private$query$range_filters
                    }
                },
                sample_filters = function(sample_filters = NULL) {
                    if (!missing(sample_filters)) {
                        private$query$sample_filters <-
                            bool_expressions_to_strings(
                                rlang::enexpr(sample_filters)
                            )
                        invisible(self)
                    } else {
                        private$query$sample_filters
                    }
                },
                sids = function(sids = NULL) {
                    if (!missing(sids)) {
                        private$query$sids <- sids
                        invisible(self)
                    } else {
                        private$query$sids
                    }
                },
                coordinate_modifier = function(coordinate_modifier = NULL) {
                    if (!missing(coordinate_modifier)) {
                        private$query$coordinate_modifier <- coordinate_modifier
                        invisible(self)
                    } else {
                        private$query$coordinate_modifier
                    }
                },
                query_jx = function(return_rse = TRUE) {
                    private$call(
                                "query_jx",
                                c(list(return_rse = return_rse), private$query)
                            )
                },
                query_exon = function(return_rse = TRUE) {
                    private$call(
                                "query_exon",
                                c(list(return_rse = return_rse), private$query)
                            )
                },
                query_gene = function(return_rse = TRUE) {
                    private$call(
                                "query_gene",
                                c(list(return_rse = return_rse), private$query)
                            )
                },
                from_url = function(url) {
                    url <- httr::parse_url(url)
                    if (url$hostname != "snaptron.cs.jhu.edu") {
                        stop(
                            "URL does not point to Snaptron server",
                            call. = FALSE
                        )
                    }
                    resp <- httr::HEAD(url)
                    if (resp$status_code != 200 ||
                        httr::http_type(resp) != "text/plain") {
                        stop(
                            sprintf("%s: is not a valid URL", url),
                            call. = FALSE
                        )
                    }
                    query <- list()
                    for (i in seq_along(url$query)) {
                        name <- switch(
                            n <- names(url$query[i]),
                            rfilter = "range_filters",
                            sfilter = "sample_filters",
                            regions = "regions",
                            n
                        )

                        if (name == "sids") {
                            query[[name]] <-
                                scan(textConnection(url$query[[i]]), sep = ",")
                        } else if (name == "contains") {
                            if (url$query[[i]] == "1") {
                                query[["coordinate_modifier"]] <-
                                    Coordinates$Within
                            }
                        } else if (name == "exact") {
                            if (url$query[[i]] == "1") {
                                query[["coordinate_modifier"]] <-
                                    Coordinates$Exact
                            }
                        } else if (name == "either") {
                            if (url$query[[i]] == "1") {
                                query[["coordinate_modifier"]] <-
                                    Coordinates$StartIsExactOrWithin
                            } else if (url$query[[i]] == "2") {
                                query[["coordinate_modifier"]] <-
                                    Coordinates$EndIsExactOrWithin
                            }
                        } else {
                            query[[name]] <- c(query[[name]], url$query[[i]])
                        }
                    }

                    query$compilation <- strsplit(url$path, "/")[[1]][[1]]
                    private$query <- query

                    invisible(self)
                },
                print = function() {
                    cat("<SnaptronQueryBuilder>\n")
                    for (param in names(private$query)) {
                        if (is.null(private$query[[param]])) {
                            next
                        } else if (param == "coordinate_modifier") {
                            desc <- switch(
                                private$query[[param]],
                                Exact = "exact",
                                Within = "contains",
                                StartIsExactOrWithin = "either=1",
                                EndIsExactOrWithin = "either=2",
                                "overlaps"
                            )
                        } else {
                            desc <-
                                paste(private$query[[param]], collapse = ",")
                        }
                        cat("   ", param, ": ", desc, "\n", sep = "")
                    }
                }
            ),
            private = list(
                query = list(),
                call = function(fn_name, args) {
                    fn <- get(fn_name, parent.frame())
                    arg_names <- intersect(names(formals(fn)), names(args))
                    do.call(fn, args[arg_names])
                }
            )
            )


#' Query Junctions/Genes/Exons
#'
#' Given one or more gene names or genomic range intervals it will return a list
#' of 0 or more genes, junctions, or exons (depending on which query form is
#' used) which overlap the ranges.
#'
#' @param compilation A single string containing the name of the Snaptron
#'   datasource
#'
#' @param regions Either a list of 1 or more `HUGO` gene names `(e.g. "BRCA1")`
#'   or a GRanges-class object containing one or more genomic intervals
#'   `(e.g. "chr1:1-1000")`.  Strand information is ignored.
#'
#' @param range_filters A list of strings defining range-related contraints
#'
#' @param sample_filters A list of strings defining sample-related contraints
#'
#' @param sids A list of rail_ids (integer sample IDs) to filter results
#'   on. Only records which have been found in at least one of these samples
#'   will be returned.
#'
#' @param coordinate_modifier Snaptron coordinate modifier enum. Invariants
#'   include:
#'
#'   Coordinate$Exact - Contraints the results so that the start/end coordinates
#'   match the start/end of the specifiied range.
#'
#'   Coordinate$Within - Contraints the results so that the coordinates are
#'   within (inclusive) the specified range.
#'
#'   Coordinate$StartIsExactOrWithin - Constraints the results so that the start
#'   coorindate matches or is within the boundaries of the specified range.
#'
#'   Coorindate$EndIsExactOrWithin - Contraints the results so that that the
#'   end coordinate matches or is within the boundaries of the specified range.
#'
#'   Coordinate$Overlaps - Contraints the results so that the coorindates
#'   overlaps the specified range.
#'
#' @param return_rse Should the query data be returned as a simple data frame or
#'   converted to a RangeSummarizedExperiment.
#'
#' @param split_by_region By default the results from multiple queries will be
#'   returned in `RangedSummarizedExperiment` object with a `rowData` entry for
#'   each labeling each result row according to the query it resulted
#'   from. However, if this is set to `TRUE`, the result will be a list of
#'   RangedSummarizedExperiment objects, one per original interval/gene. This
#'   latter option may be useful but requires a separate copy of the sample
#'   metadata for each original interval/gene.
#'
#' @examples
#' query_jx(Compilation$gtex, "chr1:1-100000",
#'          range_filters = "samples_count >= 20")
#' query_gene(Compilation$gtex, "CD99", sample_filters = SMTS == "Brain")
#'
#' @return Functions will return either a RangeSummarizedexperiment or
#'   data.table depending on whether the \code{return_rse} parameter is set to
#'   \code{TRUE} or \code{FALSE}.
#' @export
query_jx <-function(compilation, regions, range_filters = NULL,
                    sample_filters = NULL, sids = NULL,
                    coordinate_modifier = NULL, return_rse = TRUE,
                    split_by_region = FALSE)
{
    strands <- NULL
    range_filters <- bool_expressions_to_strings(rlang::enexpr(range_filters))
    sample_filters <- bool_expressions_to_strings(rlang::enexpr(sample_filters))

    if (class(regions) == "GRanges") {
        strands <- extract_strands(regions)
        regions <- extract_intervals(regions)

        assert_that(length(regions) == length(strands))
    }

    should_bind <- length(regions) > 1 && !split_by_region
    res <- lapply(seq_along(regions), function(i) {
        if (!is.null(strands) && (strands[[i]] == "+" || strands[[i]] == "-")) {
            pos <- grep("strand", range_filters)
            if (!identical(pos, integer(0))) {
                range_filters[[pos]] <- paste0("strand:", strands[[i]])
            } else {
                range_filters <-
                    c(range_filters, paste0("strand:", strands[[i]]))
            }
        }

        run_query(compilation = compilation,
                  regions = regions[[i]],
                  range_filters = range_filters,
                  sample_filters = sample_filters,
                  coordinate_modifier = coordinate_modifier,
                  sids = sids,
                  return_rse = return_rse)
    })

    if (length(res) == 1) {
        res <- res[[1]]
    }
    if (should_bind) {
        rbind_func <- if (return_rse) SummarizedExperiment::rbind else rbind
        res <- do.call(rbind_func, res)
    }

    res
}

#' @rdname query_jx
#' @export
query_gene <- function(compilation, regions,
                       range_filters = NULL, sample_filters = NULL,
                       sids = NULL, coordinate_modifier = NULL,
                       return_rse = TRUE, split_by_region = FALSE)
{
    strands <- NULL
    range_filters <- bool_expressions_to_strings(rlang::enexpr(range_filters))
    sample_filters <- bool_expressions_to_strings(rlang::enexpr(sample_filters))

    if (class(regions) == "GRanges") {
        strands <- extract_strands(regions)
        regions <- extract_intervals(regions)

        assert_that(length(regions) == length(strands))
    }

    should_bind <- length(regions) > 1 && !split_by_region
    res <- lapply(seq_along(regions), function(i) {
        if (!is.null(strands) && (strands[i] == "+" || strands[[i]] == "-")) {
            pos <- grep("strand", range_filters)
            if (!identical(pos, integer(0))) {
                range_filters[pos] <- paste0("strand:", strands[[i]])
            } else {
                range_filters <-
                    c(range_filters, paste0("strand:", strands[[i]]))
            }
        }

        run_query(compilation = compilation,
                  regions = regions[[i]],
                  range_filters = range_filters,
                  sample_filters = sample_filters,
                  coordinate_modifier = coordinate_modifier,
                  sids = sids,
                  endpoint = "genes",
                  return_rse = return_rse)
    })
    if (length(res) == 1) {
        res <- res[[1]]
    }
    if (should_bind) {
        rbind_func <- if (return_rse) SummarizedExperiment::rbind else rbind
        res <- do.call(rbind_func, res)
    }

    res
}

#' @rdname query_jx
#' @export
query_exon <- function(compilation, regions,
                       range_filters = NULL, sample_filters = NULL,
                       sids = NULL, coordinate_modifier = NULL,
                       return_rse = TRUE, split_by_region = FALSE)
{
    strands <- NULL
    range_filters <- bool_expressions_to_strings(rlang::enexpr(range_filters))
    sample_filters <- bool_expressions_to_strings(rlang::enexpr(sample_filters))

    if (class(regions) == "GRanges") {
        strands <- extract_strands(regions)
        regions <- extract_intervals(regions)

        assert_that(length(regions) == length(strands))
    }

    should_bind <- length(regions) > 1 && !split_by_region
    res <- lapply(seq_along(regions), function(i) {
        if (!is.null(strands) && (strands[[i]] == "+" || strands[[i]] == "-")) {
            pos <- grep("strand", range_filters)
            if (!identical(pos, integer(0))) {
                range_filters[[pos]] <- paste0("strand:", strands[[i]])
            } else {
                range_filters <-
                    c(range_filters, paste0("strand:", strands[[i]]))
            }
        }

        run_query(
            compilation = compilation,
            regions = regions[i],
            range_filters = range_filters,
            sample_filters = sample_filters,
            coordinate_modifier = coordinate_modifier,
            sids = sids,
            endpoint = "exons",
            return_rse = return_rse
        )
    })
    if (length(res) == 1) {
        res <- res[[1]]
    }
    if (should_bind) {
        rbind_func <- if (return_rse) SummarizedExperiment::rbind else rbind
        res <- do.call(rbind_func, res)
    }

    res
}

#' Return the URI of the last successful request to Snaptron
#'
#' @description This function can be paired with the \code{from_url} method from
#'   the SnaptronQueryBuilder class, allowing users to share sources of data
#'   from Snaptron.
#' @return URI of last successful request to Snaptron or \code{NULL} if there
#'   have not been any successful requests.
#'
#' @export
#' @examples
#' query_jx(compilation = "gtex", regions = "CD99")
#' uri_of_last_successful_request()
uri_of_last_successful_request <- function() {
    pkg_globals$last_uri_accessed
}

get_compilation_metadata <- function(compilation) {
    assert_that(compilation %in% names(Compilation),
                msg = "Invalid compilation")

    if (is.null(pkg_globals$metadata[[compilation]])) {
        uri <- sprintf(
            "%s/%s/samples?all=1", pkg_globals$snaptron_host, compilation
        )
        tsv <- submit_query(uri)
        pkg_globals$metadata[[compilation]] <-
            data.table::fread(tsv, sep = "\t", quote = "")
    }

    pkg_globals$metadata[[compilation]]
}

is_logical_op <- function(op) {
    logical_ops <-
        c(
            as.symbol("=="),
            as.symbol("<="),
            as.symbol("<"),
            as.symbol(">"),
            as.symbol(">=")
        )
    any(vapply(logical_ops, identical, FUN.VALUE = logical(1), op))
}

tidy_filters <- function(filters) {
    filters <- gsub("==", ":", filters)
    filters <- gsub("=", ":", filters)
    filters <- gsub("\\s+", "", filters)

    filters
}

run_query <- function(compilation, regions, endpoint = "snaptron",
                      range_filters = NULL, sample_filters = NULL, sids = NULL,
                      coordinate_modifier = NULL, construct_rse = TRUE,
                      return_rse = TRUE) {
    uri <- generate_snaptron_uri(
            compilation = compilation,
            regions = regions,
            endpoint = endpoint,
            range_filters = range_filters,
            sample_filters = sample_filters,
            coordinate_modifier = coordinate_modifier,
            sids = sids
    )

    if (!is.null(tc <- getOption("test_context")) && tc == TRUE) {
        assign("last_uri_accessed", uri, pkg_globals)
        return(NULL)
    } else {
        tsv <- submit_query(uri)
        assign("last_uri_accessed", uri, pkg_globals)
    }

    query_data <- data.table::fread(tsv, sep = "\t", header = TRUE)
    if (nrow(query_data) == 0) {
        return(NULL)
    }

    if (!return_rse) {
        return(query_data)
    }

    metadata <- get_compilation_metadata(compilation)

    ## if (!is.null(sids)) {
    ##     metadata <- metadata[metadata$rail_id %in% sids,]
    ## }

    if (construct_rse == FALSE) {
        return(list(query_data = query_data, metadata = metadata))
    }

    rse(query_data, metadata)
}

generate_snaptron_uri <- function(compilation, regions,
                                  endpoint = "snaptron", range_filters = NULL,
                                  sample_filters = NULL,
                                  coordinate_modifier = NULL, sids = NULL) {
    assert_that(
        compilation %in% names(Compilation),
        msg = "Invalid compilation"
    )

    query <- ""
    path <- paste(compilation, paste0(endpoint, "?"), sep = "/")

    assert_that(is_hugo_gene(regions) || is_chromosome_interval(regions))
    if (!missing(regions)) {
        query <- paste("regions", regions, sep = "=")
    } else {
        stop("please specify either a gene or an interval", call. = FALSE)
    }

    if (!is.null(range_filters)) {
        query <-
            c(query, paste("rfilter", tidy_filters(range_filters), sep = "="))
    }

    if (!is.null(sample_filters)) {
        sample_filters <- tidy_filters(sample_filters)
        errors <- lapply(sample_filters, function(filter) {
            fields <- stringr::str_split(filter, "\\W", n = 2)[[1]]
            validate_sample_filter(
                compilation,
                name = fields[[1]],
                value = fields[[2]]
            )
        }) %>% purrr::compact()

        if (length(errors) > 0) {
            error_string <- paste(errors, collapse = "\n")
            stop(error_string, call. = FALSE)
        }
        query <-
            c(query, paste("sfilter", tidy_filters(sample_filters), sep = "="))
    }

    if (!is.null(coordinate_modifier)) {
        if (coordinate_modifier == Coordinates$Exact) {
            query <- c(query, paste("exact", "1", sep = "="))
        } else if (coordinate_modifier == Coordinates$Within) {
            query <- c(query, paste("contains", "1", sep = "="))
        } else if (coordinate_modifier == Coordinates$StartIsExactOrWithin) {
            query <- c(query, paste("either", "1", sep = "="))
        } else if (coordinate_modifier == Coordinates$EndIsExactOrWithin) {
            query <- c(query, paste("either", "2", sep = "="))
        } else {
            stop("Invalid coordinate modifier", stop. = FALSE)
        }
    }

    if (!is.null(sids)) {
        assert_that(is.wholenumber(sids), msg = "sids should be whole numbers")
        query <- c(query, paste("sids", paste(sids, collapse = ","), sep = "="))
    }

    paste0(pkg_globals$snaptron_host, path, paste(query, collapse = "&"))
}

extract_intervals <- function(g) {
    chr <- GenomicRanges::seqnames(g)
    beg <- GenomicRanges::start(g)
    end <- GenomicRanges::end(g)
    chr <- as.vector(rep(chr@values, chr@lengths))

    paste0(chr, ":", beg, "-", end)
}

extract_strands <- function(g) {
    strands <- GenomicRanges::strand(g)
    strands <- as.vector(rep(strands@values, strands@lengths))

    strands
}

submit_query <- function(uri) {
    resp <- httr::GET(uri)
    if (resp$status_code != "200"
        || resp$header[["content-type"]] != "text/plain") {
        stop("API did not return tsv", call. = FALSE)
    }
    rawToChar(resp$content)
}

convert_to_sparse_matrix <- function(rail_ids, unique_rail_ids, counts,
                                     samples_count, snaptron_ids) {
    i <- rep(seq_along(samples_count), samples_count)
    j <- match(rail_ids, unique_rail_ids)
    x <- counts

    dims <- c(length(snaptron_ids), length(unique_rail_ids))
    Matrix::sparseMatrix(
                i = i, j = j, x = x, dims = dims,
                dimnames = list(snaptron_ids, unique_rail_ids)
            )
}

get_counts <- function(query_data) {
    rail_ids_and_counts <- extract_samples(query_data) %>%
        strsplit(":", fixed = TRUE)
    rail_ids <- as.numeric(vapply(rail_ids_and_counts, `[`, 1, FUN.VALUE = ""))
    unique_rail_ids <- sort(rail_ids) %>% unique()
    counts <- as.numeric(vapply(rail_ids_and_counts, `[`, 2, FUN.VALUE = ""))

    sparse_matrix <-
        convert_to_sparse_matrix(
            rail_ids, unique_rail_ids, counts,
            query_data$samples_count, query_data$snaptron_id
        )

    list(sparse_matrix, unique_rail_ids)
}

get_col_data <- function(metadata, rail_ids = NULL) {
    if (!is.null(rail_ids)) {
        metadata[metadata$rail_id %in% rail_ids]
    } else {
        metadata
    }
}

get_row_ranges <- function(query_data) {
    cols <- c("chromosome", "start", "end", "length", "strand", "samples")
    mcols <- query_data[, !cols, with = FALSE]

    GRanges(
        seqnames = query_data$chromosome,
        IRanges::IRanges(query_data$start, query_data$end),
        strand = query_data$strand,
        mcols
    )
}

rse <- function(query_data, metadata) {
    row_ranges <- get_row_ranges(query_data)
    count_data <- get_counts(query_data)
    col_data <- get_col_data(metadata, rail_ids = count_data[[2]])

    SummarizedExperiment(
        assays = list(counts = count_data[[1]]),
        rowRanges = row_ranges,
        colData = col_data
    )
}
