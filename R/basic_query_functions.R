pkg_env <- new.env(parent = emptyenv())

pkg_env$URI <- NULL


# Compilations <- enum(a, b, c, d)

# Compilations

#' A Reference Class for building Snaptron queries
#'
#'
#' @section Usage:
#' \preformatted{
#' sb <- SnaptronQueryBuilder$new()
#'
#' sb$compilation("tcga")
#' sb$range_filters(exprs(samples_count <= 10))$sids(500:600)
#'
#' sb$from_url("http://snaptron.cs.jhu.edu/tcga/snaptron?region=CD99")$query_jx()
#'
#' print(sb)
#'}
#'
#' @section Methods:
#' \code{compilation} Get/set snaptron
#'
#' \code{genes_or_intervals} Get/set genes or intervals
#'
#' \code{range_filters} Get or set range filters
#'
#' \code{sample_filters} Get or set sample filters
#'
#' \code{sids} Get or set sample ids
#'
#' \code{query_jx} call query_jx function
#'
#' \code{query_gene} call \code{\link{query_gene}} function
#'
#' \code{query_exon} call \code{\link{query_exon}} function
#'
#' \code{query_coverage} call \code{\link{query_coverage}} function
#'
#' \code{print} print output the state of the object
#'
#' \code{from_url} set attributes from url
#'
#' @seealso \code{\link{query_jx}}
#'
#' @export
#'
#' @examples
#' sb <- SnaptronQueryBuilder$new()
#' sb$compilation("srav2")$genes_or_intervals("CD99")$query_jx()
SnaptronQueryBuilder <- R6::R6Class("SnaptronQueryBuilder",
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
        genes_or_intervals = function(genes_or_intervals = NULL) {
            if (!missing(genes_or_intervals)) {
                private$query$genes_or_intervals <- genes_or_intervals
                invisible(self)
            } else {
                private$query$genes_or_intervals
            }
        },
        range_filters = function(range_filters = NULL) {
            if (!missing(range_filters)) {
                private$query$range_filters <- range_filters
                invisible(self)
            } else {
                private$query$range_filters
            }
        },
        sample_filters = function(sample_filters = NULL) {
            if (!missing(sample_filters)) {
                private$query$sample_filters <- sample_filters
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
            if (!missing(exact)) {
                private$query$coordinate_modifier <- coordinate_modifier
                invisible(self)
            } else {
                private$query$coordinate_modifier
            }
        },
        query_jx = function(return_rse = TRUE) {
            private$call("query_jx", c(list(return_rse = return_rse), private$query))
        },
        query_exon = function(return_rse = TRUE) {
            private$call("query_exon", c(list(return_rse = return_rse), private$query))
        },
        query_gene = function(return_rse = TRUE) {
            private$call("query_gene", c(list(return_rse = return_rse), private$query))
        },
        query_coverage = function() {
            private$call("query_coverage", private$query)
        },
        from_url = function(url) {
            url <- httr::parse_url(url)
            if (url$hostname != "snaptron.cs.jhu.edu") {
                stop("URL does not point to Snaptron server")
            }

            if (is.null(getOption("test_ctx"))) {
                resp <- httr::HEAD(url)
                if (resp$status_code != 200 || httr::http_type(resp) != "text/plain") {
                    stop(sprintf("%s: is not a valid URL", url))
                }
            }

            query <- list()
            for (i in 1:(length(url$query))) {
                name <- switch(n <- names(url$query[i]),
                    rfilter = "range_filters",
                    sfilter = "sample_filters",
                    regions = "genes_or_intervals",
                    n)

                if (name == "sids") {
                    query[[name]] <- scan(textConnection(url$query[[i]]), sep = ",")
                } else if (name == "contains") {
                    if (url$query[[i]] == "1") {
                        query[["coordinate_modifier"]] <- Coordinates$Within
                    }
                } else if (name == "exact") {
                    if (url$query[[i]] == "1") {
                        query[["coordinate_modifier"]] <- Coordinates$Exact
                    }
                } else if (name == "either") {
                    if (url$query[[i]] == "1") {
                        query[["coordinate_modifier"]] <- Coordinates$StartIsExactOrWithin
                    } else if (url$query[[i]] == "2") {
                        query[["coordinate_modifier"]] <- Coordinates$EndIsExactOrWithin
                    }
                } else {
                    query[[name]] <- c(query[[name]], url$query[[i]])
                }
            }

            query$compilation <- strsplit(url$path, '/')[[1]][1]
            private$query <- query

            invisible(self)
        },
        print = function() {
            cat("<SnaptronQueryBuilder>\n")
            for (param in names(private$query)) {
                if (is.null(private$query[[param]])) {
                    next
                }

                if (param == "coordinate_modifier") {
                    desc <- switch(private$query[[param]],
                                   `Coordinates$Exact` = "exact",
                                   `Coordinates$Within` = "contains",
                                   `Coordinates$StartIsExactOrWithin` = "either=1",
                                   `Coordinates$EndIsExactOrWithin` = "either=2",
                                   "overlaps")
                } else {
                    desc <- paste(private$query[[param]], collapse = ",")
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
#' Given one or more gene names or genomic range
#' intervals it will return a list of 0 or more genes, junctions, or exons
#' (depending on which query form is used) which overlap the ranges.
#' @param compilation A single string containing the name of the Snaptron datasource
#' @param genes_or_intervals Either a list of >=1 `HUGO` gene names `(e.g. "BRCA1")` or a
#'   GRanges-class object containing one or more genomic intervals `(e.g. "chr1:1-1000")`.
#'   Strand information is ignored.
#' @param range_filters A list of strings defining range-related contraints
#' @param sample_filters A list of strings defining sample-related contraints
#' @param sids A list of rail_ids (integer sample IDs) to filter results on. Only
#'   records which have been found in at least one of these samples will be returned.
#' @param either Contraints the results so that the start (either=1) or end
#'   (either=2) coordinate matches or is within the boundaries of the specified range.
#' @param contains Constraints the results so that the coordinates within (inclusive)
#'   the specified range.
#' @param exact Constraints the results so that the start/end coordinates
#'   match the start/end of the specified range.
#' @param return_rse Should the query data be returned as a simple data frame or
#    converted to a RangeSummarizedExperiment.

#' @export
query_jx <- function(compilation, genes_or_intervals, range_filters = NULL,
                sample_filters = NULL, sids = NULL, coordinate_modifier = NULL, return_rse = TRUE, split_by_region = FALSE)
{
    regions <- genes_or_intervals
    strands <- NULL

    if (class(genes_or_intervals) == "GRanges") {
        regions <- extract_intervals(genes_or_intervals)
        strands <- extract_strands(genes_or_intervals)

        stopifnot(length(regions) == length(regions))
    }

    should_bind = length(regions) > 1 && !split_by_region

    res <- lapply(1:length(regions), function(i) {
        if (!is.null(strands) && (strands[i] == "+" || strands[i] == "-")) {
            pos <- grep("strand", range_filters)
            if (!identical(pos, integer(0))) {
                range_filters[pos] <- paste0("strand:", strands[i])
            } else {
                range_filters <- c(range_filters, paste0("strand:", strands[i]))
            }
        }

        run_query(compilation = compilation,
                  genes_or_intervals = regions[i],
                  range_filters = range_filters,
                  sample_filters = sample_filters,
                  coordinate_modifier = coordinate_modifier,
                  sids = sids,
                  return_rse = return_rse)
    })

    if (should_bind) {
        rbind_func <- if (return_rse) SummarizedExperiment::rbind else rbind
        res <- do.call(rbind_func, res)
    }

    res
}

#' @rdname query_jx
#' @export
query_gene <- function(compilation, genes_or_intervals,
    range_filters = NULL, sample_filters = NULL, sids = NULL, coordinate_modifier = NULL, return_rse = TRUE, split_by_region = FALSE)
{
    regions <- genes_or_intervals
    strands <- NULL

    if (class(genes_or_intervals) == "GRanges") {
        regions <- extract_intervals(genes_or_intervals)
        strands <- extract_strands(genes_or_intervals)

        stopifnot(length(regions) == length(regions))
    }

    should_bind = length(regions) > 1 && !split_by_region

    res <- lapply(1:length(regions), function(i) {
        if (!is.null(strands) && (strands[i] == "+" || strands[i] == "-")) {
            pos <- grep("strand", range_filters)
            if (!identical(pos, integer(0))) {
                range_filters[pos] <- paste0("strand:", strands[i])
            } else {
                range_filters <- c(range_filters, paste0("strand:", strands[i]))
            }
        }

        run_query(compilation = compilation,
                  genes_or_intervals = regions[i],
                  range_filters = range_filters,
                  sample_filters = sample_filters,
                  coordinate_modifier = coordinate_modifier,
                  sids = sids,
                  endpoint = "genes",
                  return_rse = return_rse)
    })

    if (should_bind) {
        rbind_func <- if (return_rse) SummarizedExperiment::rbind else rbind
        res <- do.call(rbind_func, res)
    }

    res
}

#' @rdname query_jx
#' @export
query_exon <- function(compilation, genes_or_intervals,
    range_filters = NULL, sample_filters = NULL, sids = NULL, coordinate_modifier = NULL, return_rse = TRUE, split_by_region = FALSE)
{
    regions <- genes_or_intervals
    strands <- NULL

    if (class(genes_or_intervals) == "GRanges") {
        regions <- extract_intervals(genes_or_intervals)
        strands <- extract_strands(genes_or_intervals)

        stopifnot(length(regions) == length(regions))
    }

    should_bind = length(regions) > 1 && !split_by_region

    res <- lapply(1:length(regions), function(i) {
        if (!is.null(strands) && (strands[i] == "+" || strands[i] == "-")) {
            pos <- grep("strand", range_filters)
            if (!identical(pos, integer(0))) {
                range_filters[pos] <- paste0("strand:", strands[i])
            } else {
                range_filters <- c(range_filters, paste0("strand:", strands[i]))
            }
        }

        run_query(compilation = compilation,
                  genes_or_intervals = regions[i],
                  range_filters = range_filters,
                  sample_filters = sample_filters,
                  coordinate_modifier = coordinate_modifier,
                  sids = sids,
                  endpoint = "exons",
                  return_rse = return_rse)
    })

    if (should_bind) {
        rbind_func <- if (return_rse) SummarizedExperiment::rbind else rbind
        res <- do.call(rbind_func, res)
    }

    res
}

#' Query Coverage data
#'
#' This is the basic coverage query function.  Given one or more
#' gene names or genomic range intervals it will return a list of 1
#' or more bases (in 1-base pair intervals) with their list of coverage
#' counts across all samples `(default)` or a sub-selection of sample
#' columns.  This form doesn’t support any filters (except sample IDs `sids`) or modifiers.
#' NOTE: coordinates in this query form always half-open intervals, so their left
#' coordinate starts at 0 while the right coordinate starts at 1.
#' @inheritParams query_jx
#' @param group_names A list of one or more labels of the same length as the list of
#'   `genes_or_intervals`. These labels serve as demarcation sentinels for the output
#'   list of the bases since any one query will split over many output records `(typically)`.
#'   Not required, but highly recommended.
#' @param bulk Use the `Snaptron` bulk query interface. This will perform better when
#'   running many queries at once. There is a limit of 50 queries per unit. More than 50
#'   queries can be submitted but they will be broken up into units of 50.
#' @param split_by_region By default the results from multiple queries will be returned
#'   in `RangedSummarizedExperiment` object with a `rowData` entry for each labeling each
#'   result row according to the query it resulted from. However, if this is set to `TRUE`,
#'   the result will be a list of RangedSummarizedExperiment objects, one per original
#'   interval/gene. This latter option may be useful but reqires metadata for each original
#'   interval/gene.

#' @export
#' @examples
#' query_coverage("gtex", "BRCA1", sids = c(50099,50102,50113))
query_coverage <- function(compilation, genes_or_intervals, group_names = NULL,
    sids = NULL, bulk = FALSE, split_by_region = FALSE)
{
    data <- run_query(compilation = compilation, genes_or_intervals = genes_or_intervals,
        endpoint = "bases", sids = sids, construct_rse = FALSE)

    if (is.null(getOption("test_ctx"))) {
        rse(
            data$query_data,
            data$metadata,
            extract_row_ranges = coverage_row_ranges,
            extract_counts = coverage_counts
        )
    }
}
#' Return the URI of the last successful request to Snaptron
#'
#' @description
#' This function can be paired with the \code{from_url} method from the
#' SnaptronQueryBuilder class, allowing users to share sources of
#' data from Snaptron.
#' @return
#' URI of last successful request to Snaptron or \code{NULL} if there have
#' not been any successful requests.
#'
#' @export
#' @examples
#' query_jx(compilation = "gtex", genes_or_intervals = "CD99")
#' uri_of_last_successful_request()
uri_of_last_successful_request <- function() {
    pkg_env$URI
}

get_compilation_metadata <- function(compilation) {
    stopifnot(compilation %in% c("tcga", "srav2", "srav1", "gtex"))

    if (is.null(pkg_env$metadata[[compilation]])) {
        uri <- sprintf("http://snaptron.cs.jhu.edu/%s/samples?all=1", compilation)
        tsv <- submit_query(uri)
        pkg_env$metadata[[compilation]] <- data.table::fread(tsv, sep = '\t', quote = "")
    }

    pkg_env$metadata[[compilation]]
}

#' Utility function for formatting range and sample filters
#'
#' @description
#' This function makes it easier to compose sample and range filters for Snaptron queries.
#' Each expression is expected to be a boolean expression.
#'
#' @param ... A list of expression to be converted to strings
#' @param frame An environment or list object
#'
#' @export
#' @examples
#' a <- 10
#' exprs(sample_filters < a)
#' exprs(samples_count == a + 10)
exprs <- function(..., frame = as.list(parent.frame())) {
    filters <- lapply(rlang::exprs(...),
        function(e) {
            if (length(e) != 3) {
                err_msg <- paste0(deparse(e), ": is not a valid expression")
                stop(err_msg)
            }

            stopifnot(is_logical_op(e[[1]]))

            if (is.character(e[[3]])) {
                e[[3]] <- as.symbol(e[[3]])
            } else {
                e[[3]] <- eval(methods::substituteDirect(e[[3]], frame = frame))
            }

            deparse(e, backtick = FALSE)
        }
    )
    simplify2array(filters)
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

run_query <- function(compilation, genes_or_intervals, endpoint = "snaptron", range_filters = NULL,
    sample_filters = NULL, sids = NULL, coordinate_modifier = NULL, construct_rse = TRUE, return_rse = TRUE) {

    uri <-
        generate_snaptron_uri(
            compilation = compilation,
            genes_or_intervals = genes_or_intervals,
            endpoint = endpoint,
            range_filters = range_filters,
            sample_filters = sample_filters,
            coordinate_modifier = coordinate_modifier,
            sids = sids
        )

    if (!is.null(getOption("test_ctx"))) {
        pkg_env$URI <- uri
        return(NULL)
    } else {
        tsv <- submit_query(uri)
        pkg_env$URI <- uri
    }

    query_data <- data.table::fread(tsv, sep = '\t')
    if (nrow(query_data) == 0) {
        stop("Query returned 0 rows", call. = FALSE)
    }

    if (!return_rse) {
        return(query_data)
    }

    metadata <- get_compilation_metadata(compilation)

    if (!is.null(sids)) {
        metadata <- metadata[metadata$rail_id %in% sids,]
    }

    if (construct_rse == FALSE) {
        return(list(query_data = query_data, metadata = metadata))
    }

    rse(query_data, metadata)
}

generate_snaptron_uri <- function(compilation, genes_or_intervals, endpoint = "snaptron", range_filters = NULL,
                     sample_filters = NULL, coordinate_modifier = NULL, sids = NULL)
{
    url <- "http://snaptron.cs.jhu.edu/"

    if (missing(compilation)) {
        stop("compilation is a required argument")
    }

    stopifnot(compilation %in% c("tcga", "srav2", "srav1", "gtex"))

    path <- paste(compilation, paste0(endpoint, "?"), sep = '/')
    query <- ""

    if (!missing(genes_or_intervals)) {
        query <- paste("regions", genes_or_intervals, sep = '=')
    } else {
        stop("please specify either a gene or an interval")
    }

    if (!is.null(range_filters)) {
        stopifnot(is.character(range_filters))

        query <- c(query, paste("rfilter", tidy_filters(range_filters), sep = '='))
    }

    if (!is.null(sample_filters)) {
        stopifnot(is.character(sample_filters))

        query <- c(query, paste("sfilter", tidy_filters(sample_filters), sep = '='))
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
            stop("Invalid coordinate filter")
        }
    }

    if (!is.null(sids)) {
        stopifnot(is.wholenumber(sids))
        query <- c(query, paste("sids", paste(sids, collapse = ','), sep = '='))
    }

    paste0(url, path, paste(query, collapse = '&'))
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

extract_samples <- function(query_data) {
    unlist(lapply(strsplit(query_data$samples, ","), `[`, -1))
}

convert_to_sparse_matrix <- function(samples, samples_count, snaptron_ids, compilation_rail_ids) {
    rail_ids_and_counts <- strsplit(samples, ':', fixed = TRUE)

    rail_ids <- as.numeric(vapply(rail_ids_and_counts, `[`, 1, FUN.VALUE = ""))

    i <- rep(seq_along(samples_count), samples_count)
    j <- match(rail_ids, compilation_rail_ids)
    x <- as.numeric(vapply(rail_ids_and_counts, `[`, 2, FUN.VALUE = ""))

    dims <- c(length(snaptron_ids), length(compilation_rail_ids))

    Matrix::sparseMatrix(i = i, j = j, x = x, dimnames = list(snaptron_ids, compilation_rail_ids), dims = dims)
}

counts <- function(query_data, metadata) {
    samples <- extract_samples(query_data)
    convert_to_sparse_matrix(samples, query_data$samples_count, query_data$snaptron_id, metadata$rail_id)
}

col_data <- function(metadata, sids = NULL) {
    metadata
}

#' @import data.table
row_ranges <- function(query_data) {
    cols <- c("chromosome", "start", "end", "length", "strand", "samples")
    mcols <- query_data[, !cols, with = FALSE]

    GenomicRanges::GRanges(
        seqnames = query_data$chromosome,
        IRanges::IRanges(query_data$start, query_data$end),
        strand = query_data$strand,
        mcols
    )
}

coverage_row_ranges <- function(query_data) {
    GenomicRanges::GRanges(
        seqnames = query_data$chromosome,
        IRanges::IRanges(query_data$start, query_data$end)
    )
}

coverage_counts <- function(query_data, metadata) {
    data <- query_data[, -c("DataSource:Type", "chromosome", "start", "end")]
    rail_ids <- as.numeric(colnames(data))
    smallest_rail_id <- metadata$rail_id[1]

    i <- rep(1:nrow(data), ncol(data))

    # did the user specify a list of sids?
    if (nrow(metadata) == ncol(data)) {
        j <- rep(1:ncol(data), each = nrow(data))
    } else {
        j <- rep((rail_ids - smallest_rail_id) + 1, each = nrow(data))
    }

    m <- base::as.matrix(data)
    dim(m) <- c(nrow(m) * ncol(m), 1)

    dims <- c(nrow(query_data), length(metadata$rail_id))
    Matrix::sparseMatrix(i = i , j = j, x = as.numeric(m[, 1]), dimnames = list(NULL, metadata$rail_id), dims = dims)
}

rse <- function(query_data, metadata, extract_counts = counts, extract_row_ranges = row_ranges, extract_col_data = col_data) {
    row_ranges <- extract_row_ranges(query_data)
    counts <- extract_counts(query_data, metadata)
    col_data <- extract_col_data(metadata)

    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        rowRanges = row_ranges,
        colData = col_data
    )
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
}
