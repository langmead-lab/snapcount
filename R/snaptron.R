pkg_env <- new.env(parent = emptyenv())

pkg_env$URI <- NULL

## enum <- function(...) {
##     values <- sapply(match.call(expand.dots = TRUE)[-1L], deparse)

##     stopifnot(identical(unique(values), values))

##     res <- setNames(seq_along(values), values)
##     res <- as.environment(as.list(res))
##     lockEnvironment(res, bindings = TRUE)
##     res
## }

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
        query_jx = function() {
            private$call("query_jx", private$query)
        },
        query_exon = function() {
            private$call("query_exon", private$query)
        },
        query_gene = function() {
            private$call("query_gene", private$query)
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
                } else {
                    query[[name]] <- c(query[[name]], url$query[[i]])
                }
            }

            query$compilation <- strsplit(url$path, '/')[[1]][1]
            private$query <- query
        },
        print = function() {
            cat("<SnaptronQueryBuilder>\n")
            for (param in names(private$query)) {
                if (is.null(private$query[[param]])) {
                    next
                }

                cat("  ", param, ": ", paste(private$query[[param]], collapse = ','), '\n', sep = '')
            }
        }
    ),
    private = list(
        query = list(),
        call = function(fn_name, args) {
            fn <- get(fn_name, parent.frame())
            args <- intersect(names(formals(fn)), names(private$query))
            do.call(fn, private$query[args])
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

#' @export
query_jx <- function(compilation, genes_or_intervals, range_filters = NULL,
                sample_filters = NULL, sids = NULL)
{
    run_query(compilation = compilation, genes_or_intervals = genes_or_intervals,
        range_filters = range_filters, sample_filters = sample_filters, sids = sids)
}

#' @rdname query_jx
#' @export
query_gene <- function(compilation, genes_or_intervals,
    range_filters = NULL, sample_filters = NULL, sids = NULL)
{
    run_query(compilation = compilation, genes_or_intervals = genes_or_intervals,
        endpoint = "genes", range_filters = range_filters, sample_filters = sample_filters, sids = sids)
}

#' @rdname query_jx
#' @export
query_exon <- function(compilation, genes_or_intervals,
    range_filters = NULL, sample_filters = NULL, sids = NULL)
{
    run_query(compilation = compilation, genes_or_intervals = genes_or_intervals,
        endpoint = "exons", range_filters = range_filters, sample_filters = sample_filters, sids = sids)
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
    sample_filters = NULL, contains = FALSE, exact = FALSE, either = FALSE, sids = NULL, construct_rse = TRUE) {

    uri <-
        generate_snaptron_uri(
            compilation = compilation,
            genes_or_intervals = genes_or_intervals,
            endpoint = endpoint,
            range_filters = range_filters,
            sample_filters = sample_filters,
            contains = contains,
            exact = exact,
            either = either,
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
    metadata <- get_compilation_metadata(compilation)

    if (construct_rse == FALSE) {
        return(list(query_data = query_data, metadata = metadata))
    }

    rse(query_data, metadata)
}

generate_snaptron_uri <- function(compilation, genes_or_intervals, endpoint = "snaptron", range_filters = NULL,
                     sample_filters = NULL, contains = FALSE, exact = FALSE, either = FALSE, sids = NULL)
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

    if (contains) {
        query <- c(query, paste("contains", "1", sep = '='))
    }

    if (exact) {
        query <- c(query, paste("exact", "1", sep = '='))
    }

    if (either) {
        query <- c(query, paste("either", "1", sep = '='))
    }

    if (!is.null(sids)) {
        stopifnot(is.wholenumber(sids))
        query <- c(query, paste("sids", paste(sids, collapse = ','), sep = '='))
    }

    paste0(url, path, paste(query, collapse = '&'))
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
    data <- query_data[, -(1:4)]
    rail_ids <- as.numeric(colnames(data))
    smallest_rail_id <- metadata$rail_id[1]

    i <- rep(1:nrow(data), ncol(data))
    j <- rep((rail_ids - smallest_rail_id) + 1, each = nrow(data))

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
