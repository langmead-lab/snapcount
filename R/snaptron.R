metadata <- new.env()

URI <- NULL

#' A Reference Class for building Snaptron queries
#'
#' @method compilation Get/set snaptron
#' @method genes
#' @method intervals
#' @method range_filters
#' @method sample_filters
#' @method sids
#' @method query_jx
#' @method query_gene
#' @method query_exon
#' @method query_coverage
#' @export
#'
#' @examples
#' sb <- SnaptronQueryBuilder$new()
#' sb$compilation("srav2")$genes("CD99")$query_jx()
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
                private$query$gene_or_intervals
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
            if (!missing(range_filters)) {
                private$query$sample_filters <- sample_filters_filters
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

            resp <- httr::HEAD(url)
            if (resp$status_code != 200 || httr::http_type(resp) != "text/plain") {
                stop(sprintf("%s: is not a valid URL", url))
            }

            query <- list()
            for (i in 1:(length(url$query))) {
                name <- switch(names(url$query[i]),
                    rfilter = "range_filters",
                    sfilter = "sample_filters",
                    regions = "genes_or_intervals",
                    name)

                if (name == "sids") {
                    query[[name]] <- scan(url$query[[i]], sep = ",")
                }
                query[[name]] <- c(query[[name]], url$query[[i]])
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
            do.call(fn, private$query)
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
    uri <-
        generate_snaptron_uri(
            compilation = compilation,
            genes_or_intervals = genes_or_intervals,
            range_filters = range_filters,
            sample_filters = sample_filters,
            sids = sids
        )

    tsv <- submit_query(uri)
    URI <<- uri
    query_data <- data.table::fread(tsv, sep = '\t')
    metadata <- get_compilation_metadata(compilation)
    rse(query_data, metadata)
}

#' @rdname query_jx
#' @export
query_gene <- function(compilation, genes_or_intervals,
    range_filters = NULL, sample_filters = NULL, sids = NULL)
{
    uri <-
        generate_snaptron_uri(
            compilation = compilation,
            genes_or_intervals = genes_or_intervals,
            endpoint = "genes",
            range_filters = range_filters,
            sample_filters = sample_filters,
            sids = sids
        )

    tsv <- submit_query(uri)
    URI <<- uri
    query_data <- data.table::fread(tsv, sep = '\t')
    metadata <- get_compilation_metadata(compilation)
    rse(query_data, metadata)
}

#' @rdname query_jx
#' @export
query_exon <- function(compilation, genes_or_intervals,
    range_filters = NULL, sample_filters = NULL, sids = NULL)
{
    uri <-
        generate_snaptron_uri(
            compilation = compilation,
            genes_or_intervals = genes_or_intervals,
            endpoint = "exons",
            range_filters = range_filters,
            sample_filters = sample_filters,
            sids = sids
        )

    tsv <- submit_query(uri)
    URI <<- uri
    query_data <- data.table::fread(tsv, sep = '\t')
    metadata <- get_compilation_metadata(compilation)
    rse(query_data, metadata)
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
    uri <- generate_snaptron_uri(
        compilation = compilation,
        gene = genes_or_intervals,
        endpoint = "bases",
        sids = sids)

    tsv <- submit_query(uri)
    URI <<- uri
    query_data <- data.table::fread(tsv, sep = '\t')
    metadata <- get_compilation_metadata(compilation)

    rse(
        query_data,
        metadata,
        extract_col_data = coverage_col_data,
        extract_row_ranges = coverage_row_ranges,
        extract_counts = coverage_counts
    )
}

#' @export
uri_of_last_successful_request <- function() {
    URI
}

get_compilation_metadata <- function(compilation) {
    stopifnot(compilation %in% c("tcga", "srav2", "srav1", "gtex"))

    if (is.null(metadata[[compilation]])) {
        uri <- sprintf("http://snaptron.cs.jhu.edu/%s/samples?all=1", compilation)
        tsv <- submit_query(uri)
        metadata[[compilation]] <- data.table::fread(tsv, sep = '\t', quote = "")
    }

    metadata[[compilation]]
}

#' @export
exprs <- function(..., frame = as.list(parent.frame())) {
    filters <- lapply(rlang::exprs(...),
        function(e) {
            if (length(e) != 3) {
                err_msg <- paste0(deparse(e), ": is not a valid expression")
                stop(err_msg)
            }
            stopifnot(is_logical_op(e[[1]]))
            e[[3]] <- eval(methods::substituteDirect(e[[3]], frame = frame))
            deparse(e)
        }
    )
    simplify2array(filters)
}

is_logical_op <- function(op) {
    logical_ops <-
        c(
            as.symbol("=="),
            as.symbol("="),
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
    sorted_rail_ids <- sort(unique(rail_ids))

    i <- rep(seq_along(samples_count), samples_count)
    j <- match(rail_ids, sorted_rail_ids)
    x <- as.numeric(vapply(rail_ids_and_counts, `[`, 2, FUN.VALUE = ""))

    dims <- c(length(snaptron_ids), length(compilation_rail_ids))

    Matrix::sparseMatrix(i = i, j = j, x = x, dimnames = list(snaptron_ids, compilation_rail_ids), dims = dims)
}

counts <- function(query_data, metadata) {
    samples <- extract_samples(query_data)
    compilation_rail_ids <- metadata$rail_id
    convert_to_sparse_matrix(samples, query_data$samples_count, query_data$snaptron_id, compilation_rail_ids)
}

col_data <- function(metadata, sids) {
    metadata[metadata$rail_id %in% sids, ]
}

row_ranges <- function(query_data) {
    print(query_data)
    mcols <- subset(query_data,
        select = -c(chromosome, start, end, length, strand, samples))

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

coverage_counts <- function(query_data) {
    data <- query_data[, -(1:4)]
    Matrix::Matrix(as.matrix(data), sparse = TRUE)
}

coverage_col_data <- function(metadata, sids) {
    metadata[metadata$rail_id %in% sids, ]
}

rse <- function(query_data, metadata, extract_counts = counts, extract_row_ranges = row_ranges, extract_col_data = col_data) {
    row_ranges <- extract_row_ranges(query_data)
    counts <- extract_counts(query_data, metadata)
    # col_data <- extract_col_data(metadata, colnames(counts))
    col_data <- metadata

    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        rowRanges = row_ranges,
        colData = col_data
    )
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
}
