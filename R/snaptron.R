metadata <- new.env()

#' @export
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
        genes = function(genes = NULL) {
            if (!missing(genes)) {
                private$query$gene <- genes
                invisible(self)
            } else {
                private$query$gene
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
            jx_query_fn <- get("query_jx", parent.frame())
            do.call(jx_query_fn, private$query)
        }
    ),
    private = list(
        query = list()
    )
)

#' @export
query_jx <- function(compilation, gene, interval, range_filters = NULL,
                sample_filters = NULL, sids = NULL)
{
    uri <-
        generate_snaptron_uri(
            compilation = compilation,
            gene = gene,
            interval = interval,
            range_filters = range_filters,
            sample_filters = sample_filters,
            sids = sids
        )

    tsv <- submit_query(uri)
    query_data <- data.table::fread(tsv, sep = '\t')
    metadata <- get_compilation_metadata(compilation)
    rse(query_data, metadata)
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
            deparse(methods::substituteDirect(e, frame = frame))
        }
    )
    tidy_filters(simplify2array(filters))
}

tidy_filters <- function(filters) {
    filters <- gsub("==", ":", filters)
    filters <- gsub("=", ":", filters)
    filters <- gsub("\\s+", "", filters)

    filters
}

generate_snaptron_uri <- function(compilation, gene, interval, range_filters = NULL,
                     sample_filters = NULL, contains = FALSE, exact = FALSE, either = FALSE, sids = NULL)
{
    url <- "http://snaptron.cs.jhu.edu/"

    if (missing(compilation)) {
        stop("compilation is a required argument")
    }

    stopifnot(compilation %in% c("tcga", "srav2", "srav1", "gtex"))

    path <- paste(compilation, "snaptron?", sep = '/')
    query <- ""

    if (!missing(gene)) {
        query <- paste("regions", gene, sep = '=')
    } else if (!missing(interval)) {
        query <- paste("regions", gene, sep = '=')
    } else {
        stop("please specify either a gene or an interval")
    }

    if (!is.null(range_filters)) {
        stopifnot(is.character(range_filters))

        query <- c(query, paste("rfilter", gsub('=', ':', range_filters), sep = '='))
    }

    if (!is.null(sample_filters)) {
        stopifnot(is.character(sample_filters))

        query <- c(query, paste("sfilter", gsub('=', ':', sample_filters), sep = '='))
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
        stopifnot(is.integer(sids))
        query <- c(query, paste("sids", paste(sids, collapse = ','), sep = '='))
    }

    paste0(url, path, paste(query, collapse = '&'))
}

submit_query <- function(uri) {
    resp <- curl::curl_fetch_memory(uri)
    if (resp$status_code != "200"
        && curl::parse_headers_list(resp$header)[["content-type"]] != "text/plain") {
        stop("API did not return tsv", call. = FALSE)
    }
    rawToChar(resp$content)
}

extract_samples <- function(query_data) {
    unlist(lapply(strsplit(query_data$samples, ","), `[`, -1))
}

convert_to_sparse_matrix <- function(samples, samples_count, snaptron_ids) {
    rail_ids_and_counts <- strsplit(samples, ':', fixed = TRUE)

    rail_ids <- as.numeric(vapply(rail_ids_and_counts, `[`, 1, FUN.VALUE = ""))
    sorted_rail_ids <- sort(unique(rail_ids))

    i <- rep(seq_along(samples_count), samples_count)
    j <- match(rail_ids, sorted_rail_ids)
    x <- as.numeric(vapply(rail_ids_and_counts, `[`, 2, FUN.VALUE = ""))

    Matrix::sparseMatrix(i = i, j = j, x = x, dimnames = list(snaptron_ids, sorted_rail_ids))
}

counts <- function(query_data) {
    samples <- extract_samples(query_data)
    convert_to_sparse_matrix(samples, query_data$samples_count, query_data$snaptron_id)
}

rse <- function(query_data, metadata) {
    mcols <- subset(query_data,
        select = -c(chromosome, start, end, length, strand, samples))

    row_ranges <- GenomicRanges::GRanges(
        seqnames = query_data$chromosome,
        IRanges::IRanges(query_data$start, query_data$end),
        strand = query_data$strand,
        mcols
    )

    counts <- counts(query_data)

    col_data <- metadata[metadata$rail_id %in% colnames(counts), ]

    SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts), rowRanges = row_ranges, colData = col_data)
}
