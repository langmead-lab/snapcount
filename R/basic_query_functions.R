#' Query Junctions/Genes/Exons
#'
#' Given one or more gene names or genomic range
#' intervals it will return a list of 0 or more genes, junctions, or exons
#' (depending on which query form is used) which overlap the ranges.
#'
#' @param sb A SnaptonQueryBuilder object
#' @param return_rse Should the query data be returned as a simple data frame or
#'   converted to a RangeSummarizedExperiment.
#'
#' @param split_by_region By default the results from multiple queries will be returned
#'   in `RangedSummarizedExperiment` object with a `rowData` entry for each labeling each
#'   result row according to the query it resulted from. However, if this is set to `TRUE`,
#'   the result will be a list of RangedSummarizedExperiment objects, one per original
#'   interval/gene. This latter option may be useful but requires a separate copy of the
#'   sample metadata for each original interval/gene.
#'
#' @examples
#' # Contruct SnaptronQueryBuilder object using wrapper functions
#' qb <- QueryBuilder(compilation = "gtex", regions = "chr1:1-100000")
#' qb <- set_row_filters(samples_count >= 20)
#' query_jx(qb)
#'
#' qb <- set_row_filters(NULL)
#' qb <- set_column_filters(SMTS == "Brain")
#' query_gene(qb)
#'
#' # or directly using R6
#' sb <- SnaptronQueryBuilder$new()
#' sb$compilation("gtex")
#' sb$regions("chr1:1-100000")
#' sb$row_filters("samples_count >= 20")
#' query_jx(sb)
#'
#' sb$regions("CD99")
#' sb$row_filters(NULL)
#' sb$column_filters(SMTS == "Brain")
#' query_gene(sb)
#'
#' @return Functions will return either a RangeSummarizedexperiment or data.table depending
#'   on whether the \code{return_rse} parameter is set to \code{TRUE} or \code{FALSE}.

#' @export
query_jx <- function(sb, return_rse = TRUE, split_by_region = FALSE)
{
    strands <- NULL
    regions <- sb$regions()
    row_filters <- sb$row_filters()
    column_filters <- sb$column_filters()

    if (class(sb$regions()) == "GRanges") {
        strands <- extract_strands(sb$regions())
        regions <- extract_intervals(sb$regions())

        assert_that(length(regions) == length(strands))
    }

    should_bind <- length(regions) > 1 && !split_by_region
    res <- lapply(seq_along(regions), function(i) {
        if (!is.null(strands) && (strands[[i]] == "+" || strands[[i]] == "-")) {
            pos <- grep("strand", row_filters)
            if (!identical(pos, integer(0))) {
                row_filters[[pos]] <- paste0("strand:", strands[[i]])
            } else {
                row_filters <-
                    c(row_filters, paste0("strand:", strands[[i]]))
            }
        }

        run_query(compilation = sb$compilation(),
                  regions = regions[[i]],
                  row_filters = row_filters,
                  column_filters = column_filters,
                  coordinate_modifier = sb$coordinate_modifier(),
                  sids = sb$sids(),
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
query_gene <- function(sb, return_rse = TRUE, split_by_region = FALSE)
{
    strands <- NULL
    regions <- sb$regions()
    row_filters <- sb$row_filters()
    column_filters <- sb$column_filters()

    if (class(regions) == "GRanges") {
        strands <- extract_strands(sb$regions())
        regions <- extract_intervals(sb$regions())

        assert_that(length(regions) == length(strands))
    }

    should_bind <- length(regions) > 1 && !split_by_region
    res <- lapply(seq_along(regions), function(i) {
        if (!is.null(strands) && (strands[i] == "+" || strands[[i]] == "-")) {
            pos <- grep("strand", row_filters)
            if (!identical(pos, integer(0))) {
                row_filters[pos] <- paste0("strand:", strands[[i]])
            } else {
                row_filters <- c(row_filters, paste0("strand:", strands[[i]]))
            }
        }

        run_query(compilation = sb$compilation(),
                  regions = regions[[i]],
                  row_filters = row_filters,
                  column_filters = column_filters,
                  coordinate_modifier = sb$coordinate_modifier(),
                  sids = sb$sids(),
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
query_exon <- function(sb, return_rse = TRUE, split_by_region = FALSE)
{
    strands <- NULL
    regions <- sb$regions()
    row_filters <- sb$row_filters()
    column_filters <- sb$column_filters()

    if (class(regions) == "GRanges") {
        strands <- extract_strands(sb$regions())
        regions <- extract_intervals(sb$regions())


        assert_that(length(regions) == length(strands))
    }

    should_bind <- length(regions) > 1 && !split_by_region
    res <- lapply(seq_along(regions), function(i) {
        if (!is.null(strands) && (strands[i] == "+" || strands[i] == "-")) {
            pos <- grep("strand", row_filters)
            if (!identical(pos, integer(0))) {
                row_filters[pos] <- paste0("strand:", strands[i])
            } else {
                row_filters <- c(row_filters, paste0("strand:", strands[i]))
            }
        }

        run_query(compilation = sb$compilation(),
                  regions = regions[i],
                  row_filters = row_filters,
                  column_filters = column_filters,
                  coordinate_modifier = sb$coordinate_modifier(),
                  sids = sb$sids(),
                  endpoint = "exons",
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
#' query_jx(compilation = "gtex", regions = "CD99")
#' uri_of_last_successful_request()
uri_of_last_successful_request <- function() {
    pkg_globals$last_uri_accessed
}

get_compilation_metadata <- function(compilation) {
    assert_that(compilation %in% names(Compilation),
                msg = "Invalid compilation")

    if (is.null(pkg_globals$metadata[[compilation]])) {
        uri <- sprintf("%s/%s/samples?all=1",
                       pkg_globals$snaptron_host, compilation)
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
                      row_filters = NULL, column_filters = NULL, sids = NULL,
                      coordinate_modifier = NULL, construct_rse = TRUE,
                      return_rse = TRUE) {
    uri <-
        generate_snaptron_uri(
            compilation = compilation,
            regions = regions,
            endpoint = endpoint,
            row_filters = row_filters,
            column_filters = column_filters,
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
        warning(sprintf("query with uri: %s, returned no data.", uri))
        return(NULL)
    }

    if (!return_rse) {
        return(query_data)
    }

    metadata <- get_compilation_metadata(compilation)

    if (construct_rse == FALSE) {
        return(list(query_data = query_data, metadata = metadata))
    }

    rse(query_data, metadata)
}

generate_snaptron_uri <- function(compilation, regions,
                                  endpoint = "snaptron", row_filters = NULL,
                                  column_filters = NULL,
                                  coordinate_modifier = NULL, sids = NULL) {
    assert_that(compilation %in% names(Compilation),
                msg = "Invalid compilation")

    query <- ""
    path <- paste(compilation, paste0(endpoint, "?"), sep = "/")

    if (!missing(regions)) {
        query <- paste("regions", regions, sep = "=")
    } else {
        stop("please specify either a gene or an interval", call. = FALSE)
    }

    if (!is.null(row_filters)) {
        query <- c(query, paste("rfilter",
                                tidy_filters(row_filters), sep = "="))
    }

    if (!is.null(column_filters)) {
        column_filters <- tidy_filters(column_filters)
        errors <- lapply(column_filters, function(filter) {
            fields <- stringr::str_split(filter, "\\W", n = 2)[[1]]
            validate_sample_filter(compilation,
                                   name = fields[[1]],
                                   value = fields[[2]])
        }) %>% purrr::compact()

        if (length(errors) > 0) {
                error_string <- paste(errors, collapse = "\n")
                stop(error_string, call. = FALSE)
        }
        query <- c(query, paste("sfilter",
                                tidy_filters(column_filters),
                                sep = "="))
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
        assert_that(is.wholenumber(sids),
                    msg = "sids should be whole numbers")
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

convert_to_sparse_matrix <-
    function(rail_ids, unique_rail_ids, counts, samples_count, snaptron_ids) {
    i <- rep(seq_along(samples_count), samples_count)
    j <- match(rail_ids, unique_rail_ids)
    x <- counts

    dims <- c(length(snaptron_ids), length(unique_rail_ids))
    Matrix::sparseMatrix(i = i, j = j, x = x, dims = dims,
                         dimnames = list(snaptron_ids,
                                         paste0("rail_", unique_rail_ids)))
}

get_counts <- function(query_data) {
    rail_ids_and_counts <- extract_samples(query_data) %>%
        strsplit(":", fixed = TRUE)
    rail_ids <- as.numeric(vapply(rail_ids_and_counts, `[`, 1, FUN.VALUE = ""))
    unique_rail_ids <- sort(rail_ids) %>% unique()
    counts <- as.numeric(vapply(rail_ids_and_counts, `[`, 2, FUN.VALUE = ""))

    sparse_matrix <-
        convert_to_sparse_matrix(rail_ids, unique_rail_ids, counts,
                                 query_data$samples_count,
                                 query_data$snaptron_id)

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

    GenomicRanges::GRanges(
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

    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = count_data[[1]]),
        rowRanges = row_ranges,
        colData = col_data,
        metadata = list(uri = uri_of_last_successful_request())
    )
}
