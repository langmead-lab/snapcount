FormatForJunctionSeq <-
    R6Class(
        "FormatForJunctionSeq",
        public = list(
            initialize = function(query_builders, group_names, gene,
                                  sample_names = NULL) {
                private$query_builders <- query_builders
                private$group_names <- group_names
                private$gene <- gene
                private$sample_names <- sample_names
                private$design_ <-
                    data.frame(condition = factor(unlist(group_names)))
                private$create_gff_inputs()
            },

            format_for_junction_seq = function(env = parent.frame()) {
                obj_name <- private$get_object_name(env = env)
                if (is.null(obj_name)) {
                    obj_name <- "obj"
                }

                keys <- c("countdata", "samplenames", "design", "flat.gff.file")
                vals <- c("countdata", "samplenames",
                          "design", "flat.gff.file") %>%
                    paste(obj_name, ., sep = "$")
                args <- paste(keys, vals, sep = "=", collapse = ", ")
                f1 <- paste0("jscs <- JunctionSeq::readJunctionSeqCounts(",
                             args, ")")

                keys <- c("geneID", "jscs", "plot.type", "displayTranscripts",
                          "plot.exon.results", "plot.junction.results",
                          "plot.novel.junction.results",
                          "plot.untestable.results")
                vals <- c(paste(obj_name, "geneID", sep = "$"),
                          "jscs", "\"rawCounts\"", "TRUE",
                          "TRUE", "TRUE", "TRUE", "TRUE")
                args <- paste(keys, vals, sep = "=", collapse = ", ")
                f2 <- paste0("JunctionSeq::plotJunctionSeqResultsForGene(",
                             args, ")")

                cat(f1, f2, sep = "\n\n")
            },
            write_gff_file = function(filename = NULL) {
                if (is.null(filename)) {
                    private$gff_filename <- paste0(private$gene, ".gff")
                } else {
                    private$gff_filename <- filename
                }
                writeLines(private$gff, private$gff_filename)
            }
        ),
        active = list(
            flat.gff.file = function() {
                if (is.null(private$gff_filename)) {
                    stop(paste("please call 'write_gff_file' method",
                               "first with an optional filename"))
                }
                private$gff_filename
            },
            gff_data = function() {
                private$gff
            },
            countdata = function() {
                private$count_dataframes
            },

            design = function() {
                private$design_
            },
            samplenames = function() {
                if (is.null(private$sample_names)) {
                    private$sample_names <-
                        paste("SAMP", seq_along(private$group_names), sep = "")
                }
                private$sample_names
            },
            geneID = function() {
                private$gene
            }
        ),
        private = list(
            idx = 0,
            design_ = NULL,

            query_builders = list(),
            group_names = list(),
            gene = NULL,
            sample_names = NULL,

            junctions = list(),
            exons = list(),
            genes = list(),

            key_to_index = list(),
            key_to_type = list(),
            key_to_group = list(),

            group_counts = list(),
            count_dataframes = list(),

            gff = NULL,
            gff_filename = NULL,
            gff_format_str =
                paste(
                    "%s",
                    # chromosome
                    "ScalaUtils",
                    "%s",
                    # strand_type
                    "%d",
                    # chromosome start
                    "%d",
                    # chromosome end
                    ".",
                    "%s",
                    # strand
                    ".",
                    "gene_id %s; tx_set %s; num %03d; gene_set %s",
                    sep = "\t"
                ),
            counts_format_str =
                paste("%s:%s%03d",
                    "%d",
                    sep = "\t"),

            create_gff_inputs = function() {
                private$run_queries()
                private$process_genes()
                private$process_exons()
                private$process_junctions()
                private$fill_in_missing_counts()
                private$counts_to_dataframe()
            },

            run_queries = function() {
                lapply(seq_along(private$group_names), function(i) {
                    name <- private$group_names[[i]]
                    private$junctions[[name]] <-
                        private$query_builders[[i]]$query_jx(return_rse = FALSE)
                    if (is.null(private$junctions[[name]])) {
                        stop(paste("junction query with uri: ",
                                   uri_of_last_successful_request(),
                                   ", returned NULL"))
                    }
                    query_builder <-
                        private$query_builders[[i]]$clone()$row_filters(NULL)
                    private$exons[[name]] <-
                        query_builder$query_exon(return_rse = FALSE)
                    if (is.null(private$exons[[name]])) {
                        stop(
                            paste(
                                "exon query with uri: ",
                                uri_of_last_successful_request(),
                                ", returned NULL"
                            )
                        )
                    }
                    query_builder <-
                        private$query_builders[[i]]$clone()$row_filters(NULL)
                    private$genes[[name]] <-
                        query_builder$query_gene(return_rse = FALSE)
                    if (is.null(private$genes[[name]])) {
                        stop(
                            paste(
                                "gene query with uri: ",
                                uri_of_last_successful_request(),
                                ", returned NULL"
                            )
                        )
                    }
                })
            },

            process_junctions = function() {
                for (group in private$group_names) {
                    apply(private$junctions[[group]], 1, function(row) {
                        key <-
                            private$create_key(
                                        row[["chromosome"]], row[["start"]],
                                        row[["end"]], row[["strand"]])
                        if (key %in% names(private$key_to_group)) {
                            idx <- private$key_to_index[[key]]
                            type <- private$key_to_type[[key]]
                        } else {
                            private$idx <- private$idx + 1
                            private$key_to_index[[key]] <- private$idx
                            private$key_to_group[[key]] <- NULL
                            idx <- private$key_to_index[[key]]
                            type <- "J"
                            jx_type <- "splice_site"
                            tx_set <- "tx1"
                            gene_set <- private$gene

                            if (row[["annotated"]] == "0") {
                                type <- "N"
                                jx_type <-
                                    paste("novel", jx_type, sep = "_")
                                tx_set <- "UNKNOWN_TX"
                                gene_set <- "UNKNOWN_GENE_SET"
                            }

                            private$write_gff(
                                row[["chromosome"]],
                                jx_type,
                                row[["start"]],
                                row[["end"]],
                                row[["strand"]],
                                private$gene,
                                tx_set,
                                private$idx,
                                gene_set
                            )
                        }

                        private$write_count(
                                    group,
                                    private$gene,
                                    type,
                                    idx,
                                    row[["coverage_sum"]])
                        private$assign_group_to_key(group, key)
                        private$key_to_type[[key]] <- type
                    })
                }
            },

            process_exons = function() {
                for (group in private$group_names) {
                    apply(private$exons[[group]], 1, function(row) {
                        gene_name <-
                            strsplit(
                                row[["gene_id:gene_name:gene_type:bp_length"]],
                                split = ":"
                            )[[1]][[2]]
                        if (gene_name != private$gene) {
                            return()
                        }

                        key <-
                            private$create_key(
                                        row[["chromosome"]], row[["start"]],
                                        row[["end"]], row[["strand"]])
                        if (key %in% names(private$key_to_group)) {
                            idx <- private$key_to_index[[key]]
                            type <- private$key_to_type[[key]]
                        } else {
                            private$idx <- private$idx + 1
                            private$key_to_index[[key]] <- private$idx
                            private$key_to_group[[key]] <- NULL
                            idx <- private$key_to_index[[key]]
                            str_type <- "exonic_part"
                            type <- "E"

                            private$write_gff(
                                row[["chromosome"]],
                                str_type,
                                row[["start"]],
                                row[["end"]],
                                row[["strand"]],
                                private$gene,
                                "tx1",
                                idx,
                                private$gene
                            )
                        }

                        private$write_count(
                                    group,
                                    private$gene,
                                    type,
                                    idx,
                                    row[["coverage_sum"]])
                        private$assign_group_to_key(group, key)
                        private$key_to_type[[key]] <- type
                    })
                }
            },

            process_genes = function() {
                for (group in private$group_names) {
                    apply(private$genes[[group]], 1, function(row) {
                        gene_name <-
                            strsplit(
                                row[["gene_id:gene_name:gene_type:bp_length"]],
                                split = ":"
                            )[[1]][[2]]
                        if (gene_name != private$gene) {
                            return()
                        }

                        key <-
                            private$create_key(
                                        row[["chromosome"]],
                                        row[["start"]],
                                        row[["end"]],
                                        row[["strand"]])
                        if (key %in% names(private$key_to_group)) {
                            idx <- private$key_to_index[[key]]
                            type <- private$key_to_type[[key]]
                        } else {
                            private$key_to_index[[key]] <- private$idx
                            private$key_to_group[[key]] <- NULL
                            idx <- private$key_to_index[[key]]

                            type <- "A"
                            str_type <- "aggregate_gene"

                            private$write_gff(
                                row[["chromosome"]],
                                str_type,
                                row[["start"]],
                                row[["end"]],
                                row[["strand"]],
                                private$gene,
                                "tx1",
                                idx,
                                private$gene
                            )
                        }

                        private$write_count(
                                    group,
                                    private$gene,
                                    type,
                                    idx,
                                    row[["coverage_sum"]])
                        private$assign_group_to_key(group, key)
                        private$key_to_type[[key]] <- type
                    })
                }
            },

            write_gff = function(chromosome,
                type,
                start,
                end,
                strand,
                gene_id,
                tx_set,
                idx,
                gene_set) {
                gff_entry <- sprintf(
                    private$gff_format_str,
                    chromosome,
                    type,
                    as.numeric(start),
                    as.numeric(end),
                    strand,
                    gene_id,
                    tx_set,
                    idx,
                    gene_set
                )
                private$gff <-
                    c(private$gff, gff_entry)
            },

            write_count = function(group, gene_id, type, idx, coverage_sum) {
                private$group_counts[[group]] <- c(
                    private$group_counts[[group]],
                    sprintf(
                        private$counts_format_str,
                        gene_id,
                        type,
                        idx,
                        as.numeric(coverage_sum)
                    )
                )
            },

            create_key = function(...) {
                paste(c(...), collapse = "_")
            },

            assign_group_to_key = function(group, key) {
                if (is.null(private$key_to_group) ||
                    !any(is.element(private$key_to_group[[key]], group))) {
                    private$key_to_group[[key]] <-
                        c(private$key_to_group[[key]], group)
                }
            },

            counts_to_dataframe = function() {
                private$count_dataframes <-
                    lapply(private$group_counts, function(group) {
                        tsv <- paste(group, collapse = "\n")
                        dt <-
                            data.table::fread(tsv, sep = "\t", header = FALSE)
                        as.data.frame(data.table::setorder(dt, V1))
                    })
            },

            fill_in_missing_counts = function() {
                for (key in names(private$key_to_group)) {
                    groups <- private$key_to_group[[key]]
                    idx <- private$key_to_index[[key]]
                    type <- private$key_to_type[[key]]
                    for (group in setdiff(private$group_names, groups)) {
                        private$write_count(group, private$gene, type, idx, 0)
                    }
                }
            },

            get_object_name = function(env) {
                res <- ls(envir = env) %>%
                    purrr::keep(private$FormatForJunctionSeq_objs, env) %>%
                    purrr::keep(private$same_address_as_self, env)

                return(if (length(res) == 0) NULL else res[[1]])
            },

            FormatForJunctionSeq_objs = function(obj_name, env) {
                obj <- get(obj_name, envir = env)
                is(obj, "FormatForJunctionSeq")
            },

            same_address_as_self = function(obj_name, env) {
                obj <- get(obj_name, envir = env)
                address(obj) == address(self)
            }
        )
    )

#' Formats Snapcount results in a form that can be easily passed into
#' the third-party package, JunctionSeq.
#'
#' This allows for the visualization of junction counts along with
#' their associated gene and exon counts in the context of a gene
#' model. Performance is dependent on the total number of junctions
#' and the length of the gene model i.e. larger gene regions may
#' either completely fail to draw or take too long to be visualized.
#'
#' @param query_builders A list of 1 of more QueryBuilder objects.
#' @param group_names A vector of strings representing tissue groups.
#' @param gene The name of the gene to match with the exon and gene
#'     query results.
#' @param sample_names A vector of strings representing sample names.
#' @param js_params An object returned from call to `get_JunctionSeq_params`.
#'
#' @return
#' \code{get_JunctionSeq_params} returns an object containing the necessary
#'   parameters needed by calls to \code{JunctionSeq::readJunctionSeqCounts} and
#'   \code{JunctionSeq::plotJunctionSeqResultsForGene}.
#'
#' \code{output_example_function_calls} outputs the necessary JunctionSeq
#'   calls needed to produce a plot.
#'
#' @examples
#' \dontrun{
#' sb1 <- QueryBuilder(compilation = "gtex", regions = "chr7:128393029-128394277")
#' sb1 <- set_row_filters(sb1, contains == 1, coverage_sum >= 1000)
#' sb1 <- set_column_filters(sb1, SMTS == "Brain")
#'
#' sb2 <- set_column_filters(sb1, SMTS == "Pituitary")
#'
#' sb3 <- set_column_filters(sb2, SMTS == "Spleen")
#'
#' js <- get_JunctionSeq_params(
#'     query_builders = list(sb1, sb2, sb3),
#'     gene = "IMPDH1",
#'     group_names = list("Brain", "Pituitary", "Spleen")
#' )
#'
#' output_example_function_calls(js)
#' }
#' @export
get_JunctionSeq_params <- function(query_builders, group_names, gene) {
    js <- FormatForJunctionSeq$new(query_builders, group_names, gene)
    js$write_gff_file()

    js
}

#' @export
#' @rdname get_JunctionSeq_params
output_example_function_calls <- function(js_params) {
    assert_that(is(js_params, "FormatForJunctionSeq"))

    js_params$format_for_junction_seq(env = parent.frame())
}
