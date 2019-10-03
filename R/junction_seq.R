#' @export
FormatForJunctionSeq <-
    R6::R6Class(
        "FormatForJunctionSeq",
        public = list(
            initialize = function(query_builders, group_names, gene) {
                private$query_builders <- query_builders
                private$group_names <- group_names
                private$gene <- gene
                private$design <-
                    data.frame(condition = factor(unlist(group_names)))
                private$create_gff_inputs()
            },

            ## format_for_junction_seq = function() {
            ##     f1 <- paste("js_counts <- readJunctionSeqCounts(countdata = get_count_data(), samplenames = get_sample_names(), design = get_design(), flat.gff_file = get_gff_filename())", sep = ", ")
            ##     f2 <- paste("plotJunctionSeqResultsForGene(geneID = )", private)
            ## },

            get_gff_filename = function() {
                if (is.null(private$gff_filename)) {
                    stop(
                        "GFF flat file has not been saved to disk, please run 'write_gff_file' with an optional filename"
                    )
                }
                private$gff_filename
            },

            get_gff_data = function() {
                private$gff
            },

            write_gff_file = function(filename = NULL) {
                if (is.null(filename)) {
                    private$gff_filename <- paste0(private$gene, ".gff")
                } else {
                    private$gff_filename <- filename
                }
                writeLines(private$gff, private$gff_filename)
            },

            get_count_data = function() {
                private$count_dataframes
            },

            get_design = function() {
                private$design
            }
        ),
        private = list(
            idx = 0,
            design = NULL,

            query_builders = list(),
            group_names = list(),
            gene = NULL,

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
                    private$exons[[name]] <-
                        private$query_builders[[i]]$clone()$range_filters(NULL)$query_exon(return_rse = FALSE)
                    private$genes[[name]] <-
                        private$query_builders[[i]]$clone()$range_filters(NULL)$query_gene(return_rse = FALSE)
                })
            },

            process_junctions = function() {
                for (group in private$group_names) {
                    apply(private$junctions[[group]], 1, function(row) {
                        key <-
                            private$create_key(row[["chromosome"]], row[["start"]],
                                row[["end"]], row[["strand"]])
                        if (key %in% names(private$key_to_group)) {
                            idx <- private$key_to_index[[key]]
                            type <- private$key_to_type[[key]]
                        } else {
                            private$idx <- private$idx + 1
                            private$key_to_index[[key]] <-
                                private$idx
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

                        private$write_count(group,
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
                            base::strsplit(row[["gene_id:gene_name:gene_type:bp_length"]],
                                split = ":")[[1]][2]
                        if (gene_name != private$gene) {
                            return()

                        }

                        key <-
                            private$create_key(row[["chromosome"]], row[["start"]],
                                row[["end"]], row[["strand"]])
                        if (key %in% names(private$key_to_group)) {
                            idx <- private$key_to_index[[key]]
                            type <- private$key_to_type[[key]]
                        } else {
                            private$idx <- private$idx + 1
                            private$key_to_index[[key]] <-
                                private$idx
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

                        private$write_count(group,
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
                            base::strsplit(row[["gene_id:gene_name:gene_type:bp_length"]],
                                split = ":")[[1]][2]
                        if (gene_name != private$gene) {
                            return()
                        }

                        key <-
                            private$create_key(row[["chromosome"]], row[["start"]],
                                row[["end"]], row[["strand"]])
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

                        private$write_count(group,
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
                if (is.null(private$key_to_group)
                    ||
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
                            data.table::fread(
                                tsv,
                                sep = "\t",
                                header = FALSE
                             )
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
            }
        )
    )
