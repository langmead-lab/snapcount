SnaptronQueryBuilder <- R6Class("SnaptronQueryBuilder",
    public = list(
        initialize = function(...) {
            private$query <- list(...)
        },
        compilation = function(compilation = NULL) {
            if (!missing(compilation)) {
                assert_that(compilation %in% names(Compilation),
                            msg = paste0(compilation, ": is not a valid compilation"))
                private$query$compilation <- compilation
                invisible(self)
            } else {
                private$query$compilation
            }
        },
        regions = function(regions = NULL) {
            if (!missing(regions)) {
                if (is(regions, "GRanges")) {
                    private$query$regions <- regions
                } else if (is_genes_or_intervals(regions)) {
                    private$query$regions <- regions
                } else {
                    message <-
                        paste("regions must contain strings representing",
                              "HUGO genes or chromosome intervals, of the form",
                              "chr:start-end or chr:start-end:strand,",
                              "or GRanges object.")
                    private$query$regions <-
                        tryCatch(GenomicRanges::GRanges(regions),
                                 error = function(e) {
                                     stop(message)
                                 })
                }
                invisible(self)
            } else {
                private$query$regions
            }
        },
        row_filters = function(...) {
            if (!missing(...)) {
                private$query$row_filters <-
                    bool_expressions_to_strings(rlang::enquos(...))
                invisible(self)
            } else {
                private$query$row_filters
            }
        },
        column_filters = function(...) {
            if (!missing(...)) {
                private$query$column_filters <-
                    bool_expressions_to_strings(rlang::enquos(...))
                invisible(self)
            } else {
                private$query$column_filters
            }
        },
        sids = function(sids = NULL) {
            if (!missing(sids)) {
                assert_that(is.wholenumber(sids),
                    msg = "sids should be whole numbers")
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
            private$call("query_jx", list(sb = self, return_rse = return_rse))
        },
        query_exon = function(return_rse = TRUE) {
            private$call("query_exon", list(sb = self, return_rse = return_rse))
        },
        query_gene = function(return_rse = TRUE) {
            private$call("query_gene", list(sb = self, return_rse = return_rse))
        },
        from_url = function(url) {
            url <- httr::parse_url(url)
            if (url$hostname != "snaptron.cs.jhu.edu") {
                stop("URL does not point to Snaptron server", stop. = FALSE)
            }
            resp <- httr::HEAD(url)
            if (resp$status_code != 200 ||
                httr::http_type(resp) != "text/plain") {
                stop(sprintf("%s: is not a valid URL", url), call. = FALSE)
            }
            query <- list()
            for (i in seq_along(url$query)) {
                name <- switch(n <- names(url$query[i]),
                    rfilter = "row_filters",
                    sfilter = "column_filters",
                    regions = "regions",
                    n)

                if (name == "sids") {
                    query[[name]] <-
                        scan(textConnection(url$query[[i]]), sep = ",")
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

            query$compilation <- strsplit(url$path, "/")[[1]][1]
            private$query <- query

            invisible(self)
        },
        print = function() {
            cat("<SnaptronQueryBuilder>\n")
            for (param in names(private$query)) {
                if (is.null(private$query[[param]])) {
                    next
                } else if (param == "coordinate_modifier") {
                    desc <- switch(private$query[[param]],
                                   Exact = "exact",
                                   Within = "contains",
                                   StartIsExactOrWithin = "either=1",
                                   EndIsExactOrWithin = "either=2",
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
            if (is.null(self$compilation())) {
                stop(
                    paste("Please set a compilation before running", fn_name),
                    call. = FALSE
                )
            }
            if (is.null(self$regions())) {
                stop(
                    paste("Please specify query regions before running", fn_name),
                    call. = FALSE
                )
            }
            fn <- get(fn_name, parent.frame())
            do.call(fn, args)
        }
    )
)
