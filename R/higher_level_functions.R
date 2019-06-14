#' @import data.table

`%>%` <- magrittr::`%>%`

#' @export
junction_union <- function(...) {
    merge_compilations(..., all = TRUE)
}

#' @export
junction_intersection <- function(...) {
    merge_compilations(..., all = FALSE)
}

merge_compilations <- function(..., all) {
    compilations <- lapply(list(...), function(sb) {
        sb$query_jx(return_rse = FALSE)
    })

    if (length(compilations) == 1) {
        return(compilations[[1]])
    }

    res <- compilations[[1]]
    for (i in 2:length(compilations)) {
        res <- merge(res, compilations[[i]], all = all,
                     by = c("chromosome", "start", "end", "strand")) %>%
            finalize_merge(col_names = names(compilations[[1]]))
    }

    res
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
    args <- list(...)
    strings <- lapply(args, stringr::str_replace_na, replacement = "")

    paste(strings[[1]], strings[[2]], sep = sep)
}

choose_non_na <- function(x, y) {
    ifelse(is.na(x), y, x)
}

calculate_coverage_median <- function(samples) {
    samples[c(FALSE, TRUE)] %>% as.numeric() %>% median()
}

#' @export
junction_inclusion_ratio <- function(group1, group2, group_names = NULL) {
    stopifnot(is.list(group1), is.list(group2))

    c(s1, s2) %<-% run_queries(group1, group2)

    jir <- merge(s1, s2, by = "sample_id", all = TRUE) %>%
        replace_na(0)

    jir[, jir := calc_jir(coverage.y, coverage.x)]

    if (is.null(group_names)) {
        group_names <- c("Group1 Coverage", "Group2 Coverage")
    } else {
        group_names <- paste(group_names, "Coverage")
    }

    data.table::setnames(jir, old = "jir", new = "JIR")
    data.table::setnames(jir, old = "coverage.x", new = group_names[1])
    data.table::setnames(jir, old = "coverage.y", new = group_names[2])

    jir[order(-JIR)]
}

calc_jir <- function(a, b) {
    (a - b)/(a + b + 1)
}

#' @export
percent_spliced_in <- function(inclusion_group1, inclusion_group2, exclusion_group, min_count = 20, group_names = NULL) {
    c(g1, g2, ex) %<-% run_queries(inclusion_group1, inclusion_group2, exclusion_group)

    psi <- merge(g1, g2, by = "sample_id", all = TRUE) %>%
        merge(ex, by = "sample_id", all = TRUE) %>%
        replace_na(0L)

    psi[, psi := calc_psi(coverage.x, coverage.y, coverage, min_count)][]
}

calc_psi <- function(inclusion1, inclusion2, exclusion, min_count) {
    mean_inclusion <- (inclusion1 + inclusion2) / 2.
    total <- mean_inclusion + exclusion
    psi <- mean_inclusion / total

    ifelse(inclusion1 == 0 | inclusion2 == 0 | total < min_count, -1, psi)
}

#' @export
tissue_specificity <- function(..., group_names = NULL) {
    list_of_groups <- list(...)
    num_groups <- length(list_of_groups)

    if (is.null(group_names)) {
        group_names <- paste0("g", 1:num_groups)
    }

    dfs <- lapply(1:num_groups, function(i) {
        g <- list_of_groups[[i]]
        name <- group_names[[i]]
        tissue_specificity_per_group(g[[1]], g[[2]], name)
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
    c(res1, res2) %<-% run_queries(group1, group2, summarize = FALSE)
    res <- merge(res1, res2, by = "sample_id", all = TRUE) %>%
        replace_na(0)
    res <- res[, shared := shared(coverage.x, coverage.y)][, !c("coverage.x", "coverage.y")]

    metadata <- get_compilation_metadata(group1[[1]]$compilation())
    metadata <- metadata[, .(rail_id, SMTS)]

    res <- merge(res, metadata, by.x = "sample_id", by.y = "rail_id", all = TRUE)
    res <- replace_na(res, 0, colnames = setdiff(names(res), "SMTS"))

    res$group <- rep(group_name, nrow(res))
    data.table::setnames(res, old = "SMTS", new = "tissue")

    unique(res)
}

shared <- function(cov1, cov2) {
    as.integer(cov1 != 0 & cov2 != 0)
}

#' @export
shared_sample_counts <- function(..., group_names = NULL) {
    list_of_groups <- list(...)
    num_groups <- length(list_of_groups)

    counts <- lapply(list_of_groups, function(g) {
        shared_sample_count(g[[1]], g[[2]])
    })

    if (is.null(group_names)) {
        group_names <- paste0("g", 1:num_groups)
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

    c(g1, g2) %<-% run_queries(group1, group2)
    intersect(g1$sample_id, g2$sample_id) %>% length()
}

run_queries <- function(..., summarize = TRUE) {
    lapply(list(...), count_samples, summarize = summarize)
}

count_samples <- function(group, summarize = TRUE) {
    dfs <- lapply(group, function(e) {
        q <- e$query_jx(return_rse = FALSE)

        data.table::data.table(sample = extract_samples(q)) %>%
            tidyr::separate("sample", into = c("sample_id", "coverage"), convert = TRUE)
    })

    res <- do.call(rbind, dfs)

    if (summarize) {
        res[, .(coverage = sum(coverage)), by = .(sample_id)]
    } else {
        res
    }
}

`%<-%` <- function(bindings, values) {
    values <- force(values)
    bindings <- substitute(bindings)

    stopifnot(length(bindings) == length(values) + 1)

    for (i in 1:length(values)) {
        var_name <- deparse(bindings[[i+1]])
        assign(var_name, values[[i]], pos = parent.frame())
    }
}

is_query_builder <- function(object) {
    return("SnaptronQueryBuilder" %in% class(object))
}

replace_na <- function(dt, replacement, colnames = NULL) {
    if (!is.null(colnames)) {
        for (name in colnames) {
            set(dt, which(is.na(dt[[name]])), name, 0)
        }
    } else {
        for (i in seq_along(dt)) {
            set(dt, which(is.na(dt[[i]])), i, 0)
        }
    }

    dt
}
