num_cores <- parallel::detectCores()
`%>%` <- magrittr::`%>%`

#' @export
junction_inclusion_ratio <- function(group1, group2, group_names = NULL) {
    stopifnot(is.list(group1), is.list(group2))

    c(s1, s2) %<-% run_queries(group1, group2)

    jir <- s1[s2, on = .(sample_id), all = TRUE]
    jir[is.na(jir)] <- 0
    jir[, jir := calc_jir(coverage, i.coverage)]

    if (is.null(group_names)) {
        group_names <- c("Group1 Coverage", "Group2 Coverage")
    } else {
        group_names <- paste(group_names, "Coverage")
    }

    data.table::setnames(jir, old = "jir", new = "JIR")
    data.table::setnames(jir, old = "coverage", new = group_names[1])
    data.table::setnames(jir, old = "i.coverage", new = group_names[2])

    jir[order(-JIR)]
}

calc_jir <- function(a, b) {
    (a - b)/(a + b + 1)
}

#' @export
percent_spliced_in <- function(inclusion_group1, inclusion_group2, exclusion_group, group_names = NULL) {
    c(g1, g2, ex) %<-% run_queries(inclusion_group1, inclusion_group2, exclusion_group)

    psi <- merge(g1, g2, by = "sample_id", all = TRUE) %>%
        merge(ex, by = "sample_id", all = TRUE)
    psi[is.na(psi)] <- 0

    psi[, psi := calc_psi(coverage.x, coverage.y, coverage)]
}

calc_psi <- function(inclusion1, inclusion2, exclusion) {
    mean_inclusion = (inclusion1 + inclusion2) / 2.
    total = mean_inclusion + exclusion

    if (inclusion1 == 0 || inclusion2 == 0) {
        return(-1)
    }
    mean_inclusion / total
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
    res <- rbind(res1, res2)
    res <- res[, .N, by = sample_id]
    res$N <- ifelse(res$N > 1, 1, 0)

    metadata <- get_compilation_metadata(group1[[1]]$compilation())
    metadata <- metadata[, .(rail_id, SMTS)]

    res <- merge(res, metadata, by.x = "sample_id", by.y = "rail_id", all = TRUE)
    res[is.na(res)] <- 0
    res$group <- rep(group_name, nrow(res))
    data.table::setnames(res, old = "SMTS", new = "tissue")
    data.table::setnames(res, old = "N", new = "shared")

    res
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
    groups <- list(...)
    suggested_cores <- min(num_cores, length(groups))
    parallel::mclapply(groups, count_samples, summarize = summarize, mc.cores = suggested_cores)
}

count_samples <- function(group, summarize = TRUE) {
    dfs <- lapply(group, function(e) {
        q <- e$query_jx()

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
