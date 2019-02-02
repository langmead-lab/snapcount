num_cores <- parallel::detectCores()
`%>%` <- magrittr::`%>%`

#' @export
junction_inclusion_ratio <- function(group1, group2, group_names = NULL) {
    stopifnot(is.list(group1), is.list(group2))

    c(s1, s2) %<-% parallel::mclapply(list(group1, group2), process_queries, mc.cores = min(num_cores, 2))

    s1 <- s1[, .(coverage = sum(coverage)), by = .(sample_id)]
    s2 <- s2[, .(coverage = sum(coverage)), by = .(sample_id)]

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

process_queries <- function(group) {
    dfs <- lapply(group, function(e) {
        q <- e$query_jx()

        data.table::data.table(sample = extract_samples(q)) %>%
            tidyr::separate("sample", into = c("sample_id", "coverage"), convert = TRUE)
    })

    do.call(rbind, dfs)
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
