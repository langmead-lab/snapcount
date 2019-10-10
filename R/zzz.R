pkg_globals <- new.env(parent = emptyenv())

.onLoad <- function(...) {
    snaptron_host <- "http://snaptron.cs.jhu.edu"

    if (!is.null(host <- getOption("snapcount.host"))) {
        snaptron_host <- host
    }
    if (!is.null(port <- getOption("snapcount.port"))) {
        if (is.wholenumber(port) && port != 80) {
            snaptron_host <- paste(snaptron_host, port, sep = ":")
        }
    }
    snaptron_host <- paste0(snaptron_host, "/")
    assign("snaptron_host", snaptron_host, envir = pkg_globals)
    json_data <-
        paste(pkg_globals$snaptron_host, "snaptron", "registry", sep = "/") %>%
        httr::GET()
    if (httr::http_error(json_data)) {
        return()
    }
    start <- 8
    end <- length(json_data$content)
    registry <- json_data$content[start:end] %>%
        rawToChar() %>% jsonlite::fromJSON()
    compilation_names <- names(registry) %>% as.list()

    if (length(compilation_names) == 0) {
        return()
    }
    assign("registry", registry, pkg_globals)
    assign("metadata", list(), pkg_globals)
    Compilation <<- do.call(enum, compilation_names)
}
