`%>%` <- magrittr::`%>%`

pkg_globals <- new.env(parent = emptyenv())

.onLoad <- function(...) {
    json_data <- httr::GET("http://snaptron.cs.jhu.edu/snaptron/registry")
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
