`%>%` <- magrittr::`%>%`

Registry <- NULL

.onLoad <- function(...) {
    registry <- httr::GET("http://snaptron.cs.jhu.edu/snaptron/registry")
    start <- 8
    end <- length(registry$content)
    Registry <<- registry$content[start:end] %>%
        rawToChar() %>% jsonlite::fromJSON()
    variants <- names(Registry) %>% as.list()

    Compilation <<- do.call(enum, variants)
}
