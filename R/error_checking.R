`%>%` <- magrittr::`%>%`

sample_filter_validator <- function() {
    registry <- httr::GET("http://snaptron.cs.jhu.edu/snaptron/registry")
    start <- 8
    end <- length(registry$content)
    table <- registry$content[start:end] %>%
        rawToChar() %>% jsonlite::fromJSON()

    function(compilation, name, value){
        msg <- NULL
        if (is.null(type <- table[[compilation]][[name]])) {
            field_names <- names(table[[compilation]])
            msg <- paste0("Error: ", "`", name, "'",
                          " is not a valid sample filter.")
            match <- agrep(name, field_names,
                           ignore.case = TRUE, value = TRUE, max.distance = 0.1)

            if (!identical(match, character(0))) {
                msg <- paste0(msg, " Perhaps you meant: ",
                              match[1], "?")
            }
        } else {
            if (value_to_snaptron_type(value) != type) {
                msg <- paste0("Error: ", "`", name, "'",
                              " filter expects value of type ",
                              make_verbose(type), ", but got ",
                              make_verbose(value_to_snaptron_type(value)))
            }
        }

        msg
    }
}

value_to_snaptron_type <- function(value) {
    n <- suppressWarnings(as.numeric(value))

    if (is.na(n)) {
        return("t")
    }

    if (is.wholenumber(n)) {
        return("i")
    }

    return("f")
}

make_verbose <- function(type) {
    switch(type,
           f = "Float",
           t = "String",
           i = "Integer")
}

validate_range_filters <- function() {

}
