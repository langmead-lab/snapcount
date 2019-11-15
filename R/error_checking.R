validate_sample_filter <- function(compilation, name, value) {
    if (is.null(registry <- pkg_globals$registry)) {
        return()
    }
    if (is.null(expected_type <- registry[[compilation]][[name]])) {
        field_names <- names(registry[[compilation]]) %>% sort()
        msg <- paste0(
            "Error: ", "`", name, "'",
            " is not a valid sample filter."
        )
        sorted_indexes <-
            utils::adist(field_names, name) %>% as.vector() %>% order()
        closest_match <- field_names[sorted_indexes][1]
        if (!identical(match, character(0))) {
            msg <-
                paste0(msg, " Perhaps you meant: ", closest_match, "?")
        }
    } else {
        type <- value_to_snaptron_type(value)
        if (expected_type == type ||
            (expected_type == "f" && type == "i")) {
            return(NULL)
        }
        else {
            msg <- paste0(
                "Error: ", "`", name, "'",
                " filter expects value of type ", make_verbose(expected_type),
                ", but got ", make_verbose(type)
            )
        }
    }

    msg
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
    switch(
        type,
        f = "Float",
        t = "String",
        i = "Integer"
    )
}

validate_range_filters <- function() {

}
