extract_samples <- function(query_data) {
    unlist(lapply(strsplit(query_data$samples, ","), `[`, -1))
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    assert_that(is.numeric(x))
    all(abs(x - round(x)) < tol)
}

## renaming second arg from values to value to
## appease R checks
`%<-%` <- function(bindings, value) {
    values <- force(value)
    bindings <- substitute(bindings)

    assert_that(length(bindings) == length(values) + 1)

    for (i in seq_along(values)) {
        var_name <- deparse(bindings[[i + 1]])
        assign(var_name, values[[i]], pos = parent.frame())
    }
}

bool_expressions_to_strings <- function(exprs) {
    if (rlang::is_bare_character(exprs)) {
        return(exprs)
    }

    res <- NULL
    if (rlang::is_call(exprs)) {
        if (exprs[[1]] == rlang::as_name("c") ||
            exprs[[1]] == rlang::as_name("list")) {
            res <- lapply(exprs[-1], expression_to_string_helper)
        } else {
            res <- expression_to_string_helper(exprs)
        }
    }

    res
}

expression_to_string_helper <- function(expr) {
    assert_that(
        length(expr) == 3,
        msg = paste(deparse(expr), ": is not a valid filter")
    )
    assert_that(
        is_logical_op(expr[[1]]),
        msg = paste(deparse(expr), ": is not a bool expression")
    )
    assert_that(
        rlang::is_symbol(expr[[2]]),
        msg = paste(deparse(expr[[2]]), ": is not a valid name")
    )

    if (rlang::is_call(expr[[3]]) || rlang::is_symbol(expr[[3]])) {
        res <- rlang::eval_tidy(expr[[3]])
        assert_that(
            rlang::is_syntactic_literal(res),
            msg = paste(
                deparse(expr[[3]]),
                ": does not evaluate to a basic type"
            )
        )
        expr[[3]] <- res
    } else if (rlang::is_bare_character(expr[[3]])) {
        ## convert a string a symbol to prevent deparse from
        ## escaping quotes
        expr[[3]] <- rlang::sym(expr[[3]])
    }

    deparse(expr, backtick = FALSE)
}
