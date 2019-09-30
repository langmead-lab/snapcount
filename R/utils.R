extract_samples <- function(query_data) {
    unlist(lapply(strsplit(query_data$samples, ","), `[`, -1))
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    assert_that(is.numeric(x))
    all(abs(x - round(x)) < tol)
}

`%<-%` <- function(bindings, values) {
    values <- force(values)
    bindings <- substitute(bindings)

    assert_that(length(bindings) == length(values) + 1)

    for (i in seq_along(values)) {
        var_name <- deparse(bindings[[i + 1]])
        assign(var_name, values[[i]], pos = parent.frame())
    }
}

create_conjunction <- function(expressions) {
    if (length(expressions) == 0) {
        return(NULL)
    }
    if (length(expressions) == 1) {
        return(expressions[[1]])
    }
    for (i in seq_along(expressions)) {
        if (i == 1) {
            predicate <- expressions[[1]]
        } else {
            predicate <- rlang::expr(!!predicate & !!expressions[[i]])
        }
    }

    predicate
}

string_to_bool_expression <- function(sample_filters) {
    sample_filters <- gsub("\\s*<:\\s*:", "<=", sample_filters)
    sample_filters <- gsub("\\s*>:\\s*:", ">=", sample_filters)
    sample_filters <- gsub("\\s*:\\s*", "==", sample_filters)

    exprs <- lapply(sample_filters, string_to_expression_helper)
    create_conjunction(exprs)
}

string_to_expression_helper <- function(string) {
    expr <- rlang::parse_expr(string)

    if (!rlang::is_syntactic_literal(expr[[3]])) {
        if (rlang::is_call(expr[[3]])) {
            operator <- rlang::as_string(expr[[1]])
            c(lhs, rhs) %<-% stringr::str_split(string, operator, n = 2)[[1]]
            expr[[3]] <- rhs
        } else {
            expr[[3]] <- rlang::as_string(expr[[3]])
        }
    }

    expr
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
    assert_that(length(expr) == 3,
                msg = paste(deparse(expr), ": is not a valid filter"))
    assert_that(is_logical_op(expr[[1]]),
                msg = paste(deparse(expr), ": is not a bool expression"))
    assert_that(rlang::is_symbol(expr[[2]]),
                msg = paste(deparse(expr[[2]]), ": is not a valid name"))

    if (rlang::is_call(expr[[3]]) || rlang::is_symbol(expr[[3]])) {
        res <- rlang::eval_tidy(expr[[3]])
        assert_that(rlang::is_syntactic_literal(res),
                    msg = paste(deparse(expr[[3]]),
                                ": does not evaluate to a basic type"))
        expr[[3]] <- res
    } else if (rlang::is_bare_character(expr[[3]])) {
        ## convert a string a symbol to prevent deparse from
        ## escaping quotes
        expr[[3]] <- rlang::sym(expr[[3]])
    }

    deparse(expr, backtick = FALSE)
}

apply_sample_filters_to_metadata <- function(metadata, sample_filters) {
    if (is.null(sample_filters)) {
        return(metadata)
    }

    predicate_expression <- string_to_bool_expression(sample_filters)
    eval(rlang::expr(metadata[!!predicate_expression]))
}
