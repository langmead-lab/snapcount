enum <- function(...) {
    variants <- rlang::ensyms(...)
    values <- purrr::map(variants, rlang::as_string)

    res <- setNames(values, values)
    res <- as.environment(as.list(res))

    lockEnvironment(res, bindings = TRUE)

    res
}

#' @export
Coordinates <- enum(Exact, Within, StartIsExactOrWithin, EndIsExactOrWithin, Overlaps)

#' @export
Compilation <- enum(gtex, tcga, srav2, sra)
