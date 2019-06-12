# https://stackoverflow.com/a/44152358
enum <- function(..., class) {
    values <- sapply(match.call(expand.dots = TRUE)[-1L], deparse)

    stopifnot(identical(unique(values), values))

    res <- setNames(seq_along(values), values)
    res <- as.environment(as.list(res))

    class(res) <- c("enum", "environment")
    lockEnvironment(res, bindings = TRUE)

    res
}

#' @export
Coordinates <- enum(Exact, Within, StartIsExactOrWithin, EndIsExactOrWithin, Overlaps)
