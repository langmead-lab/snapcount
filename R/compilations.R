foo1 <- function(x) {
    foo2(rlang::enexpr(x))
}

foo2 <- function(x) {
    print(x)
}
