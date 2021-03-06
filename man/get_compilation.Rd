% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query_builder_wrappers.R
\name{get_compilation}
\alias{get_compilation}
\alias{set_compilation}
\title{Get or set query compilation}
\usage{
get_compilation(qb)

set_compilation(qb, compilation)
}
\arguments{
\item{qb}{A QueryBuilder object constructed using the
\code{\link{QueryBuilder}} function.}

\item{compilation}{A single string containing the name of the Snaptron data
source. Any variant of the \code{Compilation} enum can also be passed as an
argument.}
}
\value{
\code{get_compilation} returns the current compilation as string.
\code{set_compilation} returns a new \code{QueryBuilder} object with
the compilation set to the value of \code{compilation}.
}
\description{
Get or set query compilation
}
\examples{
qb <- QueryBuilder(compilation = "gtex", regions = "CD99")
get_compilation(qb)
qb <- set_compilation(qb, Compilation$tcga)
get_compilation(qb)
}
