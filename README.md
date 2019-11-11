<!-- README.md is generated from README.Rmd. Please edit that file -->



# snapr
[![Build Status](https://travis-ci.com/langmead-lab/snapr.svg?token=vEUBb2QKjox3PdRAssp8&branch=master)](https://travis-ci.com/langmead-lab/snapr)
[![Coverage Status](https://img.shields.io/codecov/c/github/langmead-lab/snapr/master.svg)](https://codecov.io/gh/langmead-lab/snapr?branch=master)

## Overview
`snapr` is a package for interfacing with Snaptron's REST API.

## Installation


```r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("langmead-lab/snapr")
```

## Usage

##### snapr can be used with either a procedural interface

```r
library(snapcount)

query_jx(compilation = "gtex", regions = "CD99")
#> class: RangedSummarizedExperiment 
#> dim: 3485 9662 
#> metadata(0):
#> assays(1): counts
#> rownames(3485): 28340058 28340273 ... 28352407 28352408
#> rowData names(12): DataSource:Type snaptron_id ... coverage_median
#>   source_dataset_id
#> colnames(9662): 50099 50100 ... 59759 59760
#> colData names(322): rail_id Run ... junction_coverage
#>   junction_avg_coverage
query_jx(compilation = "gtex", regions = "CD99", range_filters = samples_count == 10)
#> class: RangedSummarizedExperiment 
#> dim: 25 226 
#> metadata(0):
#> assays(1): counts
#> rownames(25): 28342638 28343346 ... 28352394 28352402
#> rowData names(12): DataSource:Type snaptron_id ... coverage_median
#>   source_dataset_id
#> colnames(226): 50117 50145 ... 59629 59713
#> colData names(322): rail_id Run ... junction_coverage
#>   junction_avg_coverage
```

##### Or using the query-builder class

```r
sb <- SnaptronQueryBuilder$new()
sb$compilation("gtex")$regions("CD99")$query_jx()
#> class: RangedSummarizedExperiment 
#> dim: 3485 9662 
#> metadata(0):
#> assays(1): counts
#> rownames(3485): 28340058 28340273 ... 28352407 28352408
#> rowData names(12): DataSource:Type snaptron_id ... coverage_median
#>   source_dataset_id
#> colnames(9662): 50099 50100 ... 59759 59760
#> colData names(322): rail_id Run ... junction_coverage
#>   junction_avg_coverage

# sb$from_url("http://snaptron.cs.jhu.edu/gtex/bases?regions=BRCA1&sids=50099,50102,50113")$query_gene()
```

## Reference
For more information on Snaptron please visit: http://snaptron.cs.jhu.edu/index.html.
