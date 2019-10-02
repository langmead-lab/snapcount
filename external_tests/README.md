
This will serve for now for smoke testing the output of `snapcount` against Snaptron directly.  

It takes ~1m to run (on Stingray at least):

`cd external_tests && /bin/bash -x check.sh`

Then check the `*.diff` files for any that are not 0-sized.

The R script requires the [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) package to be installed (via BioConductor) for its `rowRanges` function.
