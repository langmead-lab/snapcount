context("FormatForJunctionSeq")

quote_string <- function(string) {
    stringr::str_replace_all(string, "(\\W)", "\\\\\\1")
}

test_that("test FormatForJunctionSeq class", {
    sb1 <- SnaptronQueryBuilder$new()
    sb1$compilation("gtex")
    sb1$regions("chr7:128393029-128394277")
    sb1$row_filters(contains == 1, coverage_sum >= 1000)
    sb1$column_filters(SMTS == "Brain")

    sb2 <- sb1$clone(deep = TRUE)
    sb2$column_filters(SMTS == "Pituitary")

    sb3 <- sb2$clone(deep = TRUE)
    sb3$column_filters(SMTS == "Brain")

    expect_silent(
        js <- get_JunctionSeq_params(query_builders = list(sb1, sb2, sb3),
                                     group_names = list("Brain", "Pituitary", "Spleen"),
                                     gene = "IMPDH1"))

    expect_output(output_example_function_calls(js),
                  'jscs <- JunctionSeq::readJunctionSeqCounts\\(countdata=js\\$countdata, samplenames=js\\$samplenames, design=js\\$design, flat\\.gff\\.file=js\\$flat.gff.file\\)\\n\\nJunctionSeq::plotJunctionSeqResultsForGene\\(geneID=js\\$geneID, jscs=jscs, plot\\.type="rawCounts", displayTranscripts=TRUE, plot\\.exon\\.results=TRUE, plot\\.junction\\.results=TRUE, plot\\.novel\\.junction.results=TRUE, plot\\.untestable\\.results=TRUE\\)')

    # clean up side-effect of IMPDH1.gff file being written in the working
    # directory, triggers a check warning otherwise
    unlink(js$flat.gff.file)
})
