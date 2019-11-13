context("FormatForJunctionSeq")


test_that("test FormatForJunctionSeq class", {
    sb1 <- SnaptronQueryBuilder$new()
    sb1$compilation("gtex")
    sb1$regions("chr7:128393029-128394277")
    sb1$range_filters(c(contains == 1, coverage_sum >= 1000))
    sb1$sample_filters(SMTS == "Brain")

    sb2 <- sb1$clone(deep = TRUE)
    sb2$sample_filters(SMTS == "Pituitary")

    sb3 <- sb2$clone(deep = TRUE)
    sb3$sample_filters(SMTS == "Brain")

    expect_silent(js <- FormatForJunctionSeq$new(list(sb1, sb2, sb3),
                                                  group_names = list("Brain", "Pituitary", "Spleen"),
                                                  gene = "IMPDH1"))
    expect_error(js$get_gff_filename(),
                 "please call 'write_gff_file' method first")
    expect_silent(js$write_gff_file())
    expect_output(js$format_for_junction_seq(),
                  'js_counts <- JunctionSeq::readJunctionSeqCounts\\(countdata=js\\$get_count_data\\(\\), samplenames=js\\$get_sample_names\\(\\), design=js\\$get_design\\(\\), flat\\.gff_file=js\\$get_gff_filename\\(\\)\\)\\n\\nJunctionSeq::plotJunctionSeqResultsForGene\\(geneID=js\\$get_gene_name\\(\\), jscs=js_counts, plot\\.type="rawCounts", displayTranscripts=TRUE\\)')

    #clean up side-effect of IMPDH1.gff file being written in the working directory, triggers a check warning otherwise
    file.remove(js$get_gff_filename())
})
