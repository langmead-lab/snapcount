#test higher level functions in Snapcount

devtools::install_github("langmead-lab/snapr")
library(snapcount)
library(recount)
setwd("D:/snapr/hitests/")


######Merge testing
sb <- SnaptronQueryBuilder$new()

urls1=list("http://snaptron.cs.jhu.edu/encode1159/snaptron?regions=chr1:1879786-1879786&either=2&rfilter=strand:-",
           "http://snaptron.cs.jhu.edu/rpc/snaptron?regions=chr1:1879786-1879786&either=2&rfilter=strand:-")

len<-length(urls1)
sbs1<-lapply(urls1, function(g) { sb$from_url(g)$clone(deep=TRUE) })
             
test_junction_union_output<-junction_union(sbs1[[1]],sbs1[[2]])
save(test_junction_union_output,file="test_junction_union_output.rda")

#write out the counts
write.table(as.matrix(assays(test_junction_union_output)$counts),file="test_junction_union_output.tsv",sep="\t",row.names=FALSE,quote=FALSE)
#now write out the metadata
write.table(colData(test_junction_union_output),file="test_junction_union_output.md.tsv",sep="\t",row.names=FALSE,quote=FALSE)
#write out the coordinates (rows), includes extra snaptron row level metadata
write.table(rowRanges(test_junction_union_output),file="test_junction_union_output.coords.tsv",sep="\t",row.names=FALSE,quote=FALSE)


####SSC testing
#Shared Sample Count (SSC) high level function fails
sb <- SnaptronQueryBuilder$new()
urls1=list("http://snaptron.cs.jhu.edu/gtex/snaptron?regions=chr1:1879786-1879786&either=2&rfilter=strand:-",
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=chr1:9664595-9664595&either=2&rfilter=strand:+",
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=chr6:32831148-32831148&either=2&rfilter=strand:-")
urls2=c("http://snaptron.cs.jhu.edu/gtex/snaptron?regions=chr1:1879903-1879903&either=1&rfilter=strand:-",
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=chr1:9664759-9664759&either=1&rfilter=strand:+",
        "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=chr6:32831182-32831182&either=1&rfilter=strand:-")

len<-length(urls1)
sbs1<-lapply(urls1, function(g) { sb$from_url(g)$clone(deep=TRUE) })
sbs2<-lapply(urls2, function(g) { sb$from_url(g)$clone(deep=TRUE) })
ssc_inputs<-lapply(1:length(sbs1), function(g) { list(sbs1[[g]], sbs2[[g]])})
test_ssc_output <- shared_sample_counts(ssc_inputs[[1]], ssc_inputs[[2]], ssc_inputs[[3]])

save(test_ssc_output,file="test_ssc_output.rda")


########JIR testing
sb1 <- SnaptronQueryBuilder$new()
sb1<-sb1$from_url("http://snaptron.cs.jhu.edu/srav2/snaptron?regions=chr2:29446395-30142858&contains=1&rfilter=strand:-")
sb2 <- SnaptronQueryBuilder$new()
sb2<-sb2$from_url("http://snaptron.cs.jhu.edu/srav2/snaptron?regions=chr2:29416789-29446394&contains=1&rfilter=strand:-")
test_jir_output<-junction_inclusion_ratio(list(sb1),list(sb2))
save(test_jir_output,file="test_jir_output.rda")

write.table(test_jir_output,file="test_jir_output.tsv",sep="\t",row.names=FALSE,quote=FALSE)

######PSI testing
inclusion_group1 <- SnaptronQueryBuilder$new()
inclusion_group1 <- inclusion_group1$from_url("http://snaptron.cs.jhu.edu/srav2/snaptron?regions=chr1:94468008-94472172&exact=1&rfilter=strand:+")
inclusion_group2 <- SnaptronQueryBuilder$new()
inclusion_group2 <- inclusion_group2$from_url("http://snaptron.cs.jhu.edu/srav2/snaptron?regions=chr1:94472243-94475142&exact=1&rfilter=strand:+")
exclusion_group <- SnaptronQueryBuilder$new()
exclusion_group <- exclusion_group$from_url("http://snaptron.cs.jhu.edu/srav2/snaptron?regions=chr1:94468008-94475142&exact=1&rfilter=strand:+")

test_psi_output<-percent_spliced_in(list(inclusion_group1), list(inclusion_group2), list(exclusion_group), min_count=1)
save(test_psi_output,file="test_psi_output.rda")

write.table(test_psi_output,file="test_psi_output.tsv",sep="\t",row.names=FALSE,quote=FALSE)


######TS testing
inclusion_group1 <- SnaptronQueryBuilder$new()
inclusion_group1<-inclusion_group1$from_url("http://snaptron.cs.jhu.edu/gtex/snaptron?regions=chr4:20763023-20763023&either=2&rfilter=strand:-")
inclusion_group2 <- SnaptronQueryBuilder$new()
inclusion_group2<-inclusion_group2$from_url("http://snaptron.cs.jhu.edu/gtex/snaptron?regions=chr4:20763098-20763098&either=1&rfilter=strand:-")

test_ts_output<-tissue_specificity(list(inclusion_group1, inclusion_group2))
save(test_ts_output,file="test_ts_output.rda")

write.table(test_ts_output,file="test_ts_output.tsv",sep="\t",row.names=FALSE,quote=FALSE)
