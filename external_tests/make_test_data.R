#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(devtools)
#assume we're being run from within external_tests
path='../'
load_all(path)

sb <- SnaptronQueryBuilder$new()
sb$compilation("srav2")$regions("CD99")$row_filters(annotated == 1)$column_filters(description == "cortex")
cd99.jx <- query_jx(sb)
sb$row_filters(NULL)
cd99.gene<-query_gene(sb)
sb$compilation("gtex")$range_filters(SMTS == "Brain")
cd99.gene2<-query_gene(sb)
cd99.exon<-query_exon(sb)
gtex.metadata<-get_compilation_metadata('gtex')
for(z in c('cd99.jx','cd99.exon','cd99.gene','cd99.gene2'))
{
    o <- get(z)
    write.table(rowRanges(o),file=paste0(z,'.row.r.tsv'),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(colData(o),file=paste0(z,'.col.r.tsv'),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(as.matrix(assays(o)$counts),file=paste0(z,'.count.r.tsv'),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
}

for(z in c('cd99.jx','cd99.exon','cd99.gene','cd99.gene2'))
{
    o <- get(z)
    save(o,file=paste0(z,".rda"))
}
