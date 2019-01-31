#!/usr/bin/env Rscript
library(SummarizedExperiment)
library(devtools)
#assume we're being run from within external_tests
path='../'
load_all(path)
cd99.jx<-query_jx(compilation='srav2', genes_or_intervals='CD99', range_filters='annotated:1', sample_filters='description:cortex')
cd99.gene<-query_gene(compilation='srav2', genes_or_intervals='CD99', sample_filters='description:cortex')
cd99.gene2<-query_gene(compilation='gtex', genes_or_intervals='CD99', sample_filters='SMTS:Brain')
cd99.exon<-query_exon(compilation='gtex', genes_or_intervals='CD99', sample_filters='SMTS:Brain')
cd99.base<-query_coverage(compilation='gtex', genes_or_intervals='chrX:2691179-2691188')
cd99.base_sids<-query_coverage(compilation='gtex', genes_or_intervals='chrX:2691179-2691188', sids=c(50099,50200,59759))
gtex.metadata<-get_compilation_metadata('gtex')
for(z in c('cd99.jx','cd99.exon','cd99.gene','cd99.gene2','cd99.base','cd99.base_sids'))
{
    o <- get(z)
    write.table(rowRanges(o),file=paste0(z,'.row.r.tsv'),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(colData(o),file=paste0(z,'.col.r.tsv'),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(as.matrix(assays(o)$counts),file=paste0(z,'.count.r.tsv'),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
}

for(z in c('cd99.jx','cd99.exon','cd99.gene','cd99.gene2','cd99.base','cd99.base_sids'))
{
    o <- get(z)
    save(o,file=paste0(z,".rda"))
}
