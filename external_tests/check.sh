#!/usr/bin/env bash
set -o pipefail -o nounset -o errexit
#checks output of snapR/snaptronR package against directly run Snaptron web service queries
#at the end runs diffs, so ls -ltr *.diff

#first get the data from Snaptron directly
cat all.urls | perl -ne 'chomp; ($fn,$u)=split(/\t/,$_); `curl "$u" | tail -n +2 > $fn`;'

#now get it through the R module
Rscript ./make_test_data.R

for f in gene gene2 exon jx; do
    c=`perl -e '$c="'${f}'"; if($c eq "jx" || $c eq "gene") { print "srav2"; } else { print "gtex";}'`
    /bin/bash -x diff_counts.sh ${f}.tsv $c cd99.${f}.count.r.tsv > ${f}.count.diff
    f1=${f}.tsv
    paste <(cut -f 3-7 $f1) <(cut -f 1-2 $f1) <(cut -f 8-12,14- $f1) | perl -ne 'chomp; $f=$_; $f=~s/\t\t\t\t/\tNA\tNA\tNA\t/g; $f=~s/\t\t/\tNA\t/; $f=~s/\.0\t/\t/g; print "$f\n";' > ${f1}.cut
    diff ${f1}.cut cd99.${f}.row.r.tsv > ${f}.row.diff
    cat /data/snaptron_data/${c}/samples.tsv | tail -n +2 | perl -ne 'chomp; $f=$_; $f=~s/\t(\d+)\.0\t/\t$1\t/g; $f=~s/\t(\d+)\.(\d+)0\t/\t$1\.$2\t/g; $f=~s/\t\.(\d+)\t/\t0\.$1\t/g; $f=~s/\ttrue\t/\tTRUE\t/g; $f=~s/\tfalse\t/\tFALSE\t/g; $f=~s/\t0\.0\t/\t0\t/g; $f=~s/\t0\.0$/\t0/g; while($f=~/\t\t/) { $f=~s/\t\t/\tNA\t/; } print "$f\n";' > ${c}.samples 
    cat cd99.${f}.col.r.tsv | perl -ne 'chomp; $f=$_; $f=~s/\t(\d+)\.(\d+)0\t/\t$1\.$2\t/g; while($f=~/\t\t/) { $f=~s/\t\t/\tNA\t/; } print "$f\n";' > cd99.${f}.col.r.tsv.samples
    diff ${c}.samples cd99.${f}.col.r.tsv.samples > ${f}.col.diff
done
c='gtex'
for f in base base_sids; do
    cut -f 5- ${f}.tsv > ${f}.tsv.counts
    diff ${f}.tsv.counts cd99.${f}.count.r.tsv > ${f}.count.diff
    cut -f 2-4 ${f}.tsv | perl -ne 'chomp; print "$_\t2\t*\n";' > ${f}.tsv.cut
    diff ${f}.tsv.cut cd99.${f}.row.r.tsv > ${f}.row.diff
    cat cd99.${f}.col.r.tsv | perl -ne 'chomp; $f=$_; $f=~s/\t(\d+)\.(\d+)0\t/\t$1\.$2\t/g; while($f=~/\t\t/) { $f=~s/\t\t/\tNA\t/; } print "$f\n";' > cd99.${f}.col.r.tsv.samples
    diff ${c}.samples cd99.${f}.col.r.tsv.samples > ${f}.col.diff
done
