#!/usr/bin/env bash
set -o pipefail -o nounset -o errexit

cut -f 13 $1 | perl -ne 'BEGIN { open(IN,"</data/snaptron_data/'${2}'/samples.tsv"); while($line=<IN>) { next if($line=~/^rail_id/); chomp($line); @f=split(/\t/,$line); $sid=shift(@f); push(@a,$sid); } close(IN); } chomp; $f=$_; @f=split(/,/,$f); shift(@f); %h=(); map { ($a,$b)=split(/:/,$_); $h{$a}=$b; } @f; $first=1; for $e (@a) { print "\t"  if(!$first); $c=$h{$e}; print "$c" if($c); print "0" if(!$c); $first=undef; } print "\n";' > ${1}.counts
diff ${1}.counts ${3}
