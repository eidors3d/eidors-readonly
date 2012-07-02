#!/bin/sh
find ../../htdocs/tutorial/ -name *.m | sort | \
perl -ne'$tuts{$1}++ if /(.*)[0-9][0-9].*\.m/;' \
   -e'END{foreach (keys %tuts){' \
   -e' $files = `ls "$_*"`;' \
   -e' print "$_\n";'  \
   -e'};}'

 # `ls $_*.m`;'
