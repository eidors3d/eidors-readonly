#!/bin/bash
# Sequentially checkout versions 
REPO=${1:-"$HOME/docs/eidors-clean"}
echo "Look at $REPO";
OUTF=datafile.m
LASTREV=$(cd $REPO && svn info --show-item revision)
echo "LastRev is $LASTREV";
echo 'loc = [ % all eidors  htdocs   dev   version  date' > $OUTF
for ver in `seq 1 $LASTREV`; do 
   echo "$REPO: VER=$ver";
   (cd $REPO && svn up -r$ver )
   DATE=`(cd $REPO && svn info --show-item last-changed-date) | perl -MDate::Parse -pe'$_=str2time($_)'`;
   for PTH in $REPO $REPO/eidors $REPO/htdocs $REPO/dev ; do
       find $PTH -type f -name \*.m -print0 2>/dev/null | xargs -0 perl \
           -e  'BEGIN{my $lines = 0};' \
           -ne 'next if /^\s*%/; next if /^\s*$/; $lines++;' \
           -e  'END{printf"  %d,",$lines;}' >> $OUTF ; 
   done
   echo "  $ver,    $DATE; " >> $OUTF;
done
echo "];" >> $OUTF
