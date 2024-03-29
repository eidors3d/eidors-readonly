#!/bin/bash
# Sequentially checkout versions 
REPO=${1:-"$HOME/docs/eidors-git"}
echo "Look at $REPO";
OUTF=datafile.m
LASTREV=$(cd $REPO && git svn log --oneline HEAD^^^^^..HEAD | head -1 | cut -f1 -d'|' | sed 's/r//')
echo "LastRev is $LASTREV";
echo 'loc = [ % all(/)  htdocs   dev   version  date' > $OUTF
for ver in `seq 1 $LASTREV`; do 
   echo "$REPO: VER=$ver";
   GITVER=`(cd $REPO && git svn find-rev r$ver)`;
   (cd $REPO && git checkout $GITVER);
   DATE=`(cd $REPO && git show -s --format=%ct $GITVER) | perl -pe'chomp'`;
   find $REPO -type f -name \*.m -print0 | xargs -0 perl \
       -e  'BEGIN{my $lines = 0};' \
       -ne 'next if /^\s*%/; next if /^\s*$/; $lines++;' \
       -e  'END{printf"  %d,",$lines;}' >> $OUTF ; # all (/)
   find $REPO/htdocs -path $REPO/htdocs/doc -prune -o -type f -name \*.m -print0 | xargs -0 perl \
       -e  'BEGIN{my $lines = 0};' \
       -ne 'next if /^\s*%/; next if /^\s*$/; $lines++;' \
       -e  'END{printf"  %d,",$lines;}' >> $OUTF ; # htdocs
   find $REPO/dev -type f -name \*.m -print0 | xargs -0 perl \
       -e  'BEGIN{my $lines = 0};' \
       -ne 'next if /^\s*%/; next if /^\s*$/; $lines++;' \
       -e  'END{printf"  %d,",$lines;}' >> $OUTF ; # dev
   echo "  $ver,    $DATE; " >> $OUTF;
done
echo "];" >> $OUTF
