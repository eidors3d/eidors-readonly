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
   (cd $REPO && git checkout $GITVER && git clean -d -x -f);
   DATE=`(cd $REPO && git show -s --format=%ct $GITVER) | perl -pe'chomp'`;
   for PTH in $REPO $REPO/htdocs $REPO/dev ; do
       find $PTH -type f -name \*.m -print0 2>/dev/null | xargs -0 perl \
           -e  'BEGIN{my $lines = 0};' \
           -ne 'next if /^\s*%/; next if /^\s*$/; $lines++;' \
           -e  'END{printf"  %d,",$lines;}' >> $OUTF ; 
   done
   echo "  $ver,    $DATE; " >> $OUTF;
done
echo "];" >> $OUTF
