#!/bin/bash
# Sequentially checkout versions 
REPO=${1:-"$HOME/docs/eidors-git"}
OUTF=datafile.m
LASTREV=$(cd $REPO && git svn log --oneline HEAD^^^^^..HEAD | head -1 | cut -f1 -d'|' | sed 's/r//')
echo 'loc = [ % all(/)  htdocs   dev   version  date' > $OUTF
for ver in `seq 1 $LASTREV`; do 
   echo "$REPO: VER=$ver";
   GITVER=`(cd $REPO && git svn find-rev r$ver)`;
   (cd $REPO && git checkout $GITVER);
   DATE=`(cd $REPO && git show -s --format=%ct $GITVER) | perl -pe'chomp'`;
   find $REPO -type f -name \*.m -print0 | xargs -0 perl \
       -e  'BEGIN{my $lines = 0};' \
       -ne 'next if /^\s*%/; next if /^\s*$/; $lines++;' \
       -e  'END{print "  $lines,";}' >> $OUTF ; # all (/)
   find $REPO/htdocs -type f -name \*.m -print0 | xargs -0 perl \
       -e  'BEGIN{my $lines = 0};' \
       -ne 'next if /^\s*%/; next if /^\s*$/; $lines++;' \
       -e  'END{print "  $lines,";}' >> $OUTF ; # htdocs
   find $REPO/dev -type f -name \*.m -print0 | xargs -0 perl \
       -e  'BEGIN{my $lines = 0};' \
       -ne 'next if /^\s*%/; next if /^\s*$/; $lines++;' \
       -e  'END{print "  $lines,";}' >> $OUTF ; # dev
   echo "  $ver,    $DATE; " >> $OUTF;
done
echo "];" >> $OUTF
