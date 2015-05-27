#!/bin/sh
# Sequentially checkout versions 
REPO="$HOME/docs/eidors-git";
OUTF=datafile.txt
#echo > $OUTF
for ver in `seq 4951 4970`; do 
   echo "$REPO: VER=$ver";
   GITVER=`(cd $REPO && git svn find-rev r$ver)`;
   (cd $REPO && git checkout $GITVER);
   DATE=`(cd $REPO && git show -s --format=%ct $GITVER) | perl -pe'chomp'`;
   echo -n "$REPO ($ver,$DATE) [" >> $OUTF;
   find $REPO -type f -name \*.m -print0 | xargs -0 perl \
       -e  'BEGIN{my $lines = 0};' \
       -ne 'next if /^\s*%/; next if /^\s*$/; $lines++;' \
       -e  'END{print "  $lines(/)";}' >> $OUTF
   find $REPO/htdocs -type f -name \*.m -print0 | xargs -0 perl \
       -e  'BEGIN{my $lines = 0};' \
       -ne 'next if /^\s*%/; next if /^\s*$/; $lines++;' \
       -e  'END{print "  $lines(htdocs)";}' >> $OUTF
   find $REPO/dev -type f -name \*.m -print0 | xargs -0 perl \
       -e  'BEGIN{my $lines = 0};' \
       -ne 'next if /^\s*%/; next if /^\s*$/; $lines++;' \
       -e  'END{print "  $lines(dev)";}' >> $OUTF
   echo "]" >> $OUTF
done
