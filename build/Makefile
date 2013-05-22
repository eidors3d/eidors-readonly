# Build script for EIDORS Packages

# Need to have a local apache pointing to eidors

SZLIM= 20k
WKDIR = cd ..
REL= ../../../release/
INCLDIRS =  data_contrib \
            examples \
            paper-cheating \
            programming \
            tutorial 
WEBLINK= http://eidors3d.sourceforge.net/
#EIDORSDL=http://localhost:80/~adler/eidors/index.shtml
#CUTDL=1
 EIDORSDL=http://www.sce.carleton.ca/faculty/adler/eidors/index.shtml
 CUTDL=3

#all: $(REL) shtml other
all: $(REL)/htdocs

copy-to-mirror:
	@echo "Copy to Mirror ## Must run on adler-server"
	chgrp adler-res -R ../htdocs/
	chmod g+rw -R ../htdocs
	rsync --exclude=.svn -av ../htdocs/ /home/adler/htdocs/eidors/ 
	@echo "It's OK if some warnings"

$(REL):
	mkdir -p $(REL)
#	$(WKDIR) && mkdir -p $@

shtml:
	$(WKDIR) && find . -name \*.shtml | while read ff ; do \
	   dest="$(REL)/`dirname $$ff`/`basename $$ff .shtml`-s.html" ; \
	   echo "SHTML:$$ff ->$$dest"; \
	   if [  ! -f $$dest  ] ; then \
	      mkdir -p `dirname $$dest`; \
	      perl -pe's{\.shtml}{-s\.html}' $$ff > $$dest ; \
	   fi \
	done

other: 
	$(WKDIR) && find . \( -size -$(SZLIM) \) -o  \
	        \( -name \*.html -o -name \*.png \) \
           | while read ff ; do \
	   dest="$(REL)/$$ff"; \
	   echo "OTHER:$$ff ->$$dest"; \
	   if [  ! -f $$dest  ] ; then \
	      mkdir -p `dirname $$dest`; \
	      cp $$ff $$dest ; \
	   fi \
	done

#   if [  ! -f $$dest  -o  $$ff -nt $$dest  ] ; then \

jnk:

#         -e'print STDERR "$$2-$$_" if $$m{$$2};' \
#  -e'for ( keys %m) {print $$_,"=>",$$m{$$_} ."\n"}' 
#         -e'next unless m{href \s* = \s* . (^[$$q]*) \s* .}x;' \
          -e'print STDERR "$$r-$$_";}' \

#	mv $@/localhost*/eidors/* $@
#	rm -rf $@/localhost*
#         -e'next unless m{(href|src) \s* = \s* ($$q)(.*/)?([^\2]*?)(\2) }x;' \
#         -e'if (my $$r= $$m{$$4}) {s/$$2$$3$$4$$5/"$$r"/g};' \

$(REL)/htdocs: $(REL) Makefile
	( mkdir -p $@ && cd $@ && \
        wget -r -nH -np --cut-dirs=$(CUTDL) $(EIDORSDL) && \
	echo "################### 1 ####################" )  
	find $@ -name \*.shtml -print \
         | while read a ; do \
           mv $$a `dirname $$a`/`basename $$a .shtml`-s.html; \
         done
	echo "################### 2 ####################"
	find $@ -name \*.html -exec \
	 perl -i -pe's{\.shtml\s*"}{-s.html"}' \{} \;
	echo "################### 3 ####################"
	find $@ -name \*-s.html.bak | xargs rm -f
	echo "################### 3 ####################"
	cp $@/index-s.html $@/EIDORS-docs.html
	bigfiles=`find $@ -size +$(SZLIM) -and -not -iname \*html` ; \
	find $@ -name \*-s.html -exec \
	perl -i -pe'BEGIN{' -e'print STDERR "PROC: $$ARGV[0]\n";' \
             -e"@bigf = qw{$$bigfiles};" \
             -e'for (@bigf) {s{%5C}{/}g; $$s=$$t= $$_;' \
             -e  '$$t=~ s{^.*(data_contrib|tutorial|news_pics|examples)(/.*)}{$(WEBLINK)$$1$$2};'\
             -e  '$$s=~ s{^.*/(.*)}{$$1};'\
             -e  '$$m{$$s}= $$t;'\
          -e'} } my $$q= q{['\''"]}; $$n= q{^['\''"]};' \
          -e'next unless m{(href|src) \s* = \s* ($$q)(.*/)?([^\2]*?)(\2) }x;' \
          -e'next unless $$m{$$4};' \
          -e's{(href|src) \s* = \s* ($$q)(.*/)?([^\2]*?)(\2) }' \
          -e '{$$1 = "$$m{$$4}"}x;' \
        \{} \;
	find $@ -size +$(SZLIM) -and -not -iname \*html | xargs rm -f
	find $@ -name \*-s.html.bak | xargs rm -f

# Break out the steps
P = $(HOME)/eidors-release/htdocs
step0: 
	echo $P

step1:
	( mkdir -p $P && cd $P && \
        wget -r -nH -np --cut-dirs=$(CUTDL) $(EIDORSDL) && \
	echo "################### 1 ####################" )  

step2:
	find $P -name \*.shtml -print \
         | while read a ; do \
           mv $$a `dirname $$a`/`basename $$a .shtml`-s.html; \
         done
	echo "################### 2 ####################"

step3:
	find $P -name \*.html -exec \
	 perl -i -pe's{\.shtml\s*"}{-s.html"}' \{} \;
	echo "################### 3 ####################"

step4:
	find $P -name \*-s.html.bak | xargs rm -f
	echo "################### 3 ####################"

step5:
	cp $P/index-s.html $P/EIDORS-docs.html

step6:
	bigfiles=`find $P -size +$(SZLIM) -and -not -iname \*html` ; \
	find $P -name \*-s.html -exec \
	perl -i -pe'BEGIN{' -e'print STDERR "PROC: $$ARGV[0]\n";' \
             -e"@bigf = qw{$$bigfiles};" \
             -e'for (@bigf) {s{%5C}{/}g; $$s=$$t= $$_;' \
             -e  '$$t=~ s{^.*(data_contrib|tutorial|news_pics|examples)(/.*)}{$(WEBLINK)$$1$$2};'\
             -e  '$$s=~ s{^.*/(.*)}{$$1};'\
             -e  '$$m{$$s}= $$t;'\
          -e'} } my $$q= q{['\''"]}; $$n= q{^['\''"]};' \
          -e'next unless m{(href|src) \s* = \s* ($$q)(.*/)?([^\2]*?)(\2) }x;' \
          -e'next unless $$m{$$4};' \
          -e's{(href|src) \s* = \s* ($$q)(.*/)?([^\2]*?)(\2) }' \
          -e '{$$1 = "$$m{$$4}"}x;' \
        \{} \;

step7:
	find $P -size +$(SZLIM) -and -not -iname \*html | xargs rm -f
	find $P -name \*-s.html.bak | xargs rm -f

step8:
	svn export svn://svn.code.sf.net/p/eidors3d/code/trunk/eidors $P/../eidors

step9:
	echo "COPY LOCAL EIDORS BUILD DOC"

VERSION=v3.7RC
step10:
	mkdir $(REL)/eidors-$(VERSION)
	mv $P $(REL)/eidors $(REL)/eidors-$(VERSION)
	(cd $(REL) && zip -r eidors-$(VERSION).zip eidors-$(VERSION))
	

fix-cmyk-jpgs:
	find -iname \*.jpg | \
   while read ff ; do \
      if [ -n "$$(identify -verbose $$ff | grep -i 'colorspace: 4')" ]; then \
         echo $$ff ; \
         convert -quality 80 -colorspace RGB $$ff $$ff ; \
      fi ; done


# FIX adv hyperparams in build
# for fn in *html ; do perl -i -pe's{./adv_image_reconst/TV_compare_hyperparams.*html"}{http://eidors3d.sf.net/tutorial/adv_image_reconst/TV_compare_hyperparams.shtml"}' $fn ; done
# rm *html.bak
# rm -r adv_image_reconst/tv_hp_imgs

	

# TODO FOR RELEASE
# find -name \*.shtml | while read a  ; do mv $a `dirname $a`/`basename $a .shtml`.html ; done
# find -name \*.html -exec perl -i -pe's{\.shtml\s*"}{.html"}' \{} \;
# find -name \*.html.bak | xargs rm
# FIND LINKS THAT ARE TOO BIG (>0.5MB)
#  cg_deforming_tank_phantom
#  cg_normal_breathing
#  netgen_moving_ball
