# Build script for EIDORS Packages

# Need to have a local apache pointing to eidors

SZLIM= 20k
WKDIR = cd ..
REL= $(HOME)/release/
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
	rsync --exclude=.svn -av -O ../htdocs/ /home/adler/htdocs/eidors/ 
	@echo "It's OK if some warnings"

copy-mirror-to-sf:
	rsync --exclude=.svn -avz -O /home/adler/htdocs/eidors/ aadler,eidors3d@web.sf.net:htdocs/

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
	@echo "##### ENSURE YOU HAVE COPIED TO THE MIRROR ON sce.carleton.ca ####"

# download eidors though a shtml web site to get the dyn html 
step1:
	( mkdir -p $P && cd $P && \
        wget -r -nH -np --cut-dirs=$(CUTDL) $(EIDORSDL) && \
	echo "################### 1 ####################" )  

# rename shtml to -s.html files
step2:
	find $P -name \*.shtml -print \
         | while read a ; do \
           mv $$a `dirname $$a`/`basename $$a .shtml`-s.html; \
         done
	@echo "################### 2 ####################"

# rewrite internal links so  shtml to -s.html files
step3:
	find $P -name \*.html -exec \
	 perl -i -pe's{\.shtml\s*"}{-s.html"}' \{} \;
	@echo "################### 3 ####################"

# fix rewrite for the mirror link
step3a:
	find $P -name \*.html -exec \
	 perl -i -pe's{(http://\S*)-s\.html"}{\1.shtml"}' \{} \;
	@echo "################### 3 ####################"

# get rid of the bak files that are created in step3
step4:
	find $P -name \*-s.html.bak | xargs rm -f
	@echo "################### 3 ####################"

# create and EIDORS-docs.html files
step5:
	cp $P/index-s.html $P/EIDORS-docs.html

# create web links to big files
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

# erase too big files
step7:
	find $P -type f -size +$(SZLIM) -and -not -iname \*html | xargs rm -f
	find $P -name \*-s.html.bak | xargs rm -f

# export the latest eidors svn
# export is like checkout except no .svn directories
step8:
	rm -rf $P/../eidors
	svn export svn://svn.code.sf.net/p/eidors3d/code/trunk/eidors $P/../eidors

# do not uncomment cd ~/eidors-release
step9:
	@echo "BUILD LOCAL DOCS: run build_doc.m making sure VERSION=1"

step9a:
	cp README $P/..
	cp -r ../htdocs/doc $P

VERSION=v3.10
step10:
	mkdir -p $(REL)
	rm -f $P/robots.txt
	rm -fr $(REL)/eidors-$(VERSION)
	rm -f $(REL)/eidors-$(VERSION).zip
	cp -r $P/.. $(REL)/eidors-$(VERSION)
	(cd $(REL) && zip -r eidors-$(VERSION).zip eidors-$(VERSION))

step11:
	cp $(REL)/eidors-$(VERSION).zip $(HOME)/public_html

ngversion: # Just a sketch. Needs a fix
	mv ~/release/eidors-v3.9 ~/release/eidors-v3.9-ng
	cp README+NG ~/release/eidors-v3.9-ng/README
	cp -r ~/Downloads/eidors-v3.8-ng/Netgen-5.3_x64/ ~/release/eidors-v3.9-ng
	(cd ~/release/ && zip -r eidors-v3.9-ng.zip eidors-v3.9-ng)


# update as required
tag-release:
	svn copy svn+ssh://aadler@svn.code.sf.net/p/eidors3d/code/trunk svn+ssh://aadler@svn.code.sf.net/p/eidors3d/code/tags/$(VERSION) -m"Tagging release $(VERSION)"

	

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
