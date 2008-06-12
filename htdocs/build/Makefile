# Build script for EIDORS Packages

# Need to have a local apache pointing to eidors

SZLIM= 100k
WKDIR = cd ..
REL= ../../release/
INCLDIRS =  data_contrib \
            examples \
            paper-cheating \
            programming \
            tutorial 
WEBLINK= http://eidors3d.sourceforge.net/
#EIDORSDL=http://localhost:80/eidors/index.shtml
EIDORSDL=http://www.sce.carleton.ca/faculty/adler/eidors/index.shtml
CUTDL=3

#all: $(REL) shtml other
all: $(REL)/htdocs

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

$(REL)/htdocs: $(REL) Makefile
	mkdir -p $@ && cd $@ && \
        wget -r -nH --cut-dirs=$(CUTDL) $(EIDORSDL)
	find $@ -name \*.shtml -print \
         | while read a ; do \
           mv $$a `dirname $$a`/`basename $$a .shtml`-s.html; \
         done
	find $@ -name \*.html -exec \
	 perl -i -pe's{\.shtml\s*"}{-s.html"}' \{} \;
	find $@ -name \*-s.html.bak | xargs rm
	cp $@/index-s.html $@/EIDORS-docs.html
	bigfiles=`find $@ -size +$(SZLIM)` ; \
	find $@ -name \*-s.html -exec \
	perl -i -pe'BEGIN{' -e'print STDERR "PROC: $$ARGV[0]\n";' \
             -e"@bigf = qw{$$bigfiles};" \
             -e'for (@bigf) {s{%5C}{/}g; $$s=$$t= $$_;' \
             -e  '$$t=~ s{^.*(data_contrib|tutorial|examples)(/.*)}{$(WEBLINK)$$1$$2};'\
             -e  '$$s=~ s{^.*/(.*)}{$$1};'\
             -e  '$$m{$$s}= $$t;'\
          -e'} } my $$q= q{['\''"]}; $$n= q{^['\''"]};' \
          -e'next unless m{(href|src) \s* = \s* ($$q)(.*/)?([^\2]*)(\2) }x;' \
          -e'if (my $$r= $$m{$$4}) {s/$$2$$3$$4$$5/"$$r"/g};' \
        \{} \;
	find $@ -size +$(SZLIM) | xargs rm 
	find $@ -name \*-s.html.bak | xargs rm
	

# TODO FOR RELEASE
# find -name \*.shtml | while read a  ; do mv $a `dirname $a`/`basename $a .shtml`.html ; done
# find -name \*.html -exec perl -i -pe's{\.shtml\s*"}{.html"}' \{} \;
# find -name \*.html.bak | xargs rm
# FIND LINKS THAT ARE TOO BIG (>0.5MB)
#  cg_deforming_tank_phantom
#  cg_normal_breathing
#  netgen_moving_ball
