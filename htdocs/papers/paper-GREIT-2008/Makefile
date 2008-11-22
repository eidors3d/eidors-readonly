# Makefile for GREIT paper
# $Id$

all: GREIT-algorithm.pdf

GREIT-algorithm.pdf: %.pdf: %.tex
	pdflatex $< && pdflatex $<

# EDIT the SVN version into the 
upload:
	cd ../../../ && \
	rsync --exclude=.svn -av htdocs/ $$HOME/htdocs/eidors/
	FILEVER=`svn info $< | perl -n \
          -e'print "$$1 $$2:$$3" if ' \
          -e'/^Last .* (20\d{2}\-\d{2}\-\d{2}) +(\d{1,2}):(\d{2}):(\d{2})/;'` \
	&& echo $$FILEVER && \
	perl -i -pe"s/FILEVERSION/$$FILEVER/" \
          $$HOME/htdocs/eidors/GREIT/documents.shtml
	rsync $$HOME/htdocs/eidors/ aadler,eidors3d@web.sf.net:htdocs/

