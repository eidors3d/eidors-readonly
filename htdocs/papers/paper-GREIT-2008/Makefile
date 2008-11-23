# Makefile for GREIT paper
# $Id$
#
NAME= GREIT-algorithm

all: $(NAME).pdf

$(NAME).pdf: %.pdf: %.tex
	pdflatex $< && pdflatex $<

# EDIT the SVN version into the 
upload: $(NAME).pdf
	cd ../../../ && \
	rsync --exclude=.svn -av htdocs/ $$HOME/htdocs/eidors/
	FILEVER=`svn info $(NAME).tex | perl -n \
          -e'print "$$1 $$2:$$3" if ' \
          -e'/^Last .* (20\d{2}\-\d{2}\-\d{2}) +(\d{1,2}):(\d{2}):(\d{2})/;'` \
	&& echo $$FILEVER && \
	perl -i -pe"s/FILEVERSION/$$FILEVER/" \
          $$HOME/htdocs/eidors/GREIT/documents.shtml
	rsync  -avz $$HOME/htdocs/eidors/ aadler,eidors3d@web.sf.net:htdocs/

