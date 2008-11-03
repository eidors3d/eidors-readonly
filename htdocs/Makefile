# Build EIDORS web site and move to sourceforge
#

EIDORSUSER= aadler

website: $(PDFS)
	rsync --exclude=.svn  -avz . $(EIDORSUSER),eidors3d@web.sf.net:htdocs/

PDFS= papers/paper-GREIT-2008/GREIT-algorithm.pdf
$(PDFS): %.pdf: %.tex
	(cd `dirname $<` && \
         pdflatex `basename $<` && \
         pdflatex `basename $<` )

