# Build EIDORS web site and move to sourceforge
#

EIDORSUSER= aadler

website: pdfs
	rsync --exclude=.svn  -avz . $(EIDORSUSER),eidors3d@web.sf.net:htdocs/

pdfs: $(PDFS)

PDFS= papers/paper-GREIT-2008/GREIT-algorithm.pdf
$(PDFS): %.pdf: %.tex
	pdflatex $<
	pdflatex $<

