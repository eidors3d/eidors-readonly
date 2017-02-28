# Build EIDORS web site and move to sourceforge

BUILD_DOCS = ../build-htdocs

#SHTML = $(shell find -type f -name '*.shtml')
#THTML:= $(SHTML,%.shtml=$(BUILD_DOCS)/%.html)
# Seems need to reorder variable assignment
 THTML = $(patsubst %.shtml, $(BUILD_DOCS)/%.html, \
             $(shell find -type f -name '*.shtml'))

all: $(BUILD_DOCS) $(THTML)


#tpage $< > $@ #perl -MTemplate -e'$$t->process("$<") or die Template->error();' > $@
$(THTML): $(BUILD_DOCS)/%.html: %.shtml
	@echo $< " => " $@
	mkdir -p `dirname $@`
	tpage $< > $@ 



#OLD CONTENT:
#EIDORSUSER= aadler
#
#website: $(PDFS)
#	rsync --exclude=.svn  -avz . $(EIDORSUSER)@web.sourceforge.net:/home/project-web/eidors3d/htdocs
#
#PDFS= papers/paper-GREIT-2008/GREIT-algorithm.pdf
#$(PDFS): %.pdf: %.tex
#	(cd `dirname $<` && \
#         pdflatex `basename $<` && \
#         pdflatex `basename $<` )

