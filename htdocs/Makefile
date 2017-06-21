# Build EIDORS web site and move to sourceforge
# Build options
#   all: build 

BUILD_DOCS = ../build-htdocs/

THTML = $(filter-out nav-sidebar.shtml, \
            $(patsubst %.shtml, $(BUILD_DOCS)/%.shtml, \
            $(shell find -type f -name '*.shtml')))

RAW   = $(patsubst %, $(BUILD_DOCS)/%, \
            $(shell find -type f -iname '*.png' \
                             -o  -iname '*.jpg' \
                             -o  -iname '*.jpeg' \
                             -o  -iname '*.pdf' \
                             -o  -iname '*.zip' \
                             -o  -iname '*.mat' \
          ) \
          tutorial/data_structures/tutorial020_f_solve.m \
          tutorial/workshop/code_2june2015.m \
          tutorial/workshop/code_21june2017.m \
        )

INDEX = $(BUILD_DOCS)/index.html



all: $(THTML) $(RAW) $(INDEX)

$(INDEX): $(BUILD_DOCS)/index.shtml
	cp $< $@

$(THTML): nav-sidebar.shtml

#tpage $< > $@ #perl -MTemplate -e'$$t->process("$<") or die Template->error();' > $@
$(THTML): $(BUILD_DOCS)/%.shtml: %.shtml
	@echo $< " => " $@
	mkdir -p `dirname $@`
	(cd `dirname $<` && tpage `basename $<`) > $@ 

$(RAW): $(BUILD_DOCS)/%: %
	@echo $< " => " $@
	mkdir -p `dirname $@`
	cp $< $@ 




EIDORSUSER= aadler
upload: all
	rsync --exclude=.svn -rltzuv -O --no-perms $(BUILD_DOCS) $(EIDORSUSER)@web.sourceforge.net:/home/project-web/eidors3d/htdocs/

# ONLY AFTER RUNNING DOC GENERATE SCRIPT build_doc
# no trailing slash
upload-doc: all
	rsync --exclude=.svn -rltzuv -O --no-perms ./doc $(EIDORSUSER)@web.sourceforge.net:/home/project-web/eidors3d/htdocs/

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

