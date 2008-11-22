# Makefile for GREIT paper
# $Id$

all: GREIT-algorithm.pdf

GREIT-algorithm.pdf: %.pdf: %.tex
	pdflatex $< && pdflatex $<

upload:
	cd ../../../ && \
	rsync --exclude=.svn -av htdocs/ $$HOME/htdocs/eidors/
	rsync $$HOME/htdocs/eidors/ aadler,eidors3d@web.sf.net:htdocs/

