# Makefile for adler-lionheart-2005-EIDORS.tex
# Paper describing EIDORS v3.0
#
# $Id: Makefile,v 1.2 2005-10-28 03:40:53 aadler Exp $

NAME = adler-lionheart-2005-EIDORS

FIGURES = eidors-logo.eps

all: $(NAME).pdf

$(NAME).pdf: $(NAME).dvi
	dvipdf $<

$(NAME).dvi: $(NAME).tex $(FIGURES)
	latex $<

eidors-logo.eps: ../eidors-logo.jpg
	convert $< eps2:$@
