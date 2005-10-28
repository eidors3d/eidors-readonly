# Makefile for adler-lionheart-2005-EIDORS.tex
# Paper describing EIDORS v3.0
#
# $Id: Makefile,v 1.1 2005-10-28 02:31:30 aadler Exp $

NAME = adler-lionheart-2005-EIDORS

all: $(NAME).pdf

$(NAME).pdf: $(NAME).dvi
	dvipdf $<

$(NAME).dvi: $(NAME).tex $(FIGURES)
	latex $<
