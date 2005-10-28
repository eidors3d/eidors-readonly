# Makefile for adler-lionheart-2005-EIDORS.tex
# Paper describing EIDORS v3.0
#
# $Id: Makefile,v 1.3 2005-10-28 04:29:22 aadler Exp $

NAME = adler-lionheart-2005-EIDORS


CHEAT_FIGS= deform-grille.eps deform.eps \
    happy.eps happynoise.eps sad.eps \
    laplace0-icrime.eps laplace3-icrime.eps \
    laplace3.eps tikhonov-icrime.eps \
    tikhonov.eps

FIGURES = eidors-logo.eps $(CHEAT_FIGS)

all: $(NAME).pdf

$(NAME).pdf: $(NAME).dvi
	dvipdf $<

$(NAME).dvi: $(NAME).tex $(FIGURES)
	latex $<

eidors-logo.eps: ../eidors-logo.jpg
	convert $< eps2:$@

$(CHEAT_FIGS): %.eps: ../paper-cheating/%.png
	convert $< eps2:$@

clean: 
	-rm $(NAME).dvi $(NAME).log $(NAME).aux $(NAME).pdf
	-rm $(FIGURES)
