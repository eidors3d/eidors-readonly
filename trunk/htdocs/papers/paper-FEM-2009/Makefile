# Makefile for GREIT paper
# $Id$
#
NAME= FEM-Errors
LATEXPART= iopart.cls iopart10.clo iopart12.clo

#all: $(NAME).pdf
all: $(NAME).dvi

$(LATEXPART):
	wget ftp://ftp.iop.org/pub/journals/ioplatexguidelines.zip
	unzip ioplatexguidelines.zip $^        

$(NAME).dvi: $(LATEXPART)

$(NAME).dvi: %.dvi: %.tex
	latex $< && latex $<

$(NAME).pdf: %.pdf: %.tex
	pdflatex $< && pdflatex $<

