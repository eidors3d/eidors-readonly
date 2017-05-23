pdflatex:=pdflatex -file-line-error -halt-on-error -shell-escape
latexmk:=latexmk -pdf -pdflatex="$(pdflatex)" -bibtex

all: eidors3p9-eit2017.pdf
eidors3p9-eit2017.pdf: *.png
%.pdf: %.tex
	$(latexmk) $<
