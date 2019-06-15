pdflatex:=pdflatex -file-line-error -halt-on-error -shell-escape
latexmk:=latexmk -pdf -f -pdflatex="$(pdflatex)" -bibtex

all: eidors3p10-eit2019.pdf
%.pdf: %.tex
	$(latexmk) $<

# Graph is generated as follows
#   Get a git svn clone of eidors
#   run git-loc.sh (. git-loc.sh)
#   in matlab, run graph_loc.m
