all: eit2015.pdf

FIGURES = fig_loc.pdf mesh-eidors3p8.pdf
# provide a default repository location for EIDORS (git repo)
REPO = ../../..

eit2015.pdf: eit2015.tex $(FIGURES)
	pdflatex -halt-on-error $<

clean:
	-rm -f *.log *.bbl *.blg *.aux eit2015.pdf

loc.m:
	./git-loc ${REPO} > loc.m


EXCLUDE:=-not -name \*.png -not -name \*.gif -not -name \*.jpg -not -name \*.zip  -not -name \*.pdf -not -name \*.7z -not -name \*.jpeg -not -name \*.tomel -not -name \*.mat -not -name \*.pdf -not -name \*.mex\* -not -name \*.xsl -not -name \*.xsd -not -name \*.xml  -not -name \*.txt -not -name \*.frames -not -name \*.sframes -not -name \*.jar -not -name \*.msh -not -name \*.get -not -name \*.lyx
loc-num:
	@echo "eidors/"
	@find ${REPO}/eidors ${EXCLUDE} -type f | xargs wc -l | tail -1
	@echo "htodcs/"
	@find ${REPO}/htdocs ${EXCLUDE} -type f | xargs wc -l | tail -1
	@echo "dev/"
	@find ${REPO}/dev ${EXCLUDE} -type f | xargs wc -l | tail -1
