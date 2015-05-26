N=GREIT-EIT2015

all: $(N).pdf

$(N).bbl: $(N).tex 
	pdflatex -halt-on-error $<
	bibtex $(N)

$(N).pdf: $(N).tex $(N).bbl
	pdflatex -halt-on-error $<


clean:
	rm $(N).pdf $(N).log $(N).bbl $(N).blg
