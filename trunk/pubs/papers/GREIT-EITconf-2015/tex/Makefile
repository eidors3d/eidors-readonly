N=GREIT-EIT2015
P=GREIT-EIT2015-talk

all: $(N).pdf $(P).pdf

$(N).bbl: $(N).tex 
	pdflatex -halt-on-error $<
	bibtex $(N)

$(N).pdf: $(N).tex $(N).bbl
	pdflatex -halt-on-error $<


$(P).pdf: $(P).tex
	pdflatex -halt-on-error $<

clean:
	rm $(N).pdf $(N).log $(N).bbl $(N).blg
	rm $(P).pdf $(P).log $(P).bbl $(P).blg
