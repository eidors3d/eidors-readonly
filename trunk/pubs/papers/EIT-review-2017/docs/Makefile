PDFS=EIT-review.pdf 

 all: $(PDFS)

pdflatex:=pdflatex -file-line-error -halt-on-error -shell-escape
latexmk:=latexmk -pdf -pdflatex="$(pdflatex)" -bibtex

$(PDFS): %.pdf: %.tex ../figs/skip??.png ../figs/stream?.pdf ../figs/fem.pdf ../figs/sens?.png ../figs/sens-cb.png \
	     ../figs/hp1.png ../figs/hp??.png ../figs/hp-cb.png
	$(latexmk) $<


watch:
	$(latexmk) -pvc $N.tex

clean:
	$(latexmk) -C $N.tex
	-rm -f *-diff*.tex

ifdef FIGS
../code/%1.png ../code/%21.png ../code/%22.png ../code/%23.png ../code/%24.png ../code/%25.png ../code/%31.png ../code/%32.png ../code/%33.png ../code/%34.png ../code/%35.png ../code/%41.png ../code/%42.png ../code/%43.png ../code/%44.png ../code/%45.png ../code/%-cb.png: ../code/%.m
	cd $(dir $<) && ./eidors -u $(notdir $(basename $<))
../code/%11.png ../code/%12.png ../code/%13.png ../code/%21.png ../code/%22.png ../code/%23.png: ../code/%.m
	cd $(dir $<) && ./eidors -u $(notdir $(basename $<))
../code/%1.png ../code/%2.png ../code/%3.png ../code/%-cb.png: ../code/%.m
	cd $(dir $<) && ./eidors -u $(notdir $(basename $<))
../code/%1.pdf ../code/%2.pdf: ../code/%.m
	cd $(dir $<) && ./eidors -u $(notdir $(basename $<))
../code/%.pdf: ../code/%.m
	cd $(dir $<) && ./eidors -u $(notdir $(basename $<))
../figs/%.png: ../code/%.png
	convert $< -trim $@
../figs/%.pdf: ../code/%.pdf
	pdfcrop $<
	mv $(basename $<)-crop.pdf $@
endif
