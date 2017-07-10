all: EIT-review.pdf

pdflatex:=pdflatex -file-line-error -halt-on-error
latexmk:=latexmk -pdf -pdflatex="$(pdflatex)" -dependents

%.pdf: %.tex
	$(latexmk) $<

DEP= ../figs/fem.pdf\
    ../figs/fig07-offplane-sens-1plane_F.pdf\
    ../figs/fig07-offplane-sens-2pl_sqF.pdf\
    ../figs/fig08_conductivity_contrast.pdf\
    ../figs/fig10-hypoventilation-fraction.pdf\
    ../figs/hp-cb.png\
    ../figs/hp1.png\
    ../figs/hp21.png\
    ../figs/hp22.png\
    ../figs/hp23.png\
    ../figs/hp25.png\
    ../figs/hp31.png\
    ../figs/hp32.png\
    ../figs/hp33.png\
    ../figs/hp35.png\
    ../figs/hp41.png\
    ../figs/hp42.png\
    ../figs/hp43.png\
    ../figs/hp45.png\
    ../figs/sens-cb.png\
    ../figs/sens1.png\
    ../figs/sens2.png\
    ../figs/sens3.png\
    ../figs/skip11.png\
    ../figs/skip12.png\
    ../figs/skip13.png\
    ../figs/skip21.png\
    ../figs/skip22.png\
    ../figs/skip23.png\
    ../figs/stream1.pdf\
    ../figs/stream2.pdf\

EIT-review.pdf : $(DEP)

EIT-review.zip : EIT-review.tex $(DEP)
	-rm -rf $(basename $@) $@
	mkdir $(basename $@)
	cp $^ $(basename $@)/
	zip -r $@ $(basename $@)/
	rm -rf $(basename $@)/
	unzip $@
	cd $(basename $@)/ && $(latexmk) $(basename $@).tex
	rm -rf $(basename $@)/

clean:
	$(latexmk) -C EIT-review.tex
	-rm -f *-diff*.tex
	-rm -f ../figs/*
	-rm -f ../code/*.png
	-rm -f ../code/*.pdf
	-rm -f ../code/*.eps
	-rm -f ../code/ng.opt
	-rm -f EIT-review.zip

ifndef SKIPFIGS
# hpXX.png
../code/%1.png ../code/%21.png ../code/%22.png ../code/%23.png ../code/%24.png ../code/%25.png ../code/%31.png ../code/%32.png ../code/%33.png ../code/%34.png ../code/%35.png ../code/%41.png ../code/%42.png ../code/%43.png ../code/%44.png ../code/%45.png ../code/%-cb.png: ../code/%.m
	cd $(dir $<) && ./eidors -u $(notdir $(basename $<))

# skipXX.png
../code/%11.png ../code/%12.png ../code/%13.png ../code/%21.png ../code/%22.png ../code/%23.png: ../code/%.m
	cd $(dir $<) && ./eidors -u $(notdir $(basename $<))

# sensX.png, sens-cb.png
../code/%1.png ../code/%2.png ../code/%3.png ../code/%-cb.png: ../code/%.m
	cd $(dir $<) && ./eidors -u $(notdir $(basename $<))

# streamX.pdf
../code/%1.pdf ../code/%2.pdf: ../code/%.m
	cd $(dir $<) && ./eidors -u $(notdir $(basename $<))

# fem.pdf
../code/%.pdf: ../code/%.m
	cd $(dir $<) && ./eidors -u $(notdir $(basename $<))

# fig07-offplane-sens-Yplane_X.pdf
# ../figs/fig07-offplane-sens-2pl_sqF.pdf
../code/fig07-%-sens-1plane_F.pdf ../code/fig07-%-sens-2pl_sqF.pdf: ../code/%_sensitivity.m
	cd $(dir $<) && ./eidors -u $(notdir $(basename $<))

# png code/ -> figs/
../figs/%.png: ../code/%.png
	convert $< -trim $@

# pdf code/ -> figs/
../figs/%.pdf: ../code/%.pdf
	pdfcrop $< $@

%.pdf: %.odg
	soffice --headless --convert-to pdf $< --outdir $(dir $@)
endif
