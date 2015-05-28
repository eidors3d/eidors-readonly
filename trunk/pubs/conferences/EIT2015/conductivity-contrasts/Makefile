N= eit2015-conductivity-contrasts

P= eit2015-conductivity-contrasts-pres


IMGS= contrasts_04b.pdf \
      contrasts_04p.pdf \
      contrasts_04h.pdf \
      contrasts_04v.pdf \
      signal.pdf 

all: $(N).pdf $(P).pdf

# Create images by running figure1.m
images: $(IMGS)
	for i in $(IMGS) ; do pdfcrop $$i $$i ; done

 # $(IMGS)
$(N).pdf: $(N).tex 
	pdflatex -halt-on-error $(N)


$(P).pdf: $(P).tex 
	pdflatex -halt-on-error $(P)
