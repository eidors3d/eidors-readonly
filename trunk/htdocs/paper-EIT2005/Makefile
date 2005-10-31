# Makefile for adler-lionheart-2005-EIDORS.tex
# Paper describing EIDORS v3.0
#
# $Id: Makefile,v 1.4 2005-10-31 03:36:59 aadler Exp $

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

fig_3d:
	echo \
n3r2= mk_common_model('n3r2'); 
n3r2.fwd_model.nodes = n3r2.fwd_model.nodes/3;
fmdl_n3r2 = n3r2.fwd_model;
sigma= ones( size(fmdl_n3r2.elems,1) ,1);
img_n3r2= eidors_obj('image', 'homog', 'elem_data', sigma, ...
                'fwd_model', fmdl_n3r2 );
vh_n3r2= fwd_solve( img_n3r2);

n3z= mk_common_model('n3z'); 
n3z.fwd_model.nodes = n3z.fwd_model.nodes/3;
fmdl_n3z = n3z.fwd_model;
img_n3z= eidors_obj('image', 'homog', 'elem_data', sigma, ...
                'fwd_model', fmdl_n3z );
vh_n3z= fwd_solve( img_n3z );

load('datacom.mat','A','B');
sigma(A)= 1.2; sigma(B)= 0.8;
img_n3r2.elem_data= sigma; 
vi_n3r2= fwd_solve( img_n3r2 );
noise= 0.1*std( vh_n3r2.meas - vi_n3r2.meas )*randn( size(vi_n3r2.meas) );
vi_n3r2.meas= vi_n3r2.meas + noise;

img_n3z.elem_data= sigma; 
vi_n3z= fwd_solve( img_n3z );
noise= 0.1*std( vh_n3z.meas - vi_n3z.meas )*randn( size(vi_n3z.meas) );
vi_n3z.meas= vi_n3z.meas + 0*noise;

imn3r2= inv_solve( n3r2, vh_n3r2, vi_n3r2);
imn3z= inv_solve( n3z, vh_n3z, vi_n3z);


