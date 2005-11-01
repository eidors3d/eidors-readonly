# Makefile for adler-lionheart-2005-EIDORS.tex
# Paper describing EIDORS v3.0
#
# $Id: Makefile,v 1.9 2005-11-01 02:43:04 aadler Exp $

NAME = adler-lionheart-2005-EIDORS

SET_COLOUR= global eidors_colours; \
            eidors_colours.greylev= 0.1; \
            eidors_colours.sat_adj= 0.98; \
            eidors_colours.backgnd= [1,1,1]; \
            eidors_colours.mapped_colour= 256;

CHEAT_FIGS= deform-grille.eps deform.eps \
    happy.eps happynoise.eps sad.eps \
    laplace0-icrime.eps laplace3-icrime.eps \
    laplace3.eps tikhonov-icrime.eps \
    tikhonov.eps

FIGURES = fig_tik_w-icrime.eps \
          fig_tik_wo-icrime.eps \
          fig_lapl_w-icrime.eps \
          fig_lapl_wo-icrime.eps \
	  fig_3d.eps

all: $(NAME).pdf

$(NAME).pdf: $(NAME).dvi
	dvipdf $<

$(NAME).dvi: $(NAME).tex $(FIGURES)
	latex $<

eidors-logo.eps: ../eidors-logo.jpg
	convert $< eps2:$@

#$(CHEAT_FIGS): %.eps: ../paper-cheating/%.png
#	convert $< eps2:$@

clean: 
	-rm $(NAME).dvi $(NAME).log $(NAME).aux $(NAME).pdf
	-rm $(FIGURES)

fig_noise.m: Makefile
	fname=`basename $@ .m`; \
	echo \
  "$(SET_COLOUR)" \
  "cheating_2d('1',[1702 7901 10801 11204]);" \
  "set(gcf, 'PaperPosition', [0.25,2.5, 4,2]);" \
  "print( gcf, '-depsc2','-painters', '$$fname.eps');" \
      > $@

fig_noise.eps: fig_noise.m
	matlab -r "run('../../eidors3d/startup'); run('./$<'); exit"

STARTUP = run('../../eidors3d/startup'); $(SET_COLOUR)
PRINTEPS= set(gcf, 'PaperPosition', [0.25,2.5, 4,2]); print( gcf, '-depsc2','-painters', '

fig_tik_w-icrime.eps:
	matlab -r "$(STARTUP); cheating_2d('2a'); $(PRINTEPS)$@'); exit"

fig_tik_wo-icrime.eps:
	matlab -r "$(STARTUP); cheating_2d('2b'); $(PRINTEPS)$@'); exit"

fig_lapl_w-icrime.eps:
	matlab -r "$(STARTUP); cheating_2d('3a'); $(PRINTEPS)$@'); exit"

fig_lapl_wo-icrime.eps:
	matlab -r "$(STARTUP); cheating_2d('3b'); $(PRINTEPS)$@'); exit"

fig_3d.m: Makefile
	fname=`basename $@ .m`; \
	echo \
  "$(SET_COLOUR)" \
  \
  "n3r2= mk_common_model('n3r2'); " \
  "n3r2.fwd_model.nodes = n3r2.fwd_model.nodes/3;" \
  "fmdl_n3r2 = n3r2.fwd_model;" \
  "sigma= ones( size(fmdl_n3r2.elems,1) ,1);" \
  "img_n3r2= eidors_obj('image', 'homog', 'elem_data', sigma, " \
  "                'fwd_model', fmdl_n3r2 );" \
  "vh_n3r2= fwd_solve( img_n3r2);" \
  \
  "n3z= mk_common_model('n3z'); " \
  "n3z.fwd_model.nodes = n3z.fwd_model.nodes/3;" \
  "n3z.hyperparameter.value= 1e-4;" \
  "fmdl_n3z = n3z.fwd_model;" \
  "img_n3z= eidors_obj('image', 'homog', 'elem_data', sigma, " \
  "                'fwd_model', fmdl_n3z );" \
  "vh_n3z= fwd_solve( img_n3z );" \
  \
  "load('datacom.mat','A','B');" \
  "sigma(A)= 1.2; sigma(B)= 0.8;" \
  "img_n3r2.elem_data= sigma; " \
  "vi_n3r2= fwd_solve( img_n3r2 );" \
  "noise= 0.1*std( vh_n3r2.meas - vi_n3r2.meas )*randn( size(vi_n3r2.meas) );" \
  "vi_n3r2.meas= vi_n3r2.meas + noise;" \
  \
  "img_n3z.elem_data= sigma; " \
  "vi_n3z= fwd_solve( img_n3z );" \
  "noise= 0.1*std( vh_n3z.meas - vi_n3z.meas )*randn( size(vi_n3z.meas) );" \
  "vi_n3z.meas= vi_n3z.meas + noise;" \
  \
  "imn3r2= inv_solve( n3r2, vh_n3r2, vi_n3r2);" \
  "imn3z= inv_solve( n3z, vh_n3z, vi_n3z);" \
  \
  "clf;" \
  "axes('position',[.1,.05,.5,.95]);" \
  "show_fem(imn3z); view([-28 12]);" \
  "axes('position',[.6,0,.4,1]);" \
  "show_slices(imn3z,[(1:2:7)'*[inf,inf,1/8], ones(4,1), (4:-1:1)']);" \
  "set(gcf, 'PaperPosition', [0.25,2.5, 4,3]);" \
  "print( gcf, '-depsc2','-painters', '$$fname.eps');" \
      > $@

fig_3d.eps: fig_3d.m
	matlab -r "run('../../eidors3d/startup'); run('./$<'); exit"

