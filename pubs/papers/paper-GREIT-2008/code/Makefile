# Create Simulated data for GREIT evaluations
# Usage set MATLAB 
#
# (C) Andy Adler 2008 $Id$
# Licensed under GNU GPL v2 or v3


MATLAB := $(shell which matlab 2>/dev/null)
ifdef MATLAB
   MATLAB += -nodesktop
else
   MATLAB := C:/progra~1/matlab6p5/bin/win32/matlab.exe
endif

   EIDORSTARTUP = run('../../../../eidors/startup.m')

   CM1= -r "$(EIDORSTARTUP);
   CMATLAB= fn=`basename $< .m` && $(MATLAB) -nosplash $(CM1)
   CME= ;exit"
   CALLMATLAB=        $(CMATLAB) $$fn $(CME)
   CALLMATLABPARAMS = $(CMATLAB) $$fn('$@') $(CME)
   CALLMATLABPARAMS2= $(CMATLAB) $$fn('$@','$^') $(CME)
#TRIMFIG = convert -trim $@ PNG8:$@


FIGDIR = ../figures

PRECALC= ng_cyl_mdl.mat sim_targets.mat jacobian_cyl.mat
FIGURES= $(FIGDIR)/fig_cyl_mdl_targets.png \
         $(FIGDIR)/example_images \
         $(FIGDIR)/female_t_mdl.png \
         $(FIGDIR)/male_t_mdl.png

all: $(FIGURES)


# CREATE NETGEN TANK MODEL
ng_cyl_mdl.mat: make_ng_cyl_mdl.m
	$(CALLMATLAB)

# SIMULATE HOMOG MOVEMENTS
sim_targets.mat: sim_targets.m ng_cyl_mdl.mat
	$(CALLMATLABPARAMS)


# CALCULATE JACOBIAN OF TANK
jacobian_cyl.mat: calc_jacobian_mdl.m ng_cyl_mdl.mat
	$(CALLMATLAB)

# IMAGE OF SIMULATED TARGETS
$(FIGDIR)/fig_cyl_mdl_targets.png: sim_targets.mat

$(FIGDIR)/fig_cyl_mdl_targets.png: fig_cyl_mdl_targets.m
	$(CALLMATLABPARAMS)
	$(TRIMFIG)

EXPDATA = p1130107.get 1-control.RAW 1-injury.RAW goev354005.get goev354003_512-1560.get goev354008_947-1559.get

p1130107.get:
	wget http://eidors3d.sf.net/data_contrib/if-peep-acute-lung-injury/if_data_2003.zip
	unzip if_data_2003.zip $@

1-control.RAW 1-injury.RAW:
	wget http://eidors3d.sf.net/data_contrib/if-experimental-lung-injury/if-experimental-lung-injury.zip
	unzip if-experimental-lung-injury.zip $@

goev354005.get goev354003_512-1560.get goev354008_947-1559.get:
	wget http://eidors3d.sf.net/data_contrib/gh_pleural_cavity/gh_data_2006.zip
	unzip -j gh_data_2006.zip $(patsubst %.get,*%.get,$@)

sim_testdata.mat: make_sim_testdata.m
	$(CALLMATLABPARAMS)

testdata.mat: make_testdata.m $(EXPDATA)
	$(CALLMATLABPARAMS)

ALGS = Recon_GREIT_ndiff.mat \
       Recon_NOSER_diff.mat \
       Recon_NOSER_ndiff.mat \
       Recon_Sheffield_backproj.mat

Recon_NOSER_diff.mat:  jacobian_cyl.mat
Recon_NOSER_ndiff.mat: jacobian_cyl.mat
Recon_GREIT_ndiff.mat: sim_targets.mat

$(ALGS): %.mat: %.m
	$(CALLMATLABPARAMS)

# EXAMPLE IMAGES:
$(FIGDIR)/example_images: testdata.mat sim_testdata.mat

$(FIGDIR)/example_images: example_images.m $(ALGS)
	touch $@
	$(CALLMATLABPARAMS2)

# male_t_mdl
HMT027.mat: HMT027.zip
	unzip $< 

male_t_mdl.mat: make_thorax_mdl.m HMT027.mat
	$(CALLMATLABPARAMS)

$(FIGDIR)/male_t_mdl.png: print_thorax_model.m male_t_mdl.mat
	$(CALLMATLABPARAMS2)
	convert -trim $@ png8:$@

# female_t_mdl
HFT027.mat: HFT027.zip
	unzip $< 

female_t_mdl.mat: make_thorax_mdl.m HFT027.mat
	$(CALLMATLABPARAMS)

$(FIGDIR)/female_t_mdl.png: print_thorax_model.m female_t_mdl.mat
	$(CALLMATLABPARAMS2)
	convert -trim $@ png8:$@
