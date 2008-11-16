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
         $(FIGDIR)/example_images

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

if_data_2003.zip:
	wget http://eidors3d.sf.net/data_contrib/if-peep-acute-lung-injury/if_data_2003.zip

testdata.mat: make_testdata.m if_data_2003.zip
	unzip -o if_data_2003.zip p1130107.get
	$(CALLMATLABPARAMS)

ALGS = Recon_GREIT_ndiff.mat \
       Recon_NOSER_diff.mat \
       Recon_NOSER_ndiff.mat \
       Recon_Sheffield_backproj.mat

Recon_NOSER_diff.mat Recon_NOSER_ndiff.mat: jacobian_cyl.mat
Recon_GREIT_ndiff.mat: sim_targets.mat

$(ALGS): %.mat: %.m
	$(CALLMATLABPARAMS)

# EXAMPLE IMAGES:
$(FIGDIR)/example_images: testdata.mat

$(FIGDIR)/example_images: example_images.m $(ALGS)
	touch $@
	$(CALLMATLABPARAMS2)
