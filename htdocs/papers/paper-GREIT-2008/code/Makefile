# Create Simulated data for GREIT evaluations
# Usage set MATLAB 
#
# (C) Andy Adler 2008 $Id$
# Licensed under GNU GPL v2 or v3

ifndef MATLAB
   MATLAB = C:/progra~1/matlab6p5/bin/win32/matlab.exe
endif
   EIDORSTARTUP = run('../../../../eidors/startup.m')
   ALGSLIST = algs = { \
          'GREIT_Sheffield_backproj', \
          'GREIT_NOSER_ndiff', \
          'GREIT_NOSER_diff', \
          'GREIT_test_ndiff', \
         };

   CM1= -r "$(EIDORSTARTUP);
   CMATLAB= fn=`basename $< .m` && $(MATLAB) -nosplash $(CM1)
   CME= ;exit"
   CALLMATLAB=       $(CMATLAB) $$fn $(CME)
   CALLMATLABPARAMS= $(CMATLAB) $$fn('$@') $(CME)
   CALLMATLABALGS=   $(CMATLAB) $(ALGSLIST) $$fn $(CME)
#TRIMFIG = convert -trim $@ PNG8:$@


FIGDIR = ../figures

PRECALC= ng_cyl_mdl.mat sim_targets.mat jacobian_cyl.mat
FIGURES= $(FIGDIR)/fig_cyl_mdl_targets.png

all: $(FIGURES)

# CREATE NETGEN TANK MODEL
ng_cyl_mdl.mat: make_ng_cyl_mdl.m
	$(CALLMATLAB)

# SIMULATE HOMOG MOVEMENTS
sim_targets.mat: sim_targets.m ng_cyl_mdl.mat
	$(CALLMATLAB)


# CALCULATE JACOBIAN OF TANK
jacobian_cyl.mat: calc_jacobian_mdl.m ng_cyl_mdl.mat
	$(CALLMATLAB)

# IMAGE OF SIMULATED TARGETS
$(FIGDIR)/fig_cyl_mdl_targets.png: fig_cyl_mdl_targets.m \
               $(PRECALC)
	$(CALLMATLABPARAMS)
	$(TRIMFIG)

# TEST DATA
testdata.mat: make_testdata.m sim_targets.mat
	$(CALLMATLABPARAMS)
