# Create Simulated data for GREIT evaluations
# Usage set MATLAB 
#
# (C) Andy Adler 2008 $Id$
# Licensed under GNU GPL v2 or v3

ifndef MATLAB
   MATLAB = C:/progra~1/matlab6p5/bin/win32/matlab.exe
   EIDORSTARTUP = run('../../../../eidors/startup.m')

   CM1= -r "$(EIDORSTARTUP);
   CMATLAB= fn=`basename $< .m` && $(MATLAB) -nosplash $(CM1)
   CME= ;exit"
   CALLMATLAB= $(CMATLAB) $$fn $(CME)
endif

PRECALC= ng_cyl_mdl.mat sim_targets.mat jacobian_cyl.mat
FIGURES= fig_cyl_mdl_targets.png

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
fig_cyl_mdl_targets.png: fig_cyl_mdl_targets.m \
               $(PRECALC)
	$(CALLMATLAB)
