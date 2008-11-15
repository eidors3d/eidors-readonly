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

PRECALC= ng_cyl_mdl.mat sim_targets.mat GREIT_Jacobian_cyl.mat
FIGURES= simulation_3d_test02a.png

all: sim_targets.mat

# CREATE NETGEN TANK MODEL
ng_cyl_mdl.mat: make_ng_cyl_mdl.m
	$(CALLMATLAB)

# SIMULATE HOMOG MOVEMENTS
sim_targets.mat: sim_targets.m ng_cyl_mdl.mat
	$(CALLMATLAB)

# IMAGE OF SIMULATED TARGETS
simulation_3d_test02a.png: simulation_3d_test02.m \
               ng_tank.mat sim_radmove_homog.mat
	$(CALLMATLAB)

# CALCULATE JACOBIAN OF TANK
GREIT_Jacobian_cyl.mat: GREIT_Jacobian_cyl.m ng_tank.mat
	$(CALLMATLAB)
