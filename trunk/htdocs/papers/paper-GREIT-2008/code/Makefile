# Create Simulated data for GREIT evaluations
# Usage set MATLAB 

ifndef MATLAB
   MATLAB = C:/progra~1/matlab6p5/bin/win32/matlab.exe -nosplash
   EIDORSTARTUP = run('../../../../eidors/startup.m')
   MC1= $(EIDORSTARTUP);run ./
   MC2= exit
   
endif

# CREATE NETGEN TANK MODEL
ng_tank.mat: simulation_3d_test00.m
	$(MATLAB) -r "$(EIDORSTARTUP);run ./$<;exit"
