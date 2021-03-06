function make_sim_testdata( savefile )
% Make simulated test data and save in file
% 
% (C) Andy Adler 2008 $Id$
% Licensed under GNU GPL v2 or v3

% exclude measures at electrodes
[x,y]= meshgrid(1:16,1:16); idx= abs(x-y)>1 & abs(x-y)<15;

K=1;

% SIMULATED OBJECTS
imb=  mk_common_model('c2c',16);

img= calc_jacobian_bkgnd( imb );
vv= fwd_solve( img );
test_v(K).vh= vv.meas;

img.elem_data([25,37,49:50,65:66,81:83,101:103,121:124])= 1.1;
img.elem_data([95,98:100,79,80,76,63,64,60,48,45,36,33,22])= 1.1;
vv= fwd_solve( img );
test_v(K).vi= vv.meas;

% Noise
sig = norm(test_v(1).vi - test_v(1).vh);
randn('seed',50);noise = randn(size(test_v(1).vh));

K=K+1;
test_v(K).vi = test_v(K-1).vi + noise/norm(noise)*sig * 10^(-6/20);
test_v(K).vh = test_v(K-1).vh;


params= [0.9,0.05,0.5,0.5]; %max_posn, targ_rad, z_0, z_t
stim_pat = mk_stim_patterns(16,1,'{ad}','{ad}', {'no_meas_current'}, 1);
load ng_cyl_mdl; fmdl.stimulation = stim_pat;
[vh,vi,xyzr_pt]= simulate_3d_movement(4, fmdl, params, @simulation_radmove);
%xyzr_pt= xyzr_pt([2,1,3,4],:)/15; %Change: mdl geometry at 90 deg; radius is 15
K=K+1;
test_v(K).vh = vh;
test_v(K).vi = vi(:,4);
K=K+1;
test_v(K).vh = vh;
test_v(K).vi = vi(:,3);
K=K+1;
test_v(K).vh = vh;
test_v(K).vi = vi(:,2);

save(savefile,'test_v');


% Radial Movement - $Id$  
function [xp,yp,zp]= simulation_radmove(f_frac, radius, z0,zt);
   rp= f_frac*radius; 
%  cv= 2*pi*f_frac * 73;
   cv= 0.5;
   xp= rp * cos(cv);
   yp= rp * sin(cv);
   if nargin==4; zp = mean([zt,z0]); else; zp= 0; end


