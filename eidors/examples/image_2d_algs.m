% Based on the 'bubble' data from Eidors2D, use several 
% different algorithms to image it

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

n_elec= 16;
n_rings= 1;
 options = {'no_meas_current','no_rotate_meas'};

[st, ms]= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', options, 10);
inv2d= mk_common_model('c2c2',n_elec);
inv2d.fwd_model.stimulation = st;
inv2d.fwd_model.meas_select = ms;

load eidors2d_bubble.mat
d1= bubble2(1280+(-255:0));
d2= bubble1(1280+(-255:0));

img= inv_solve( inv2d, d1, d2);
show_fem(img);
