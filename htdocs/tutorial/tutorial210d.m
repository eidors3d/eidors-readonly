% Reconstruct images with cheating Tikhonov prior
% $Id: tutorial210d.m,v 1.1 2006-08-21 20:29:59 aadler Exp $

smdl= mk_common_model('b2c');
smdl.RtR_prior= @tutorial210_cheat_tikhonov;
smdl.tutorial210_cheat_tikhonov.cheat_elements= [];
smdl.tutorial210_cheat_tikhonov.cheat_weight= .5;

im_st(1)= inv_solve(smdl, vi, vh);
levels= [0,0,0,1,1];
show_slices( inv_solve( il_g, vi_n, vhs ), levels);
