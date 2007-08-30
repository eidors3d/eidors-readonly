% Reconstruct images with cheating Tikhonov prior (large model)
% $Id: tutorial210e.m,v 1.2 2007-08-30 03:58:28 aadler Exp $

lmdl= mk_common_model('c2c');
lmdl.RtR_prior= @tutorial210_cheat_tikhonov;
lmdl.tutorial210_cheat_tikhonov.cheat_weight= .5;

% Normal Tikhonov prior
lmdl.tutorial210_cheat_tikhonov.cheat_elements= [];
im_st(1)= inv_solve(lmdl, vh, vi);

% Normal Tikhonov with sad face
lmdl.tutorial210_cheat_tikhonov.cheat_elements=  ...
     large_face.sad;
im_st(2)= inv_solve(lmdl, vh, vi);

% Normal Tikhonov with halarge_facey face
lmdl.tutorial210_cheat_tikhonov.cheat_elements=  ...
     large_face.happy;
im_st(3)= inv_solve(lmdl, vh, vi);

% Normal Tikhonov with half face
lmdl.tutorial210_cheat_tikhonov.cheat_elements=  ...
     large_face.halfy;
im_st(4)= inv_solve(lmdl, vh, vi);

show_slices( im_st, levels);
