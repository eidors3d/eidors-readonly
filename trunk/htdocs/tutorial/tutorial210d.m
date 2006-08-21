% Reconstruct images with cheating Tikhonov prior (small model)
% $Id: tutorial210d.m,v 1.2 2006-08-21 20:46:23 aadler Exp $

smdl= mk_common_model('b2c');
smdl.RtR_prior= @tutorial210_cheat_tikhonov;
smdl.tutorial210_cheat_tikhonov.cheat_weight= .5;

% Normal Tikhonov prior
smdl.tutorial210_cheat_tikhonov.cheat_elements= [];
im_st(1)= inv_solve(smdl, vi, vh);

% Normal Tikhonov with sad face
smdl.tutorial210_cheat_tikhonov.cheat_elements=  ...
    [small_face.eyes, small_face.sad];
im_st(2)= inv_solve(smdl, vi, vh);

% Normal Tikhonov with hasmall_facey face
smdl.tutorial210_cheat_tikhonov.cheat_elements=  ...
    [small_face.eyes, small_face.smile];
im_st(3)= inv_solve(smdl, vi, vh);

% Normal Tikhonov with half face
smdl.tutorial210_cheat_tikhonov.cheat_elements=  ...
    [small_face.eyes, small_face.rsmile, small_face.lsad];
im_st(4)= inv_solve(smdl, vi, vh);

show_slices( im_st, levels);
