% Reconstruct images with cheating Laplace prior (small model)
% $Id: tutorial210f.m,v 1.1 2007-06-15 18:12:07 aadler Exp $

smdl= mk_common_model('b2c');
smdl.RtR_prior= @tutorial210_cheat_laplace;
smdl.tutorial210_cheat_laplace.cheat_weight= .3;

% Normal Tikhonov prior
smdl.tutorial210_cheat_laplace.cheat_elements= [];
im_st(1)= inv_solve(smdl, vi, vh);

% Normal Tikhonov with sad face
smdl.tutorial210_cheat_laplace.cheat_elements=  ...
    [small_face.eyes, small_face.sad];
im_st(2)= inv_solve(smdl, vi, vh);

% Normal Tikhonov with hasmall_facey face
smdl.tutorial210_cheat_laplace.cheat_elements=  ...
    [small_face.eyes, small_face.smile];
im_st(3)= inv_solve(smdl, vi, vh);

% Normal Tikhonov with half face
smdl.tutorial210_cheat_laplace.cheat_elements=  ...
    [small_face.eyes, small_face.rsmile, small_face.lsad];
im_st(4)= inv_solve(smdl, vi, vh);

show_slices( im_st, levels);
