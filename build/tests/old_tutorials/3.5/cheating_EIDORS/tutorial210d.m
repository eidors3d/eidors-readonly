% Reconstruct images with cheating Tikhonov prior (small model)
% $Id$

smdl= mk_common_model('b2c');
smdl.RtR_prior= @tutorial210_cheat_tikhonov;
smdl.tutorial210_cheat_tikhonov.cheat_weight= .5;

% Normal Tikhonov prior
smdl.tutorial210_cheat_tikhonov.cheat_elements= [];
im_st(1)= inv_solve(smdl, vh, vi);

% Normal Tikhonov with sad face
smdl.tutorial210_cheat_tikhonov.cheat_elements=  ...
    [small_face.eyes, small_face.sad];
im_st(2)= inv_solve(smdl, vh, vi);

% Normal Tikhonov with hasmall_facey face
smdl.tutorial210_cheat_tikhonov.cheat_elements=  ...
    [small_face.eyes, small_face.smile];
im_st(3)= inv_solve(smdl, vh, vi);

% Normal Tikhonov with half face
smdl.tutorial210_cheat_tikhonov.cheat_elements=  ...
    [small_face.eyes, small_face.rsmile, small_face.lsad];
im_st(4)= inv_solve(smdl, vh, vi);
im_st(1).calc_colours.greylev = .01;
im_st(1).show_slices.img_cols = 4;
show_slices( im_st);
