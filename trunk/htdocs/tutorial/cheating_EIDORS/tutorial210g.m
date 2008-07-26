% Reconstruct images with cheating Laplace prior (large model)
% $Id$

lmdl= mk_common_model('c2c');
lmdl.RtR_prior= @tutorial210_cheat_laplace;
lmdl.tutorial210_cheat_laplace.cheat_weight= .3;

% Normal Tikhonov prior
lmdl.tutorial210_cheat_laplace.cheat_elements= [];
im_st(1)= inv_solve(lmdl, vh, vi);

% Normal Tikhonov with sad face
lmdl.tutorial210_cheat_laplace.cheat_elements=  ...
     large_face.sad;
im_st(2)= inv_solve(lmdl, vh, vi);

% Normal Tikhonov with halarge_facey face
lmdl.tutorial210_cheat_laplace.cheat_elements=  ...
     large_face.happy;
im_st(3)= inv_solve(lmdl, vh, vi);

% Normal Tikhonov with half face
lmdl.tutorial210_cheat_laplace.cheat_elements=  ...
     large_face.halfy;
im_st(4)= inv_solve(lmdl, vh, vi);

show_slices( im_st, levels);
