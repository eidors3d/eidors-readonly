% GREIT params eval $Id$

% Reconstruct GREIT Images
imdl_gr = mk_common_gridmdl('GREITc1');
imgc{1} = inv_solve(imdl_gr, vh, vi);

% Reconstruct backprojection Images
imdl_bp = mk_common_gridmdl('backproj');
imgc{2} = inv_solve(imdl_bp, vh, vi);

fname_base = 'GREIT_test_params04';

