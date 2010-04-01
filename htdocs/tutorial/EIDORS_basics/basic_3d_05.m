% Reconstruction Model $Id$
J = calc_jacobian( calc_jacobian_bkgnd( imdl) );
iRtR = inv(noser_image_prior( imdl ));
hp = 0.17;
iRN = hp^2 * speye(size(J,1));
RM = iRtR*J'/(J*iRtR*J' + iRN);
imdl.solve = @solve_use_matrix; 
imdl.solve_use_matrix.RM  = RM;

% reconstruct
imgr = inv_solve(imdl, vh, vi);

print -dpng -r125 basic_3d_04a.png
