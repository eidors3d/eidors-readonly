% Sheffield MKI backprojection $Id$

% Gauss Newton Solvers
inv_mdl{1} = inv_GN;
inv_mdl{1}.hyperparameter.value= 0.2;
inv_mdl{2} = inv_GN;
inv_mdl{2}.hyperparameter.value= 0.5;
% Sheffield Backprojection solver
inv_mdl{3} = mk_common_gridmdl('backproj');

for loop=1:3
   imgr= inv_solve(inv_mdl{loop}, vh,vi);
   imgr.calc_colours.ref_level=0;
   subplot(1,3,loop); show_fem(imgr);
   axis equal; axis off
end

print_convert inv_solve_backproj04a.png
