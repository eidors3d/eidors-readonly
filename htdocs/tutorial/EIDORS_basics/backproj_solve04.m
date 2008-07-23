% Sheffield MKI backprojection $Id: backproj_solve04.m,v 1.1 2008-07-23 14:27:56 aadler Exp $

% Gauss Newton Solvers
inv_mdl{1} = inv_GN;
inv_mdl{1}.hyperparameter.value= 0.03;
inv_mdl{2} = inv_GN;
inv_mdl{2}.hyperparameter.value= 0.30;
% Sheffield Backprojection solver
inv_mdl{3} = mk_common_gridmdl('backproj');

for loop=1:3
   imgr= inv_solve(inv_mdl{loop}, vh,vi);
   imgr.calc_colours.ref_level=0;
   subplot(1,3,loop); show_fem(imgr);
   axis equal; axis off
end

print -r125 -dpng backproj_solve04a.png;
