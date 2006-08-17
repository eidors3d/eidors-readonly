% Reconstruct images
% $Id: tutorial010c.m,v 1.2 2006-08-17 20:44:02 aadler Exp $

subplot(131)
show_fem(sim_img);

%reconstruct 
rec_img= inv_solve(imdl_3d, inh_data, homg_data);

% Show reconstruction as a 3D mesh
subplot(132)
show_fem(rec_img)

subplot(133)
show_slices(rec_img,[inf,inf,2.0,1,1; ...
                     inf,inf,1.0,1,2]);

pp=get(gcf,'paperposition');
set(gcf,'paperposition',pp.*[1,1,1.4,1]);
print -r75 -dpng tutorial010c.png
set(gcf,'paperposition',pp.*[1,1,1,1]);
