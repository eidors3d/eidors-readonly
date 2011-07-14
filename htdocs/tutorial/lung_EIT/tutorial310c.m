% Lung images
% $Id$

load montreal_data_1995
imdl_d.hyperparameter.value=5e-2;
img= inv_solve(imdl_d, zc_resp(:,1), zc_resp);

clf; show_slices(img)
axis equal

print_convert tutorial310c.png;
