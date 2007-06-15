% Lung images
% $Id: tutorial310c.m,v 1.1 2007-06-15 18:10:20 aadler Exp $

load montreal_data_1995
imdl_d.hyperparameter.value=2;
img= inv_solve(imdl_d, zc_resp(:,22), zc_resp);

clf; show_slices(img)
axis equal

print -r100 -dpng tutorial310c.png;
