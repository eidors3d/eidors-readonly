% Lung images
% $Id: tutorial410b.m,v 1.1 2008-03-13 21:34:55 aadler Exp $

load montreal_data_1995
imdl.hyperparameter.value=5e-2;
vh= zc_h_stomach_pre; % abdomen before fluid
vi= zc_stomach_0_5_60min; % each 5 minutes after drink
img= inv_solve(imdl_d, vh, vi);

clf; show_slices(img)
axis equal

print -r100 -dpng tutorial410b.png;
