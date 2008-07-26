% Lung images
% $Id$

load montreal_data_1995
imdl.hyperparameter.value=5e-2;
vh= zc_h_stomach_pre; % abdomen before fluid
vi= zc_stomach_0_5_60min; % each 5 minutes after drink
img= inv_solve(imdl, vh, vi);

clf; show_slices(img)
axis equal

print -r100 -dpng tutorial410b.png;
