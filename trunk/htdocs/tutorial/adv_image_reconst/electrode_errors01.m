% Sample data + Image $Id$

% Sample data
load montreal_data_1995
vi = double( zc_resp );
vh = double( zc_resp(:,1) );

% Reconstruct
imdl = mk_common_model('c2t2',16);
imdl.hyperparameter.value = .1;

img = inv_solve(imdl, vh, vi);
subplot(221)
show_slices(img);
print_convert electrode_errors01a.jpg
close % else matlab remembers image axis
subplot(221) % to preserve size
img = inv_solve(imdl, vh, vi(:,20));
show_fem(img,[0,1,0]);axis off
print_convert electrode_errors01b.jpg
