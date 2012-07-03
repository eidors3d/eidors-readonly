% Reject electrodes + Image $Id$


imdl.meas_icov = meas_icov_rm_elecs( imdl, 13);

img = inv_solve(imdl, vh, vi(:,20));
show_fem(img,[0,1,0]);
print_convert electrode_errors02a.jpg

imdl.meas_icov = meas_icov_rm_elecs( imdl, [13,5]);

img = inv_solve(imdl, vh, vi(:,20));
show_fem(img,[0,1,0]);
print_convert electrode_errors02b.jpg
