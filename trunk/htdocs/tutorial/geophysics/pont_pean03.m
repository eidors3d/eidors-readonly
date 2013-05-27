% Difference image vs simulated data
vr = data(:,9);
img = mk_image( imdl );
vh = fwd_solve(img); vh = vh.meas;

imgr = inv_solve(imdl, vh, vr);

show_fem(imgr);
print_convert pont_pean3a.png
