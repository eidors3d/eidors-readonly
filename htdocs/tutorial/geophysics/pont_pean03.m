% Difference image vs simulated data
gps = load('Mine_20FEV2004.gps');
data= load('Mine_20FEV2004_LI.tomel');
vr = data(:,9);
img = mk_image( imdl );
vh = fwd_solve(img); vh = vh.meas;

imgr = inv_solve(imdl, vh, vr);

show_fem(imgr);
