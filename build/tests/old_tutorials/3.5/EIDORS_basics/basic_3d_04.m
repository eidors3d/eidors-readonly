% Basic 3d model $Id$

% Simulate Voltages and plot them
vh= fwd_solve(img1);
vi= fwd_solve(img2);

plot([vh.meas, vi.meas]);
axis tight
print_convert('basic_3d_04a.png','-density 60',0.4);
