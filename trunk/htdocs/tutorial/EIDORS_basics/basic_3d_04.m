% Basic 3d model $Id$

% Simulate Voltages and plot them
vh= fwd_solve(img1);
vi= fwd_solve(img2);

plot([vh.meas, vi.meas]);
print -dpng -r125 basic_3d_04a.png
