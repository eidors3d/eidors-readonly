% Add noise $Id$

sig= sqrt(norm(vi.meas - vh.meas));
m= size(vi.meas,1);
vi.meas = vi.meas + .01*sig*randn(m,1);
