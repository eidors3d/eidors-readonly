% Add noise $Id: TV_hyperparams05.m,v 1.1 2008-03-15 00:12:20 aadler Exp $

sig= sqrt(norm(vi.meas - vh.meas));
m= size(vi.meas,1);
vi.meas = vi.meas + .01*sig*randn(m,1);
