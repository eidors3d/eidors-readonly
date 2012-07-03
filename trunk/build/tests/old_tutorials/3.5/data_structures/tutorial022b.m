% 2D image $Id$

img= mk_image(mdl, conduc);
fsol= fwd_solve(img)
% fsol is the voltage drop across the resistor
