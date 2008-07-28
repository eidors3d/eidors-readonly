% Radial Movement - $Id$  
function [xp,yp,zp]= simulation_radmove(f_frac, radius, z0,zt);
   rp= f_frac*radius; 
   cv= 2*pi*f_frac * 73; % move 8 times around circle
   xp= rp * cos(cv);
   yp= rp * sin(cv);
   if nargin==4; zp = mean([zt,z0]); else; zp= 0; end
