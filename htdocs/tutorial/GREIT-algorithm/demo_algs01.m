% Tank Data $Id$

% exclude measures at electrodes
[x,y]= meshgrid(1:16,1:16); idx= abs(x-y)>1 & abs(x-y)<15;

load iirc_data_2006
v(1).vh= - real(v_reference(idx,1));
v(1).vi= - real(v_rotate(idx,1));
