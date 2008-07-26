% Demo Data $Id$

% EXCLUDE MEASURES AT ELECTRODES
[x,y]= meshgrid(1:16,1:16);
idx= abs(x-y)>1 & abs(x-y)<15;

% LOAD SOME TEST DATA
load iirc_data_2006
v_reference= - real(v_reference(idx,:));
v_rotate   = - real(v_rotate(idx,:));
