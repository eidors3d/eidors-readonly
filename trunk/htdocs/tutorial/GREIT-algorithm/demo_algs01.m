% Tank Data $Id:$

% exclude measures at electrodes
[x,y]= meshgrid(1:16,1:16); idx= abs(x-y)>1 & abs(x-y)<15;

load iirc_data_2006
vh1= - real(v_reference(idx,1));
vi1= - real(v_rotate(idx,1));
