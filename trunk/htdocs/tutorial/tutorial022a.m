% Create 2D model $Id: tutorial022a.m,v 1.1 2007-03-26 17:49:05 aadler Exp $

nn= 84;     % number of nodes
ww=4;       % width = 4
conduc= 1;  % conductivity in Ohm-meters
mdl= eidors_obj('fwd_model','2D rectangle');
mdl.nodes = [floor( (0:nn-1)/ww );rem(0:nn-1,ww)]';
mdl.elems = delaunayn(mdl.nodes);
mdl.gnd_node = 1;
stim.stim_pattern= zeros(nn,1);
stim.stim_pattern(1:ww) = -1; % 1 amp out gnd
stim.stim_pattern(nn-(0:ww-1)) = 1; % 1 amp in end
stim.meas_pattern= zeros(1,nn);
stim.meas_pattern([1,nn])= [1,-1];
mdl.stimulation(1)= stim;
    img= eidors_obj('image','2D rectangle', ...
          'elem_data', ones(size(mdl.elems,1),1) * conduc, ...
          'fwd_model', mdl); 

show_fem(mdl); axis('equal'); set(gca,'Ylim',[-.5,ww-.5]);
print -r100 -dpng tutorial022a.png
