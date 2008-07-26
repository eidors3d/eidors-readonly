% Create 2D model $Id$

nn= 84;     % number of nodes
ww=4;       % width = 4
conduc= 1;  % conductivity in Ohm-meters
mdl= eidors_obj('fwd_model','2D rectangle');
mdl.nodes = [floor( (0:nn-1)/ww );rem(0:nn-1,ww)]';
mdl.elems = delaunayn(mdl.nodes);
mdl.gnd_node = 1;
elec(1).nodes= [1:ww];      elec(1).z_contact= 0;
elec(2).nodes= nn-[0:ww-1]; elec(2).z_contact= 0;
stim.stim_pattern= [-1;1];
stim.meas_pattern= [-1,1];
mdl.stimulation= stim;
mdl.electrode= elec;
mdl.solve = @aa_fwd_solve;
mdl.system_mat = @aa_calc_system_mat;

show_fem(mdl); axis('equal'); set(gca,'Ylim',[-.5,ww-.5]);
print -r100 -dpng tutorial022a.png
