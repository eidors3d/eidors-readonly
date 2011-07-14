%% Read in the data
ctrl = eidors_readdata('2-control.RAW');
inj  = eidors_readdata('2-injury.RAW');

ex_ctrl = ctrl(:,101);
in_ctrl = ctrl(:,103);

ex_inj = inj(:,99);
in_inj = inj(:,101);

%% Reconstruct
img_ctrl = inv_solve(imdl, ex_ctrl, in_ctrl);
show_fem(img_ctrl); axis off
print_convert pig_control.png '-density 60'


img_inj = inv_solve(imdl, ex_inj, in_inj);
show_fem(img_inj); axis off
print_convert pig_injury.png '-density 60'
