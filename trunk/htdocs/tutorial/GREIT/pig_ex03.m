%% Read in the data
ctrl = eidors_readdata('2-control.RAW');
inj  = eidors_readdata('2-injury.RAW');

vh_ctrl = mean(ctrl(:,80:120),2);
ex_ctrl = ctrl(:,101);
in_ctrl = ctrl(:,103);

vh_inj = mean(inj(:,80:120),2);
ex_inj = inj(:,99);
in_inj = inj(:,101);

%% Reconstruct
img_in_ctrl = inv_solve(imdl,vh_ctrl,in_ctrl);
img_ex_ctrl = inv_solve(imdl,vh_ctrl,ex_ctrl);
img_ctrl = img_in_ctrl; 
img_ctrl.elem_data = img_ctrl.elem_data - img_ex_ctrl.elem_data;
show_fem(img_ctrl);
print_convert('pig_ex_control.png','-density 60');


img_in_inj= inv_solve(imdl,vh_inj,in_inj);
img_ex_inj = inv_solve(imdl,vh_inj,ex_inj);
img_inj = img_in_inj; 
img_inj.elem_data = img_ctrl.elem_data - img_ex_inj.elem_data;
show_fem(img_inj);
print_convert('pig_ex_injury.png','-density 60');

%%
