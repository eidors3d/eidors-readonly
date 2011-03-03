%% This first bit comes largely from the other tutorial
load CT3

trunk = trunk*.01;
lung  = lung*.01; lung = flipud(lung(1:3:end,:)); % need counterclockwise shapes
elec_pos = elec_pos*.01;

% Calculate electrode angles
pp= fourier_fit(trunk); sp = linspace(0,1,51);sp(end)=[]; centroid = mean(fourier_fit(pp, sp));
elec_pos = elec_pos - ones(size(elec_pos,1),1) * centroid;
electh= atan2(elec_pos(:,2),elec_pos(:,1))*180/pi; 

% Buiold a fwd model
[stim,meas_sel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'}, 1);
[fmdl, mat_idx] = ng_mk_extruded_model({2,{trunk,lung} ,[4,50],.1},[electh,1+0*electh],[0.1]);
fmdl.name = 'trunk_and_lungs';
fmdl.stimulation = stim;
fmdl.meas_select = meas_sel;
fmdl.normalize_measurements = 1;
fmdl.electrode(2:16) =  fmdl.electrode(16:-1:2); %flip electrodes to match 
fmdl.nodes = fmdl.nodes*diag([-1,-1,1]);

img = mk_image(fmdl,1);
img.elem_data( mat_idx{2} ) = 0.25;
%% Train GREIT
opt.imgsz = [64 64]; % 64-by-64 image (yes, we can do that now)
opt.distr = 3; % non-random, uniform
opt.Nsim = 500; % 500 hundred targets to train on, seems enough
opt.target_size = 0.01; %small targets
opt.target_offset = 0;
opt.noise_figure = 0.5; % this is key!
imdl=mk_GREIT_model(img, 0.25, [], opt);

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
figure
show_fem(img_ctrl);

img_in_inj= inv_solve(imdl,vh_inj,in_inj);
img_ex_inj = inv_solve(imdl,vh_inj,ex_inj);
img_inj = img_in_inj; 
img_inj.elem_data = img_ctrl.elem_data - img_ex_inj.elem_data;
figure
show_fem(img_inj);

%%
print_convert('injury','-density 60');
print_convert('control','-density 60');
