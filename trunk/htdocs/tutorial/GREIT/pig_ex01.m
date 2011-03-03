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

show_fem(img);
print -dpng -r75 pig_ex_fmdl.png
%print_convert pig_ex_fmdl.png '-density 75';
