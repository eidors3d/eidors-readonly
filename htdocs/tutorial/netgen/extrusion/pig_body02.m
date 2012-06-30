% Calculate electrode angles
pp= fourier_fit(trunk); sp = linspace(0,1,51);sp(end)=[]; centroid = mean(fourier_fit(pp, sp));
elec_pos = elec_pos - ones(size(elec_pos,1),1) * centroid;
electh= atan2(elec_pos(:,2),elec_pos(:,1))*180/pi; % electh=electh(1:3);

fmdl = ng_mk_extruded_model({2,{trunk, lung} ,[4,50],.1},[electh,1+0*electh],[0.1]);
[stim,meas_sel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'}, 1);
fmdl.stimulation = stim;

fmdl.nodes = fmdl.nodes*diag([-1,-1,1]); % Flip x,y axis to match medical direction
img=mk_image(fmdl,1);
img.elem_data(fmdl.mat_idx{2})= 0.3;

clf; show_fem(img); view(0,70);
print_convert pig_body02.jpg
