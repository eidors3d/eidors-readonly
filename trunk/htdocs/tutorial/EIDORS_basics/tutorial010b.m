% Simulate EIT data
% $Id$

sim_img= mk_image(imdl_3d.fwd_model,1);

% set voltage and current stimulation patterns
stim =  mk_stim_patterns(16,2,[0,1],[0,1],{},1);
sim_img.fwd_model.stimulation = stim;

% set homogeneous conductivity and simulate
homg_data=fwd_solve( sim_img );

% set inhomogeneous conductivity and simulate
sim_img.elem_data([390,391,393,396,402,478,479,480,484,486, ...
                   664,665,666,667,668,670,671,672,676,677, ...
                   678,755,760,761])= 1.15;
sim_img.elem_data([318,319,321,324,330,439,440,441,445,447, ...
                   592,593,594,595,596,598,599,600,604,605, ...
                   606,716,721,722])= 0.8;
inh_data=fwd_solve( sim_img );

clf;subplot(211);

xax= 1:length(homg_data.meas);
hh= plotyy(xax,[homg_data.meas, inh_data.meas], ...
           xax, homg_data.meas- inh_data.meas );

set(hh,'Xlim',[1,max(xax)]);
print_convert('tutorial010b.png','-density 75');
