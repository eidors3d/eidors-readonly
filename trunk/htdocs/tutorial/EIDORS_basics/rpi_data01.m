% RPI tank model $Id$

Nel= 16; %Number of elecs
Zc = .001; % Contact impedance


th= linspace(0,360,Nel+1)';th(1)=[];
els = [90-th]*[1,0]; % [radius (clockwise), z=0]
elec_sz = 1/6;
fmdl= ng_mk_cyl_models([0,1,0.1],els,[elec_sz,0,0.03]);

for i=1:Nel
   fmdl.electrode(i).z_contact= z_contact;
end

% Trig stim patterns
stim = mk_stim_patterns(Nel,1,'{trigccss}','{mono}', ...
       {'meas_current', 'no_balance_meas'},1);
fmdl.stimulation = stim;
