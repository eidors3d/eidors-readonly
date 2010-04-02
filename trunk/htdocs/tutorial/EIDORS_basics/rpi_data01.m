% RPI tank model $Id$

Nel= 32; %Number of elecs
Zc = .0001; % Contact impedance
curr = 20; % applied current mA


th= linspace(0,360,Nel+1)';th(1)=[];
els = [90-th]*[1,0]; % [radius (clockwise), z=0]
elec_sz = 1/6;
fmdl= ng_mk_cyl_models([0,1,0.1],els,[elec_sz,0,0.03]);

for i=1:Nel
   fmdl.electrode(i).z_contact= Zc;
end

% Trig stim patterns
stim = mk_stim_patterns(Nel,1,'{trigccss}','{mono}', ...
       {'meas_current', 'no_balance_meas'},curr);
fmdl.stimulation = stim;

subplot(221); show_fem(fmdl,[0,1])

%print -dpng -r125 rpi_data01a.png
print -depsc2  jnk.eps;!LD_LIBRARY_PATH="" convert -density 125 jnk.eps rpi_data01a.png
