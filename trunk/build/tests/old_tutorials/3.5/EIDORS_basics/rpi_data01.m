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
th= linspace(0,2*pi,Nel+1)';th(1)=[];
for i=1:Nel-1;
   if i<=Nel/2;
      stim(i).stim_pattern = curr*cos(th*i);
   else;
      stim(i).stim_pattern = curr*sin(th*( i - Nel/2 ));
   end
   stim(i).meas_pattern= eye(Nel)-ones(Nel)/Nel;
end

fmdl.stimulation = stim;

show_fem(fmdl,[0,1])

print_convert('rpi_data01a.png','-density 60');
