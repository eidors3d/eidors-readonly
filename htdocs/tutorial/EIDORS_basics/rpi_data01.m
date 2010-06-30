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
       {'meas_current', 'balance_meas'},curr);
if 0
   Nel=32; th= linspace(0,2*pi,Nel+1)';th(1)=[];
   for i=1:Nel-1;
     stim(i).meas_pattern= eye(Nel)-ones(Nel)/Nel;
     if i<=Nel/2;
        k=i;
        stim(i).stim_pattern = cos(th*k);
     else;
        k=i-Nel/2;
        stim(i).stim_pattern = sin(th*k);
     end
   end
end
fmdl.stimulation = stim;

show_fem(fmdl,[0,1])

print_convert('rpi_data01a.png','-density 60');
