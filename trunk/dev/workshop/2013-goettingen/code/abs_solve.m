load Rensselaer_EIT_Phantom;
vh = real(ACT2006_homog);
vi = real(ACT2000_phant);

Nel= 32; %Number of elecs
Zc = .0001; % Contact impedance
curr = 20; % applied current mA


th= linspace(0,360,Nel+1)';th(1)=[];
els = [90-th]*[1,0]; % [radius (clockwise), z=0]
elec_sz = 1/6;
fmdl= ng_mk_cyl_models([0,1,0.1],els,[elec_sz,0,0.03]);
fmdl.stimulation = mk_stim_patterns(Nel,1, '{trigccss}', '{mono}', {},30);

% Trig stim patterns
th= linspace(0,2*pi,Nel+1)';th(1)=[];
for i=1:Nel-1;
   if i<=Nel/2;
      stim(i).stim_pattern = curr*cos(th*i);
   else;
      stim(i).stim_pattern = curr*sin(th*( i - Nel/2 ));
   end
   stim(i).meas_pattern= eye(Nel)-ones(Nel)/Nel;
   stim(i).stimulation = 'Amp';
end
fmdl.stimulation = stim;

for i=1:Nel
   fmdl.electrode(i).z_contact= Zc;
end


subplot(211);
imdl = select_imdl(fmdl, {'Basic GN dif'});
imdl.hyperparameter.value = 0.1;
img = inv_solve(imdl, vh, vi);
show_fem(img)
%%%%%%%%%%%% ABS SOLVERS %%%%%%%%%%%

subplot(212);
imdl = select_imdl(fmdl, {'Basic GN abs'});
%imdl.parameters.show_iterations = 1;
imdl.hyperparameter.value = 5;
imdl.RtR_prior = @prior_noser;
img = inv_solve(imdl, vi);
show_fem(img)
eidors_colourbar(img);

