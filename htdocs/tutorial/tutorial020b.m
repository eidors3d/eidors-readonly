% tutorial1_create_fwd_model
% $Id: tutorial020b.m,v 1.1 2006-08-17 17:57:33 aadler Exp $


% Define stimulation patterns
for i=1:3
    r_mdl.stimulation(i).stimulation= 'Amps';
    r_mdl.stimulation(i).stim_pattern= ( 0.001*i );
    r_mdl.stimulation(i).meas_pattern= 1; % measure electrode 1
end

r_mdl.solve=      @tutorial020_f_solve;

% Define an 'image'
img_1k = eidors_obj('image', 'resistor');
img_1k.elem_data= 1000;
img_1k.fwd_model= r_mdl;

data_1k =fwd_solve( img_1k );

