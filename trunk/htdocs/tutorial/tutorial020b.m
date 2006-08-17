% Create stimulation patterns and solve fwd_model
% $Id: tutorial020b.m,v 1.3 2006-08-17 22:11:32 aadler Exp $


% Define stimulation patterns
for i=1:3
    r_mdl.stimulation(i).stimulation= 'Amps';
    r_mdl.stimulation(i).stim_pattern= ( 0.001*i );
    r_mdl.stimulation(i).meas_pattern= 1; % measure electrode 1
end

r_mdl.solve=      @tutorial020_f_solve;

% Define an 'image'
resistor = eidors_obj('image', 'resistor');
resistor.elem_data= 1000;
resistor.fwd_model= r_mdl;

% Calculate data for 1k resistor
data_1k0 =fwd_solve( resistor );

% Now change resistor to be 1.2k
resistor.elem_data= 1200;
data_1k2 =fwd_solve( resistor );
