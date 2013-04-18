% Construct the Schlumberger stimulation pattern
spacing= [1 1 1 2 3 4 6 8 8 11 12 14 17]; 		% Spacing between electrodes (usually called the "a" parameter)
multiples= [1 2 3 2 5/3 6/4 7/6 1 10/8 1 13/12 15/14 1];% Multiples (usually called the "n" parameter)
img.fwd_model.stimulation= stim_pattern_geophys( n_elec, 'Schlumberger', {'spacings', spacing,'multiples',multiples} );

% Solve the forward problem
dd  = fwd_solve(img);

% Show the pseudo-section of the apparent resistivity
show_pseudosection( img.fwd_model, dd.meas);
print_convert stimPatt_PseudoSec03.png '-density 100'
