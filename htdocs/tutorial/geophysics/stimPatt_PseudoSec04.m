% Construct the Dipole-dipole stimulation pattern
spacing= [1 1 1 2 3 3 4 4 5 6 6 7 8 8 9 10 10 11 12 12 13 14 14 15 16 17]; % Spacing between electrodes (usually called the "a" parameter)
multiples= [1 2 3 2 1 5/3 1 2  1 1 7/6 1 1 10/8 1 1 12/10 1 1 13/12 1 1 15/14 1 1 1]; % Multiples (usually called the "n" parameter)
img.fwd_model.stimulation= stim_pattern_geophys( n_elec, 'DipoleDipole', {'spacings', spacing,'multiples',multiples} );


% Solve the forward problem
dd  = fwd_solve(img);

% Show the pseudo-section of the apparent resistivity
show_pseudosection( img.fwd_model, dd.meas);
print_convert stimPatt_PseudoSec04.png '-density 100'
