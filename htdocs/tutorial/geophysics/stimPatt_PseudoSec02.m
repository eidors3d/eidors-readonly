% Construct the Wenner stimulation pattern
img.fwd_model.stimulation= stim_pattern_geophys( n_elec, 'Wenner', {'spacings', 1:32} );

% Use apparent_resistivity
img.fwd_model.jacobian = @jacobian_apparent_resistivity;
img.fwd_model.solve    = @fwd_solve_apparent_resistivity;

% Solve the forward problem
dd  = fwd_solve(img);

% Show the pseudo-section of the apparent resistivity
show_pseudosection( img.fwd_model, dd.meas);
print_convert stimPatt_PseudoSec02.png '-density 100'
