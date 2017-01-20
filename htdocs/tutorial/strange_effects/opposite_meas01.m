%$Id$
imdl= mk_common_model('c2c2',16);
imdl.fwd_model = rmfield(imdl.fwd_model,'meas_select');


% Simulate with AD / AD
imdl.fwd_model.stimulation = mk_stim_patterns(16,1,'{ad}','{ad}',{},1);
[vi,vh] = simulate_2d_movement( 5, imdl.fwd_model, [.9,.05], 2);
imgs= inv_solve(imdl, vh, vi);
imgs.show_slices.img_cols = 5; imgs.calc_colours.ref_level = 0;
subplot(411); show_slices(imgs)

% Simulate with OP / AD
imdl.fwd_model.stimulation = mk_stim_patterns(16,1,'{op}','{ad}',{},1);
[vi,vh] = simulate_2d_movement( 5, imdl.fwd_model, [.9,.05], 2);
imgs= inv_solve(imdl, vh, vi); 
imgs.show_slices.img_cols = 5; imgs.calc_colours.ref_level = 0;
subplot(412); show_slices(imgs)

% Simulate with AD / OP
imdl.fwd_model.stimulation = mk_stim_patterns(16,1,'{ad}','{op}',{},1);
[vi,vh] = simulate_2d_movement( 5, imdl.fwd_model, [.9,.05], 2);
imgs= inv_solve(imdl, vh, vi); 
imgs.show_slices.img_cols = 5; imgs.calc_colours.ref_level = 0;
subplot(413); show_slices(imgs)

% Simulate with OP / OP
imdl.fwd_model.stimulation = mk_stim_patterns(16,1,'{op}','{op}',{},1);
[vi,vh] = simulate_2d_movement( 5, imdl.fwd_model, [.9,.05], 2);
imgs= inv_solve(imdl, vh, vi); 
imgs.show_slices.img_cols = 5; imgs.calc_colours.ref_level = 0;
subplot(414); show_slices(imgs)

print_convert opposite_meas01.png
