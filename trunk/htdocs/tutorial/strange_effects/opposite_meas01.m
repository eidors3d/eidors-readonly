%$Id$
imdl= mk_common_model('c2c2',16);
imdl.fwd_model = rmfield(imdl.fwd_model,'meas_select');


% Simulate with AD / AD
imdl.fwd_model.stimulation = mk_stim_patterns(16,1,'{ad}','{op}',{},1);
[vi,vh] = simulate_2d_movement( 5, imdl.fwd_model, [.9,.05], 2);
imgs= inv_solve(imdl, vh, vi); 
subplot(411); show_slices(imgs,[0,inf,inf,1,1])

% Simulate with OP / AD
imdl.fwd_model.stimulation = mk_stim_patterns(16,1,'{op}','{ad}',{},1);
[vi,vh] = simulate_2d_movement( 5, imdl.fwd_model, [.9,.05], 2);
imgs= inv_solve(imdl, vh, vi); 
subplot(412); show_slices(imgs,[0,inf,inf,1,1])

% Simulate with AD / OP
imdl.fwd_model.stimulation = mk_stim_patterns(16,1,'{ad}','{op}',{},1);
[vi,vh] = simulate_2d_movement( 5, imdl.fwd_model, [.9,.05], 2);
imgs= inv_solve(imdl, vh, vi); 
subplot(413); show_slices(imgs,[0,inf,inf,1,1])

% Simulate with OP / OP
imdl.fwd_model.stimulation = mk_stim_patterns(16,1,'{op}','{op}',{},1);
[vi,vh] = simulate_2d_movement( 5, imdl.fwd_model, [.9,.05], 2);
imgs= inv_solve(imdl, vh, vi); 
subplot(414); show_slices(imgs,[0,inf,inf,1,1])

print -dpng -r100 opposite_meas01.png
