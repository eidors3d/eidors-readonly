% 3D inverse model
imdl= mk_common_model('n3r2',[16,2]);
mdl= ng_mk_cyl_models([2,1,0.1],[8,0.8,1.2],[0.05]); 
mdl.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
imdl.fwd_model = mdl; 

img= mk_image( mdl, 1); J = calc_jacobian(img);

%% inverse solution (faster solution)
hp  = 0.015;
RtR = prior_noser( imdl ); P= inv(RtR);
Rn = speye( size(J,1) );
imdl.solve = @solve_use_matrix;
imdl.solve_use_matrix.RM = P*J'/(J*P*J' + hp^2*Rn);
imgr= inv_solve(imdl,vh,vi);
