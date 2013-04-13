% create 3D forward model
fmdl = ng_mk_cyl_models([2,1,0.08],[8,0.8,1.2],[0.05]); 
fmdl.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
imgs= mk_image( fmdl, 1);

show_fem(imgs);

r = 0.05; % target radius
Npos = 20; % number of positions
Xpos =  linspace(0,0.9,Npos); % positions to simulated along x-axis
Ypos = zeros(1,Npos); 
Zpos = ones(1,Npos);  %% for off-plane, adjust the level (*1.5)
xyzr = [Xpos; Ypos; Zpos; r*ones(1,Npos)];

[vh,vi] = simulate_movement(imgs, xyzr);

% 3D inverse model
imdl= mk_common_model('n3r2',32);
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

%% calculate the GREIT parameters
levels =[inf,inf,Zpos];
show_slices(imgr, levels);
imgr.calc_colours.npoints = 128;
imgr.calc_slices.levels=levels;
params = eval_GREIT_fig_merit(imgr, xyzr);

figure
p_names = {'AR','PE','RES','SD','RNG'};
for i=1:5; subplot(5,1,i);
    plot(params(i,:)); ylabel(p_names{i});
end
