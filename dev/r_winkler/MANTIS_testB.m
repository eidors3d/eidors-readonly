% (C) 2015 Robert Winkler. License: GPL version 2 or version 3
% (re)start EIDORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(getenv('EIDORSINIT'))
  setenv('NETGENDIR','/usr/share/netgen');
  setenv('PATH', [getenv('PATH') ':/usr/share/netgen']);
  run('startup.m');
  eidors_cache('cache_size', 1*1024^3 ); % 1 GB cache

  addpath(sprintf('%s/MANTIS',pwd));                                      % MANTIS functions used in EIDORS
  setenv('EIDORSINIT','1');
else
  %eidors_cache clear
end
clear; clc; close all;
eidors_msg('log_level',0);
eidors_cache('debug_off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('===========================================\n');
fprintf('=== Making Model + Forward Computation ====\n');
fprintf('===========================================\n');

%Homogeneous and inclusion conductivity
cond_h=1.0; cond_inc=2.0;

% 3D forward model without inclusion
mdl_h= ng_mk_cyl_models([1 .7],[16,.25,.75],[0.075,0.3]);


%3D forward model with inclusion 
extra={'ball','solid ball = sphere(0.2,0.2,0.8;0.2);'};
[mdl_i,mat_idx_i]= ng_mk_cyl_models([1 .7],[16,.25,.75],[0.075,0.3],extra);
mdl_i.normalize_measurements=0;

%Stimulation patterns and add to model
mdl_i.stimulation = mk_stim_patterns(16,2,'{ad}','{ad}',{'do_redundant','balance_meas'});
mdl_h.stimulation = mdl_i.stimulation;

%Create image
img_i= mk_image(mdl_i,cond_h); img_i.elem_data(mat_idx_i{2}) = cond_inc;

%Now get "real" voltages and add noise
v_i=fwd_solve(img_i); 
rng(1); % fix seed for reproducible results
v_i_n = add_noise( 30, v_i );

fprintf('Done.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inverse solution
fprintf('===========================================\n');
fprintf('=== Default GN solver 3D ==================\n');
fprintf('===========================================\n');

imdl = mdl_h;
imdl.type = 'inv_model';
imdl.fwd_model = mdl_h;
imdl.solve = @inv_solve_abs_GN; %Default Gauss Newton solvers
imdl.reconst_type = 'absolute';
imdl.jacobian_bkgnd.value= cond_h;

imdl.inv_solve_abs_GN.show_iterations=0; %Show iteration progress
imdl.inv_solve_abs_GN.max_iterations = 20; %Number of iterations
imdl.inv_solve_abs_GN.verbose = 0;

%imdl.RtR_prior=@prior_laplace; hp = 1e-10;
imdl.RtR_prior=@prior_noser; hp = 1e-4;
imdl.hyperparameter.value = hp;

tGN=tic;
img_n = inv_solve_abs_GN(imdl, v_i_n);
tGN=toc(tGN);
fprintf('Done.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse solution using MANTIS in 3D model
fprintf('===========================================\n');
fprintf('=== MANTIS solver 3D ======================\n');
fprintf('===========================================\n');
imdl3 = mdl_h;
img3D = mk_image(imdl3,1);                                                % make "image" object
img3D.elem_data = ones(size(img3D.elem_data));                            % set arbitrary (reference) conductivity

el.Npp       = 16;                                                        % set number of electrodes per plane
el.P         =  2;                                                        % set number of plane
el.Sizes     = 0.3*0.075*ones(el.Npp*el.P,1);                             % electrode sizes
el.Impedance = num2cell(0.5*el.Sizes);                                    % arbitrary (reference) contact impedance
reco.verbose = 0;

tMANTIS=tic;
[img3R,~,~,~,~] = mantis_solve(img3D,v_i_n.meas,el,reco);
tMANTIS=toc(tMANTIS);

%% now, show the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('===========================================\n');
fprintf('=== Show results ==========================\n');
fprintf('===========================================\n');
calc_colours('defaults');
calc_colours('ref_level',cond_h);                                         % set "background" color for plot
calc_colours('clim',abs(cond_h-cond_inc));                                % set range

fg1 = figure('Color','w','Position',[20, 400, 1200, 400]); hold on
subplot(1,3,1); show_fem(img_i);title('Setting');                         % 3D FEM plot of "true" conductivity
subplot(1,3,2); show_fem(img_n);title(sprintf('GN (%2.1fs)',tGN));        % 3D FEM plot of GN reconstruction
subplot(1,3,3); show_fem(img3R);title(sprintf('MANTIS (%2.1fs)',tMANTIS));% 3D FEM plot of 3D MANTIS reconstruction
drawnow;

fg2 = figure('Color','w','Position',[20, 400, 1200, 1000]); hold on
nSlices = 3; dom.Height = 1;

subplot(1,3,1); show_slices(img_i,[inf(nSlices,2), ...                    % display some horizontal slices of "true" conductivity
  dom.Height*linspace(nSlices/(nSlices+1),1/(nSlices+1),nSlices)']);
  title('Setting');
subplot(1,3,2); show_slices(img_n,[inf(nSlices,2), ...                    % display some horizontal slices of GN reconstruction
  dom.Height*linspace(nSlices/(nSlices+1),1/(nSlices+1),nSlices)']);
  title(sprintf('GN (%2.1fs)',tGN));
subplot(1,3,3); show_slices(img3R,[inf(nSlices,2), ...                    % display some horizontal slices of 3D MANTIS reco.
  dom.Height*linspace(nSlices/(nSlices+1),1/(nSlices+1),nSlices)']);
  title(sprintf('MANTIS (%2.1fs)',tMANTIS));
drawnow;

fprintf('GN     time: %8.2fs\n',tGN);
fprintf('MANTIS time: %8.2fs\n',tMANTIS);




