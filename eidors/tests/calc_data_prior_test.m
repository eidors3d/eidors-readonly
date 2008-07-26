function ok= calc_data_prior_test
% Verify dataprior:
%   for  difference EIT: dataprior should be 1
%   for  normalized EIT: dataprior should be 1 / homogeneous
% $Id$

ok= 1;

% Create simple 16 electrode 2D model
% 
n_elec= 16;
n_rings= 1;
 options = {'no_meas_current','no_rotate_meas'};
params= mk_circ_tank(8, [], n_elec);

params.stimulation= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', ...
                            options, 10);
params.solve=      'aa_fwd_solve';
params.system_mat= 'aa_calc_system_mat';
params.jacobian=   'aa_calc_jacobian';
params.normalize_measurements = 0;
mdl_2d = eidors_obj('fwd_model', params);

% create homogeneous image + simulate data
mat= ones( size(mdl_2d.elems,1) ,1);

homg_img= eidors_obj('image', 'homg image', ...
                     'elem_data', mat, 'fwd_model', mdl_2d );

homg_data=fwd_solve( homg_img);
J= calc_jacobian( homg_img );

sumdiff= 0;
delta = 2e-5;
testvec= [5,20,40,130];
for testelem = testvec
    mat= ones( size(mdl_2d.elems,1) ,1);
    mat(testelem)= 1+delta;
    inh_img= eidors_obj('image', 'inh', ...
                         'elem_data', mat, 'fwd_model', mdl_2d );
    inh_data=fwd_solve( inh_img);

    simJ= 1/delta* (inh_data.meas - homg_data.meas);
    
%   plot([J(:,testelem) simJ]);
    sumdiff = sumdiff + std( J(:,testelem) - simJ );
end

tol= 1e-4*std(J(:));
if sumdiff/length(testvec) > tol
   error('Jacobian calculation error');
   ok=0;
end


% create normalized model
params.normalize_measurements = 1;
mdl_2d_norm = eidors_obj('fwd_model', params);
% create homogeneous image with normalize_measurements
mat= ones( size(mdl_2d.elems,1) ,1);
homg_img= eidors_obj('image', 'homg image', ...
                     'elem_data', mat, 'fwd_model', mdl_2d_norm );

J= calc_jacobian( homg_img );
sumdiff= 0;
delta = 2e-5;
testvec= [5,20,40,130];
for testelem = testvec
    mat= ones( size(mdl_2d.elems,1) ,1);
    mat(testelem)= 1+delta;
    inh_img= eidors_obj('image', 'inh', ...
                         'elem_data', mat, 'fwd_model', mdl_2d );
    inh_data=fwd_solve( inh_img);

    simJ= 1/delta* (inh_data.meas ./ homg_data.meas - 1);
    
%   plot([J(:,testelem) simJ]);
    sumdiff = sumdiff + std( J(:,testelem) - simJ );
end

tol= 1e-4*std(J(:));
if sumdiff/length(testvec) > tol
   error('normalize Jacobian calculation error');
   ok=0;
end
