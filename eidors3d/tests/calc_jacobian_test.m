function ok= calc_jacobian_test
% Verify Jacobian Calculation by small derivative from forward problem
% Also calculate dataprior
%     Difference dataprior should be 1
%     normalized difference dataprior should be 1./ homg_data

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: calc_jacobian_test.m,v 1.8 2005-11-30 17:46:46 aadler Exp $

ok= 1;

mdl= make_aa_mdl;
%
%
% test Jacobian
delta = 2e-5;
testvec= [5,20,40,130];
ok=ok & run_jacobian_test(mdl, delta, testvec );

%
%
% test dataprior
DP= calc_data_prior( homg_img );
testvec= diag(DP);
if max(abs(diff( testvec ))) > 1e-12 
   ok=0;
   error('Dataprior calculation error');
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

    simJ= 1/delta* (1 - inh_data.meas ./ homg_data.meas);
    
%   plot([J(:,testelem) simJ]);
    sumdiff = sumdiff + std( J(:,testelem) - simJ );
end

tol= 1e-4*std(J(:));
if sumdiff/length(testvec) > tol
   ok=0;
   error('normalize Jacobian calculation error');
end

% test dataprior
DP= calc_data_prior( homg_img );

testvec = homg_data.meas.^2 .* diag(DP);
if max(abs(diff( testvec ))) > 1e-12 
   max(abs(diff( testvec )))
   ok=0;
   error('Dataprior calculation error');
end


% run the jacobian test 
function ok= run_jacobian_test( mdl, delta, testvec ); 
    img= eidors_obj('image', 'homg image');
    img.fwd_model= mdl;

    img.elem_data= ones( size(mdl.elems,1) ,1);
    homg_data=fwd_solve( img);

    J= calc_jacobian( img );

    sumdiff= 0;
    for testelem = testvec
        mat= ones( size(mdl.elems,1) ,1);
        mat(testelem)= 1+delta;
        img.elem_data= mat;
        inh_data=fwd_solve( img);

        simJ= 1/delta* (homg_data.meas-inh_data.meas);
        
    %   plot([J(:,testelem) simJ]);
        sumdiff = sumdiff + std( J(:,testelem) - simJ );
    end

    tol= 1e-4*std(J(:));
    if sumdiff/length(testvec) > tol
       ok=0;
       error('Jacobian calculation error');
    else
       ok=1;
    end

% Create simple 2D model to test aa_calc_jacobian
% 
function mdl_2d= make_aa_mdl;
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


