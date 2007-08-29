function ok= calc_jacobian_test
% Verify Jacobian Calculation by small derivative from forward problem
% Also calculate dataprior
%     Difference dataprior should be 1
%     normalized difference dataprior should be 1./ homg_data

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_jacobian_test.m,v 1.23 2007-08-29 09:25:18 aadler Exp $

ok= 1;
delta = 1e-4;
testvec= [5,20,40,130];

mdl= make_aa_mdl2;
%
disp('test aa_calc_jacobian (2D) for difference data')
ok=ok & run_jacobian_test( mdl, delta, testvec );
ok=ok & run_dataprior_test( mdl );

%
disp('test aa_calc_jacobian (2D) for normalized difference data');
mdl.normalize_measurements= 1;
ok=ok & run_jacobian_test( mdl, delta, testvec );
ok=ok & run_dataprior_test( mdl );

mdl= make_aa_mdl3;
%
disp('test aa_calc_jacobian (3D) for difference data')
ok=ok & run_jacobian_test( mdl, delta, testvec );
ok=ok & run_dataprior_test( mdl );

mdl.normalize_measurements= 1;
disp('test aa_calc_jacobian (3D) for normalized difference data')
ok=ok & run_jacobian_test( mdl, delta, testvec );
ok=ok & run_dataprior_test( mdl );

mdl= make_np_mdl;
%
disp('test np_calc_jacobian for difference data')
ok=ok & run_jacobian_test( mdl, delta, testvec );
ok=ok & run_dataprior_test( mdl );

mdl.normalize_measurements= 1;
disp('test np_calc_jacobian for normalized difference data')
ok=ok & run_jacobian_test( mdl, delta, testvec );
ok=ok & run_dataprior_test( mdl );


% run the jacobian test 
function ok= run_jacobian_test( mdl, delta, testvec ); 
    calc_norm = 0;
    if isfield(mdl,'normalize_measurements')
       if mdl.normalize_measurements
           calc_norm = 1;
       end
    end

    img= eidors_obj('image', 'homg image');
    img.fwd_model= mdl;

    img.elem_data= ones( size(mdl.elems,1) ,1);
    homg_data=fwd_solve( img);

    J= calc_jacobian( img );

    % J = dF/dx = [F(x+d)  - F(x)]/d
    sumdiff= 0;
    bkgnd_elem_data= img.elem_data;
    for testelem = testvec
       img.elem_data= bkgnd_elem_data;
       img.elem_data(testelem)= bkgnd_elem_data(testelem)+delta;
       inh_data=fwd_solve( img);

       if calc_norm
          simJ= 1/delta* (inh_data.meas ./ homg_data.meas - 1);
       else
          simJ= 1/delta* (inh_data.meas-homg_data.meas);
       end
       
       plot([J(:,testelem) simJ]);
       sumdiff = sumdiff + std( J(:,testelem) - simJ );
    end

    tol= 1e-4*std(J(:));
    dif= sumdiff/length(testvec);
    if sumdiff/length(testvec) > tol
       eidors_msg('Jacobian calculation error (%s) tol(%3.2g)>diff(%3.2g)', ...
                mdl.name,tol,dif, 1);
       ok=0;
    else
       ok=1;
    end

%
% test dataprior
%     Difference dataprior should be 1
%     normalized difference dataprior should be 1./ homg_data
function ok= run_dataprior_test( mdl )
    img= eidors_obj('image', 'homg image');
    img.fwd_model= mdl;

    img.elem_data= ones( size(mdl.elems,1) ,1);
    homg_data=fwd_solve( img);

    DP= calc_meas_icov( img );

    % difference dataprior
    testvec= diag(DP);
    if isfield(mdl,'normalize_measurements')
       if mdl.normalize_measurements
           testvec = homg_data.meas.^2 .* diag(DP);
       end
    end

    if max(abs(diff( testvec ))) > 1e-12 
       ok=0;
       eidors_msg('Dataprior calculation error (%s)', mdl.name, 1);
    else
       ok=1;
    end



% 2D model with point electrodes
% 
function mdl= make_aa_mdl2;
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
    mdl = eidors_obj('fwd_model', params);
    mdl.name= 'AA_1996 mdl';

function mdl= make_aa_mdl3;
    i_mdl = mk_common_model('b3cz',16);
    mdl= i_mdl.fwd_model;
    mdl.name= 'AA_1996 mdl';
    mdl.solve=      'aa_fwd_solve';
    mdl.system_mat= 'aa_calc_system_mat';
    mdl.jacobian=   'aa_calc_jacobian';
    
    

function mdl= make_np_mdl;
    i_mdl = mk_common_model('n3r2',16);
    mdl= i_mdl.fwd_model;
    mdl.name=     'NP_2003 mdl';
    mdl.solve=      'np_fwd_solve';
    mdl.system_mat= 'np_calc_system_mat';
    mdl.jacobian=   'np_calc_jacobian';

function mdl= make_ms_mdl;
    i_mdl = mk_common_model('n3r2',16);
    mdl= i_mdl.fwd_model;
    mdl.name=       'MS_2005 mdl';
    mdl.solve=      'np_fwd_solve';
    mdl.system_mat= 'np_calc_system_mat';
    mdl.jacobian=   'ms_calc_jacobian';
