% This started from the EIDORS tutorial, then modified as required
% Basic 3d model $Id: basic_3d_01.m 2161 2010-04-04 20:33:46Z aadler $

% Make sure we're only doing this in a single-threaded fashion
maxNumCompThreads(1);

% save meshes for use in NDRM or elsewhere
MESHSAVE=1;
% save system matrices, right-hand sides and solutions for use in Meagre-Crowd
MMSAVE=1;
% show figures of the generated FEMs (very slow for finer meshes)
FIGURES=0;

% change mesh density
D=[ 0.25 0.1 0.08 0.06 0.04 0.02 0.01 ];
% WARNING: This takes a lot of memory to get to 0.01, try the left-most ones first.
% D=[0.01]; % this will fail: out of memory @ 64GB

for DIMENSIONS=[2 3]
  for fine_density=D
    
    % set coarse mesh density as a function of fine mesh density
    coarse_density = fine_density * 4;
    if coarse_density > 1.0
      coarse_denisty = 1.0;
    end

    fprintf('restart: fine_density=%f (%f..%f), coarse_density=%f\n',fine_density,D(1),D(end), coarse_density);
    % FORWARD
    
    tic;
    % fine mesh
    if DIMENSIONS == 3
      height=3;
      electrodes = [15,1,1.5,2];
    else
      height=0;
      electrodes = [15];
    end
    fmdl= ng_mk_cyl_models([height,1,fine_density],electrodes,[0.1,0,fine_density/5]); 
    t_f=toc;
    % coarse mesh
    tic;
    cmdl= ng_mk_cyl_models([height,1,coarse_density],electrodes,[0.1,0,coarse_density/5]); 
    t_c=toc;

    if MESHSAVE == 1
      fprintf('saving meshes...');
      save(sprintf('../data/mesh-%dD-%g-%g.mat',DIMENSIONS,fine_density,coarse_density),'fmdl','cmdl');
      fprintf('done\n');
    end

    % record mesh info
    fmdl
    cmdl
    n_nodes_f = size(fmdl.nodes,1);
    n_elems_f = size(fmdl.elems,1);
    n_nodes_c = size(cmdl.nodes,1);
    n_elems_c = size(cmdl.elems,1);
  
    % check we got something sensible: make a figure of our FEM configuration
    if FIGURES == 1
      figure
      show_fem(fmdl);
      figure
      show_fem(cmdl);
    end

    % build stimulation/measurement protocol and initial model
    fprintf('stim/meas pattern generation\n');
    tic;
    if DIMENSIONS == 3
      fmdl.stimulation = mk_stim_patterns(45,1,[0,3],[0,1],{},1);
    else
      fmdl.stimulation = mk_stim_patterns(15,1,[0,3],[0,1],{},1);
    end
    imdl = mk_common_model('a2c2',8); % Will replace most fields
    imdl.fwd_model = fmdl;
    t2=toc;
    % map coarse to fine mesh
    fprintf('c2f mapping\n');
    tic;
    imdl.fwd_model.coarse2fine = mk_coarse_fine_mapping( fmdl, cmdl);
    t_c2f=toc;
    imdl
    clear cmdl; % free up memory

    % FORWARD SIMULATIONS w/ a TARGET
    % Add a circular object at 0.2, 0.5
    % Calculate element membership in object
    select_fcn = inline('(x-0.2).^2 + (y-0.5).^2 + (z-2).^2 < 0.3^2','x','y','z');
    memb_frac = elem_select( fmdl, select_fcn);
   
    ffmdl = imdl; 
    ffmdl.fwd_model = fmdl;
    ffmdl.fwd_model.coarse2fine = [];
    clear fmdl;
    img1 = mk_image(ffmdl);
    img2 = mk_image(img1, 1 + memb_frac );
    clear memb_frac ffmdl; % free up memory

    if FIGURES == 1
      figure
      show_fem(img1);
      figure
      img2.calc_colours.cb_shrink_move = [0.3,0.6, 0];
      show_fem(img2,1);
    end

    img2
    tic;
    vh= fwd_solve(img1); % homogeneous
    t_fwd = toc;
    vi= fwd_solve(img2); % inhomogeneous (contains the target)
    clear img1 img2;
  
    if FIGURES == 1
      figure
      plot([vh.meas, vi.meas]);
    end
    
    % time the system matrix building seperatly (this is what is under the hood)
    tic;
    iii = mk_image(imdl);
    FC=feval(imdl.fwd_model.system_mat, imdl.fwd_model, iii)
    t_sys=toc;
    
    % need to drop the ground node
    pp= aa_fwd_parameters( imdl.fwd_model );
    idx= 1:size(FC.E,1);
    idx( imdl.fwd_model.gnd_node ) = [];
    v= zeros(pp.n_node,pp.n_stim);
    tol= 1e-5;
    v(idx,:)= forward_solver( FC.E(idx,idx), pp.QQ(idx,:), tol);

    % save matrices for the sparse timing later (matrix market format)
    if(MMSAVE == 1)
      fprintf('saving matrices...');
      mmwrite(sprintf('../data/sys-%dD-A-%d-%d.mm', DIMENSIONS, size(FC.E(idx,idx),1), nnz(FC.E(idx,idx))), FC.E(idx,idx));
      mmwrite(sprintf('../data/sys-%dD-b-%d-%d.mm', DIMENSIONS, size(FC.E(idx,idx),1), nnz(FC.E(idx,idx))), pp.QQ(idx,:));
      mmwrite(sprintf('../data/sys-%dD-x-%d-%d.mm', DIMENSIONS, size(FC.E(idx,idx),1), nnz(FC.E(idx,idx))), v(idx,:));
      fprintf('done\n');
    end

    % record system matrix info
    n_nz = nnz(FC.E(idx,idx));
   
    % free up some memory
    clear iii FC pp idx v;
    
    % INVERSE
    tic;
    J = calc_jacobian( calc_jacobian_bkgnd( imdl) );
    t3 = toc;
    size(J)
    
    % REGULARIZE
    tic;
    iRtR = inv(noser_image_prior( imdl ));
    hp = 0.17;
    iRN = hp^2 * speye(size(J,1));
    RM = iRtR*J'/(J*iRtR*J' + iRN);
    imdl.solve = @solve_use_matrix; 
    imdl.solve_use_matrix.RM  = RM;
    t4 = toc;
    clear J iRTR iRN RM; % free up some memory
    
    % aa_inv_solve does the following
    %   img_bkgnd= calc_jacobian_bkgnd( inv_model );
    %   J = calc_jacobian( img_bkgnd);
    %
    %   RtR = calc_RtR_prior( inv_model );
    %   W   = calc_meas_icov( inv_model );
    %   hp  = calc_hyperparameter( inv_model );
    %
    %   RM= (J'*W*J +  hp^2*RtR)\J'*W;
    
    tic;
    imgr = inv_solve(imdl, vh, vi);
    t5= toc;
    
    if FIGURES == 1
      figure
      show_fem(imgr);
    end


    % LOGGING
    % open a file to record these results
    filename = sprintf('../data/profiling-%dD-%f.log',DIMENSIONS,fine_density);
    logf = fopen([filename],'w');
    % NOTE: first line is the one we are collating
    % nodes (fine), elements (fine), nodes (coarse), elements (coarse), meshing time, senstivity est. time (incl. FEM assembly), jacobian calcs.
    fprintf(logf,'%g, %g, %g, %g, %g, %g, %g\n',n_nodes_f,n_elems_f,n_nodes_c,n_elems_c,t_f+t_c+t_c2f,t2+t3+t_sys, t4+t5);

    % more verbose info
    fprintf(logf,'fine meshing (meshing=%g nodes, %g elements, 3D, 15x3 electrodes, density=%g): %g seconds\n',n_nodes_f,n_elems_f,fine_density,t_f);
    fprintf(logf,'coarse meshing (meshing=%g nodes, %g elements, 3D, 15x3 electrodes, density=%g): %g seconds\n',n_nodes_c,n_elems_c,coarse_density,t_c);
    fprintf(logf,'coarse-to-fine mapping %g seconds\n',t_c2f);
    fprintf(logf,'inv. model construction: %g seconds\n', t2);
    fprintf(logf,'fwd model assembly (nz=%d): %g seconds\n',n_nz,t_sys);
    fprintf(logf,'fwd sim x1: %g seconds\n',t_fwd);
    fprintf(logf,'Jacobian (sensitivity estimation): %g seconds\n', t3);
    fprintf(logf,'inv assembly: %g seconds\n', t4);
    fprintf(logf,'inv solve: %g seconds\n', t5);
    fprintf(logf,'-----------\n');
   
    fclose(logf);
   
    % feature memstats
    % profile -memory on
    clear imdl vh vi imgr; % free up memory

  end
end
