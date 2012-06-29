if 0
shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
             'solid block  = orthobrick(-4,-4,-2;4,4,2) -maxh=0.3;\n' ...
             'solid mainobj= top and block;\n'];
[fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
fmdl.stimulation = mk_stim_patterns(size(fmdl.electrode,2), 1, ... %rings
    [0,7], [0,7], {'no_meas_current'},1);
end

% inverse model: a smaller region directly over the electrodes with coarser mesh
shape_str = 'solid mainobj= orthobrick(-5,-5,-2.0;5,5,-0.2) -maxh=0.5;';
cmdl = ng_mk_gen_models(shape_str, [], [], '');

% define the mapping between the two meshes
c2f= mk_coarse_fine_mapping( fmdl, cmdl);

% set the reconstruction parameters
inv3d= select_imdl(fmdl, {'Basic GN dif'});
inv3d.solve= @inv_solve_diff_GN_one_step;
inv3d.hyperparameter.value = 1e-1;
inv3d.RtR_prior= @laplace_image_prior;
inv3d.RtR_prior= @prior_noser;
inv3d.fwd_model.coarse2fine = c2f;
inv3d.rec_model = cmdl;
