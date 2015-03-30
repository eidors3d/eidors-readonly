shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
             'solid block  = orthobrick(-4,-4,-2;4,4,2) -maxh=0.3;\n' ...
             'solid mainobj= top and block;\n'];
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
fmdl.stimulation = stim;

% inverse model: a smaller region directly over the electrodes with coarser mesh
shape_str = 'solid mainobj= orthobrick(-3,-3,-2.0;3,3,-0.2) -maxh=0.5;';
cmdl = ng_mk_gen_models(shape_str, [], [], '');

% define the mapping between the two meshes
c2f= mk_coarse_fine_mapping( fmdl, cmdl);

% set the reconstruction parameters
inv3d= select_imdl(fmdl, {'Basic GN dif'});
inv3d.solve= @inv_solve_diff_GN_one_step;
inv3d.hyperparameter.value = .03;
inv3d.RtR_prior= @prior_laplace;
inv3d.fwd_model.coarse2fine = c2f;
inv3d.rec_model = cmdl;

hh= show_fem(cmdl); set(hh,'EdgeColor',[0,0,1]);
hold on; show_fem(fmdl); hold off
view(-16,22); set(gca,'Projection','perspective')
print_convert surface_sim03a.png '-density 100'
