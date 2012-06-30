stim = mk_stim_patterns(size(fmdl.electrode,2), 1, ... %rings
    [0,7], [0,7], {'no_meas_current'},1);

shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
             'solid block  = orthobrick(-4,-4,-2;4,4,2) -maxh=0.3;\n' ...
             'solid ball   = sphere(1,-1,-1;0.2); tlo ball;\n' ...
             'solid mainobj= top and block and not ball;\n'];
[epos_x,epos_y] = meshgrid(linspace( -2,2,5),linspace(-2,2,5));
elec_pos = [epos_x(:), epos_y(:), ones(size(epos_x(:)))*[0,0,0,1] ];
elec_shape=[0.1];
elec_obj = 'top';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
fmdl.stimulation = stim;
img= mk_image(fmdl,1);
img.elem_data(fmdl.mat_idx{2}) = 2;

show_fem(img); view(-16,22); set(gca,'Projection','perspective')
print_convert surface_sim01a.png '-density 100'
