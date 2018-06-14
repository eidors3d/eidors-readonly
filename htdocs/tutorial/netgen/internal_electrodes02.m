extra={'ctr_el','solid ctr_el = sphere(0.2,0.2,0.5;0.08);'};
fmdl = ng_mk_cyl_models(1,[8,0.5],0.1,extra);


% Get electrode boundary
ctr_elec = fmdl.elems(fmdl.mat_idx{2},:);
bdy_elec = find_boundary(ctr_elec);
elec_nod = unique(bdy_elec(:));
fmdl.electrode(end+1) = struct( ...
     'nodes', elec_nod, 'z_contact',.01);
% Make sure EIDORS knows you have internal electrodes
fmdl.system_mat_fields.CEM_boundary = bdy_elec;

show_fem(fmdl); view(90,60);
print_convert internal_electrodes02a.jpg
