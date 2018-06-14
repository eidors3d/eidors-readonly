shape_str = [ ...
  'solid hole = sphere(0.2,0.2,1.5;0.08);' ... 
  'solid cyl    = cylinder (0,0,0; 0,0,1; 1) and not hole; \n', ...
  'solid mainobj= orthobrick(-2,-2,0;2,2,2) and cyl -maxh=0.3;\n'];
th = linspace(0,2*pi,7)'; th(end) = []; cs = [cos(th), sin(th)];
elec_pos = [  cs, th/2/pi + 0.5, cs, 0*th];
elec_shape=[0.1]; elec_obj = 'cyl';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

% Put an electrode in the hole
el_nodes = fmdl.nodes - ones(num_nodes(fmdl),1)*[0.2,0.2,1.5];
shim = 1e-4; % a little bit extra to catch elements
el_nodes = sum(el_nodes.^2,2) <= .08^2 + shim;
fmdl.electrode(end+1) = struct( ...
    'nodes', find(el_nodes), ...
    'z_contact',.01);

show_fem(fmdl); view(90,60);
print_convert internal_electrodes01a.jpg
