% Forward Model
shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
             'solid mainobj= top and orthobrick(-100,-100,-200;410,100,0) -maxh=50.0;\n'];
n_elec= 54; 
e0 = linspace(0,310,n_elec)';
elec_pos = [e0,0*e0,0*e0,0*e0,0*e0,1+0*e0];
elec_shape= [0.005,0,1];
elec_obj = 'top';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

% homogeneous conductivity of 0.1 Sm and insert a conductive sphere
img = mk_image(fmdl,0.1 + mk_c2f_circ_mapping(fmdl,[100;0;-50;20])*100);

hh= show_fem(img); view(-28,22); set(hh,'LineWidth',1);
print_convert stimPatt_PseudoSec01.png '-density 75'
