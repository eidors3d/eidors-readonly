% Load data and positions
gps = load('Mine_20FEV2004.gps');
data= load('Mine_20FEV2004_LI.tomel');

% Forward Model
shape_str = ['solid top    = plane(0,0,0;0,1,0);\n' ...
             'solid mainobj= top and orthobrick(-100,-200,-100;425,10,100) -maxh=20.0;\n'];
elec_pos = gps(:,2:4); e0 = elec_pos(:,1)*0;
elec_pos = [  elec_pos, e0, e0+1, e0 ]; 
elec_shape=[0.5,.5,.5];
elec_obj = 'top';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

fmdl.stimulation = stim_meas_list( data(:,3:6) - 40100);

show_fem(fmdl);
