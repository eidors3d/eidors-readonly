clear;
shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1.0); \n', ...
             'solid tank   = orthobrick(-2,-2,0;2,2,0.4) and cyl; \n', ...
             'solid obj1   = cylinder( -.4,-.4,0; -.4,-.4,1;0.1);\n' ...
             'solid obj1t  = obj1 and tank; tlo obj1t;\n' ...
             'solid fish   = ellipsoid(0.2,0.2,0.2;0.2,0,0;0,0.05,0;0,0,0.05); tlo fish;\n', ...
             'solid mainobj= tank and not fish and not obj1 -maxh=0.3;\n'];
fx=0;fy=0;fl=20/75;fh=3/75;fw=2/75;
fx=0;fy=.7;fz=0.15;fl=20/75;fh=5/75;fw=5/75;
fish = sprintf('ellipsoid(%f,%f,%f;%f,0,0;0,%f,0;0,0,%f)',fx,fy,fz,fl/2,fh/2,fw/2);
shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1.0); \n', ...
             'solid tank   = orthobrick(-2,-2,0;2,2,0.3) and cyl; \n', ...
             'solid fish   = ',fish,'; tlo fish;\n', ...
             'solid mainobj= tank and not fish -maxh=0.3;\n'];

elec_shape([1,2])=[0.01];
elec_pos(1,:) = [ fx-fl/2,fy,fz,-1,0,0]; elec_obj{1} = 'fish';
elec_pos(2,:) = [ fx+fl/2,fy,fz, 1,0,0]; elec_obj{2} = 'fish';

n_elec = 8;
th = linspace(0,2*pi,n_elec+1)'; th(end) = [];
cs = [cos(th), sin(th)];
elec_pos = [elec_pos;  cs, 0.2+0*th, cs, 0*th];
for i=2+(1:n_elec);
   elec_obj{i} = 'cyl';
   elec_shape(i)=[0.01];
end
[fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);


fmdl.stimulation(1).stim_pattern = [1;-1;zeros(n_elec,1)];
fmdl.stimulation(1).meas_pattern = sparse([ ...
    0, 0, 1, 0, 0, 0,-1, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0,-1, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0,-1, 0;
    0, 0, 0, 0, 0, 1, 0, 0, 0,-1]);

img = mk_image( fmdl, 1);
%img.elem_data( mat_idx{1} ) = 8.1;
%img.elem_data( mat_idx{3} ) = 0.001;
show_fem(img);

img.fwd_solve.get_all_meas = 1;
vv = fwd_solve(img);

img_v = rmfield(img, 'elem_data');
img_v.node_data = vv.volt(:,1);
img_v.calc_colours.npoints = 128;

subplot(121)
show_slices(img_v,[inf,inf,fz/2]);
subplot(122)
show_slices(img_v,[fl/2,inf,inf]);


