fx=0;fy=.7;fz=0.15;fl=20/75;fh=5/75;fw=5/75;
fish = sprintf('ellipsoid(%f,%f,%f;%f,0,0;0,%f,0;0,0,%f)',fx,fy,fz,fl/2,fh/2,fw/2);
shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1.0); \n', ...
             'solid tank   = orthobrick(-2,-2,0;2,2,0.3) and cyl; \n', ...
             'solid obj1   = cylinder( -.3,.5,0; -.3,.5,1;0.1);\n' ...
             'solid obj1t  = obj1 and tank; tlo obj1t;\n' ...
             'solid fish   = ',fish,'; tlo fish;\n', ...
             'solid mainobj= tank and not fish and not obj1 -maxh=0.3;\n'];

% Electrodes on the head and tail
elec_shape([1,2])=[0.01];
elec_pos(1,:) = [ fx-fl/2,fy,fz,-1,0,0]; elec_obj{1} = 'fish';
elec_pos(2,:) = [ fx+fl/2,fy,fz, 1,0,0]; elec_obj{2} = 'fish';

% Electrodes on the tank
n_elec = 8;
th = linspace(0,2*pi,n_elec+1)'; th(end) = [];
cs = [cos(th), sin(th)];
elec_pos = [elec_pos;  cs, 0.2+0*th, cs, 0*th];
for i=2+(1:n_elec);
   elec_obj{i} = 'cyl';
   elec_shape(i)=[0.05];
end
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

show_fem(fmdl);view(90,20);
print_convert fish_tank01a.png '-density 60'
show_fem(fmdl);view(90,80);
print_convert fish_tank01b.png '-density 60'
axis([-.3,.3,.2,1,0,0.2])
print_convert fish_tank01c.png '-density 60'
