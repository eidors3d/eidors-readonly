% $Id$
xl=-3; xr= 3; yb=-15; yt= 15;
np= 35;
[x,y] = meshgrid( linspace(xr,xl,np), linspace(yb,yt,61) );
vtx= [y(:),x(:)];

for i=1:np
   elec_nodes{i   }= [y(1,i);x(1,i)]';
   elec_nodes{i+np}= [y(end,i);x(end,i)]';
end

z_contact= 1e5;
fmdl= mk_fmdl_from_nodes( vtx, elec_nodes, z_contact, 'sq_m1');
fmdl.stimulation(1).stimulation='Amp';
fmdl.stimulation(1).stim_pattern=[ones(np,1);-ones(np,1)];
fmdl.stimulation(1).meas_pattern=zeros(1,2*np); % don't care

% Add non-conductive target
ctr = interp_mesh( fmdl,0); xctr= ctr(:,1); yctr= ctr(:,2);
img = mk_image( fmdl,  ones(length(xctr),1) );

% Solve and add streamlines
img.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img);

img.fwd_model.mdl_slice_mapper.npx = 200;
img.fwd_model.mdl_slice_mapper.npy = 100;
q = show_current(img,vh.volt);
hh=show_fem(img);
set(hh,'EdgeColor',[1,1,1]*.75);
hold on;
sy = linspace(-3,3,15); sx =  15+0*sy;
hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
hold off;

% deformation
th = fmdl.nodes(:,1)/yt*(pi);
y = fmdl.nodes(:,2); y = y.*(th<=0) - y.*(th>0);
th = (th-pi/2).*(th<=0) + (pi/2-th).*(th>0);
[x,y] = pol2cart(th, y+4);
y = y+8.*(fmdl.nodes(:,1)<=0);
img2 = img; img2.fwd_model.nodes = [1.5*x, y];

figure;
show_fem(img2); axis image
img2.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img2);

img2.fwd_model.mdl_slice_mapper.npx = 200;
img2.fwd_model.mdl_slice_mapper.npy = 100;
q=  show_current(img2,vh.volt);
hh=show_fem(img2);
set(hh,'EdgeColor',[1,1,1]*.75);
hold on;
sx =  linspace(1.1,6.9,15); sy =  0*sx;
sy = -linspace(1.1,6.9,15); sx =  0*sy;
hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
hold off;


%Anisotropic material in deformed image
figure;
n_elems=size(img2.fwd_model.elems,1);
for i=1:n_elems
    elem_data(i,1,1,1) =1;
    elem_data(i,1,1,2) =0; elem_data(i,1,2,1) =0;
    elem_data(i,1,2,2) =1;       
end
%These elements are roughly in the 'middle of the bend'
for i=1700:2100
    elem_data(i,1,1,1) =1;
    elem_data(i,1,1,2) =0; elem_data(i,1,2,1) =0;
    elem_data(i,1,2,2) =10;       
end

fmdl=img2.fwd_model;
fmdl.system_mat = @system_mat_higher_order_anisotropy;
fmdl.approx_type    = 'tri3'; % linear
img3 = mk_image_anisotropy(fmdl,elem_data);
img3.fwd_solve.get_all_meas = 1; %Internal voltage
v3 = fwd_solve(img3); 
v3e=v3.meas; v3all=v3.volt;

img3.fwd_model.mdl_slice_mapper.npx = 200;
img3.fwd_model.mdl_slice_mapper.npy = 100;
q=  show_current(img3,v3all);
hh=show_fem(img2);
set(hh,'EdgeColor',[1,1,1]*.75);
hold on;
sx =  linspace(1.1,6.9,15); sy =  0*sx;
sy = -linspace(1.1,6.9,15); sx =  0*sy;
hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
hold off;















