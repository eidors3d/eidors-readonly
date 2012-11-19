% Linear model $Id$
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


show_fem(img); axis image

print_convert anisotropy1_01a.png '-density 125'



% deformation
th = fmdl.nodes(:,1)/yt*(pi/2);
y = fmdl.nodes(:,2); y = y.*(th<=0) - y.*(th>0);
th = (th-pi/2).*(th<=0) + (pi/2-th).*(th>0);
[x,y] = pol2cart(th, y+4);
y = y+8.*(fmdl.nodes(:,1)<=0);
img2 = img; img2.fwd_model.nodes = [x, y*4];

show_fem(img2); axis image

print_convert anisotropy1_01b.png '-density 125'
