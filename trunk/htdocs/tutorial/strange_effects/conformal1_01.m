% Linear model $Id$
xl=-3; xr= 3; yb=-15; yt= 15;
[x,y] = meshgrid( linspace(xl,xr,7), linspace(yb,yt,31) );
vtx= [y(:),x(:)];

elec_nodes{1}= [y(1,:);x(1,:)]';
elec_nodes{2}= [y(end,:);x(end,:)]';

z_contact= 0.01;
fmdl= mk_fmdl_from_nodes( vtx, elec_nodes, z_contact, 'sq_m1');
fmdl.stimulation(1).stimulation='mA';
fmdl.stimulation(1).stim_pattern=[1;-1];
fmdl.stimulation(1).meas_pattern=[1,-1];

% Add non-conductive target
ctr = interp_mesh( fmdl,0); xctr= ctr(:,1); yctr= ctr(:,2);
img = mk_image( fmdl,  ones(length(xctr),1) );
img.elem_data( yctr>-3 & yctr<2 & xctr>8   & xctr<10 ) = 0.01;
img.elem_data( yctr> 0 & yctr<2 & xctr>-10 & xctr<10 ) = 0.01;


subplot(221)
show_fem(img); axis image

%print -dpng -r125 rpi_data01a.png
print -depsc2 jnk.eps;!LD_LIBRARY_PATH="" convert -density 125 jnk.eps conformal1_01a.png



% Confromal deformation
z= fmdl.nodes(:,1) + 1i*fmdl.nodes(:,2);
z= exp((z-(20+80i))/100).*(z+20i).*(z-10i);
img2 = img; img2.fwd_model.nodes = [real(z), imag(z)];

show_fem(img2); axis image

%print -dpng -r125 rpi_data01a.png
print -depsc2 jnk.eps;!LD_LIBRARY_PATH="" convert -density 125 jnk.eps conformal1_01b.png
