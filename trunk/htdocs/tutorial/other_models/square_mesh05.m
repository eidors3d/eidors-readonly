% $Id$

z_contact= 0.01;
nodes_per_elec= 5;
n_elec=17;
elec_width= 0.2;
elec_spacing= 1.0;

xllim=-12; xrlim= 12; ydepth=-15;
[x,y] = meshgrid( linspace(xllim,xrlim,49), linspace(ydepth,0,31) );
vtx= [x(:),y(:)];
% Refine points close to electrodes - don't worry if points overlap
[x,y] = meshgrid( -9:.25:9, -3:.25:0 );
vtx= [vtx; x(:),y(:)];
[x,y] = meshgrid( -9:.25:9, -15:.25:-11 );
vtx= [vtx; x(:),y(:)];

xgrid=  linspace(-elec_width/2, +elec_width/2, nodes_per_elec)';
x2grid= elec_width* [-5,-4,-3,3,4,5]'/4;
for i=1:n_elec
  x0= (i-1-(n_elec-1)/2)*elec_spacing;

% Top electrode
  y0  = zeros(size(xgrid));
  y0_2= zeros(size(x2grid));
% elec_nodes{2*i-1}= [x0+ xgrid, y0];
  elec_nodes{2*i-1}= [x0, 0];
  vtx= [ vtx; ...
        [x0 + x2grid    , y0_2              ];
        [x0 + xgrid*1.5 , y0  - elec_width/2];
        [x0 + x2grid*1.5, y0_2- elec_width/2];
        [x0 + xgrid*2   , y0  - elec_width  ];
        [x0 + xgrid*2   , y0  - elec_width*2]];

% Bottom electrode
  y0  = -13*ones(size(xgrid));
  y0_2= -13*ones(size(x2grid));
  elec_nodes{2*i}= [x0,-13]; % Only point electrodes insidq
  vtx= [ vtx; ...
        [x0 + x2grid    , y0_2              ];
        [x0 + xgrid*1.5 , y0  - elec_width/2];
        [x0 + x2grid*1.5, y0_2- elec_width/2];
        [x0 + xgrid*1.5 , y0  + elec_width/2];
        [x0 + x2grid*1.5, y0_2+ elec_width/2];];
end

fmdl= mk_fmdl_from_nodes( vtx, elec_nodes, z_contact, 'sq_m1');
fmdl.solve=@aa_fwd_solve;
fmdl.system_mat=@aa_calc_system_mat;
fmdl.jacobian=@aa_calc_jacobian;


subplot(121)
show_fem(fmdl); axis image
subplot(122)
show_fem(fmdl); axis image; axis([-2 2 -14 -12]);

print -dpng -r150 square_mesh05a.png
