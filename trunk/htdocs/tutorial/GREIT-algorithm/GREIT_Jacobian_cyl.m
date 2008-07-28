function [J,map,vbkgnd] = GREIT_Jacobian_cyl;
% Calculate the GREIT 32x32 Jacobian for a cylinder
% The FEM Model ng_mdl_16x1_fine must be available
%
% (C) 2008 Andy Adler. Licenced under GPL v2 or v3
% $Id$

if exist('GREIT_Jacobian_cyl.mat','file');
   load GREIT_Jacobian_cyl.mat J map vbkgnd
else
   [J,vbkgnd,map] = Jacobian_calc;
   save GREIT_Jacobian_cyl.mat J map vbkgnd
end


function [J,vbkgnd,map] = Jacobian_calc;
use_3d_model = 0;
if use_3d_model % This 3D model has some problems
   load ng_mdl_16x1_fine;  fmdl= ng_mdl_16x1_fine;
else
   imdl = mk_common_model('f2d3c',16); fmdl= imdl.fwd_model;
   fmdl.nodes = fmdl.nodes(:,[2,1]);
end

   % yvec is reversed because image yaxis is reversed
   fmdl.nodes(:,1) = -fmdl.nodes(:,1);

   pixel_grid= 32;
   nodes= fmdl.nodes;
   xyzmin= min(nodes,[],1);  xyzmax= max(nodes,[],1);
   xvec= linspace( xyzmin(1), xyzmax(1), pixel_grid+1);
   yvec= linspace( xyzmin(2), xyzmax(2), pixel_grid+1);

   % CALCULATE MODEL CORRESPONDENCES
if use_3d_model
   zvec= [0.6*xyzmin(3)+0.4*xyzmax(3), 0.4*xyzmin(3)+0.6*xyzmax(3)];
   [rmdl,c2f] = mk_grid_model(fmdl, xvec, yvec, zvec);
else
   [rmdl,c2f] = mk_grid_model(fmdl, xvec, yvec);
end

   img= eidors_obj('image','GREIT-ng_mdl');
   img.fwd_model= fmdl;
   img.fwd_model.coarse2fine = c2f;
   img.rec_model= rmdl;
   img.elem_data= ones(size(img.fwd_model,1));

   % ADJACENT STIMULATION PATTERNS
   img.fwd_model.stimulation= mk_stim_patterns(16, 1, ...
                [0,1],[0,1], {'do_redundant', 'no_meas_current'}, 1);

   % SOLVERS
   img.fwd_model.system_mat= @aa_calc_system_mat;
   img.fwd_model.solve=      @aa_fwd_solve;
   img.fwd_model.jacobian=   @aa_calc_jacobian;

   vbkgnd = fwd_solve(img);
   vbkgnd = vbkgnd.meas;
   J= calc_jacobian(img);

if use_3d_model
   map = reshape(sum(c2f,1),pixel_grid,pixel_grid)>0;
else % need to exclude some of the boundary
   [x,y]= meshgrid(linspace(-1,1,32),linspace(-1,1,32));
   map = x.^2 + y.^2 < 1.1;
end
