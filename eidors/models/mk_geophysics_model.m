function imdl = mk_geophysics_model(str, ne, opt);
% imdl = mk_geophysics_model(str, ne, [option])
%
% ne  - number of electrodes, 5 metre spacing (+5,+10,...)
%       and 1 metre diameter
%         OR
%       a list of electrode locations in the x-dimension
% str - model, x = see hmax_rec
%       h2x -   2D half-space, linear array (2d fwd, 2d rec mdl)
%       h2p5x - 2.5D half-space, linear array (3d fwd, 2d rec mdl)
%       h3x   - 3D half-space, linear array (3d fwd, 3d rec mdl) [TODO]
% opt - override default configuration options (optional cell array)
%       'hmax_rec' - coarse reconstruction mesh density, given array width xw
%                ['a' : hmax_rec=xw/5;
%                 'b' : hmax_rec=xw/10;
%                 'c' : hmax_rec=xw/15;
%                 'd' : hmax_rec=xw/20;
%                 'e' : hmax_rec=xw/25;
%                 'f' : hmax_rec=xw/30;
%                 'g' : hmax_rec=xw/35;
%                 'h' : hmax_rec=xw/40]
%       'hmax_fwd' - forward model mesh density [hmax_rec/2]
%       'elec_width' - electrode width [0.1 m]
%       'z_contact' - electrode contact impedance [0.01 \Ohm.m]
%       'elec_spacing' - distance between electrode centers [5 m]
%       'extend_x' - extra mesh in the principle axis of the
%                    electrode array, multiple of array width [1]
%       'extend_y' - extra mesh in the minor axis of the
%                    electrode array, multiple of array width
%                    (3D models only) [1]
%       'extend_z' - extra depth of model, multiple of array width [1]
%       'skip_c2f' - skip building the rec_model to fwd_model mapping [0]
%
% The linear electrode array runs in the +X direction at Z=0. For
% the 3D model, the Y-axis is perpendicular to the electrode array.
%
% (C) 2015, 2016 A. Boyle
% License: GPL version 2 or version 3

% model: 64 electrode, 2d half-space
% Once upon a time, this code started out from the following tutorial.
% model from http://eidors3d.sourceforge.net/tutorial/other_models/square_mesh.shtml

if ischar(str) && strcmp(str,'UNIT_TEST'); do_unit_test; return; end
copt.fstr = 'mk_geophysics_model';
if nargin < 3
   opt = {};
end
imdl = eidors_cache(@mk_model,{str, ne, opt}, copt);
imdl.fwd_model.stimulation = stim_pattern_geophys(ne, 'Wenner');

function imdl=mk_model(str,ne,opt);
if str(1) ~= 'h'
   error([str ': only know how to build linear half-space model: h***']);
end
if all(length(str) ~= [5 3])
   error([str ': wrong argument length']);
end

if (length(str) == 5) && (strcmp(str(2:4), '2p5'))
   FMDL_DIM=3;
   CMDL_DIM=2;
elseif str(2) == '2'
   FMDL_DIM=2;
   CMDL_DIM=2;
elseif str(2) == '3'
   FMDL_DIM=3;
   CMDL_DIM=3;
   error([str ': model TODO (3d-3d)']);
else
   error([str ': unrecognized dimensonality: hx*, for x=2,3, or 2p5']);
end
if str(end) > 'h' || str(end) < 'a'
   error([str ': unrecognized model name (hmax) ' str(1:end-1) 'x for x=a,b,c,d,e,f,g or h']);
end

skip_c2f = 0;
elec_width= 0.1;
z_contact= 0.01;
nodes_per_elec= 3; %floor(elec_width/hmax_rec*10);
elec_spacing= 5.0;

extend_x = 1;
extend_y = 1;
extend_z = 1;
if length(opt) > 0 % allow overriding the default values
   assert(round(length(opt)/2)*2 == length(opt),'option missing value?');
   expect = {'hmax_rec','hmax_fwd', 'elec_width','z_contact','elec_spacing',...
             'extend_x', 'extend_y', 'extend_z', 'skip_c2f'};
   opts = struct(opt{:})
   for i = fieldnames(opts)'
      assert(any(strcmp(i,expect)), ['unexpected option: ',i{:}]);
      eval([i{:} ' = opts.(i{:});']);
   end
end
if length(ne) == 1 % ne: number of electrodes
   xw=(ne-1)*elec_spacing; % array width
   %xs=-(ne-1)*elec_spacing/2; % array centered
   xs=+5; % array at left-most at +5
   x0l = xs+([1:ne]-1)*elec_spacing;
   save_model_to_disk=1;
else % ne: list of electrode x-coordinates
   xw=max(ne)-min(ne);
   xs=min(ne);
   x0l = ne;
   ne=length(ne);
   save_model_to_disk=0;
end
if length(opt) > 0
   save_model_to_disk=0;
end
if ~exist('hmax_rec','var')
   switch(str(end))
      case 'a'
         hmax_rec=xw/5;
      case 'b'
         hmax_rec=xw/10;
      case 'c'
         hmax_rec=xw/15;
      case 'd'
         hmax_rec=xw/20;
      case 'e'
         hmax_rec=xw/25;
      case 'f'
         hmax_rec=xw/30;
      case 'g'
         hmax_rec=xw/35;
      case 'h'
         hmax_rec=xw/40;
      otherwise
   end
end
if ~exist('hmax_fwd','var') % allow hmax_fwd to depend on configured hmax_rec
   hmax_fwd=hmax_rec/2.0; % avoid parametrization aliasing
end

if save_model_to_disk
   filename=sprintf('imdl-%s-%03del.mat',str,ne);
   if exist(filename, 'file') == 2
      tmp=load(filename);
      imdl = tmp.imdl;
      eidors_msg(sprintf('%s: %s, %d electrode model loaded from file',filename,str,ne));
      return
   end
end

assert(extend_x>0,'extend_x must be > 0');
assert(extend_y>0,'extend_y must be > 0');
assert(extend_z>-1,'extend_z must be > -1');
assert(hmax_rec*3 < xw, sprintf('(hmax_rec=%g * 3) must be < array width=%g',hmax_rec,xw));

% 2d cmdl
xllim=xs-extend_x*elec_spacing*ne;
xrlim=xs+xw+extend_x*elec_spacing*ne;
zdepth=-(xw+extend_z*elec_spacing*ne);
xr=floor((xrlim-xllim)/hmax_rec/2)*2+1; % odd number
yr=floor(-zdepth/hmax_rec/2)*2+1; % odd number
[x,y] = meshgrid( linspace(xllim,xrlim,xr), linspace(zdepth,0,yr));
vtx= [x(:),y(:)];
cmdl= mk_fmdl_from_nodes( vtx,{vtx(1,:)}, z_contact, 'sq_m2');

% fmdl refinement
xllim1=xllim+extend_x*3/5*elec_spacing*ne;
xrlim1=xrlim-extend_x*3/5*elec_spacing*ne;
zdepth1=-(xw+extend_z*3/5*elec_spacing*ne);
assert(zdepth1 > zdepth, 'zdepth: oops, inner mesh must be smaller than outer mesh');
if 0 && FMDL_DIM == 2  % 2D fmdl
   % old code using mk_fmdl_from_nodes
   % ... shit mesh around the electrodes, bad elements scattered throughout
   xr=floor((xrlim-xllim)/hmax_fwd/2)*2+1; % odd number
   yr=floor(-zdepth/hmax_fwd/2)*2+1; % odd number
   [x,y] = meshgrid( linspace(xllim,xrlim,xr), linspace(zdepth,0,yr));
   vtx= [x(:),y(:)];

   % Refine points close to electrodes - don't worry if points overlap
   xr=floor((xrlim1-xllim1)/hmax_fwd)*2+1; % odd number
   yr=floor(-zdepth1/hmax_fwd)*2+1; % odd number
   [x,y] = meshgrid( linspace(xllim1,xrlim1,xr), linspace(zdepth1,0,yr));
   vtx= [vtx; x(:),y(:)];

   % refine around the electrodes
   xgrid =  linspace(0, +3*elec_width, 3*nodes_per_elec-2)' - elec_width*3/2;
   ygrid =  linspace(0, -2*elec_width, 2*nodes_per_elec-1)';
   y0g= repelem(ygrid,length(xgrid));
   for i=1:ne
   % Electrode centre
     x0 = x0l(i); % constructed earlier
     x0g= repmat(x0 + xgrid, length(ygrid),1);
     vtx= [ vtx; [x0g, y0g]];
     ve=find((vtx(:,2)==0) & (vtx(:,1) >= x0-elec_width/2) & (vtx(:,1) <= x0+elec_width/2));
     elec_nodes{i}= vtx(ve,:);
   end

   fmdl= mk_fmdl_from_nodes( vtx, elec_nodes, z_contact, 'sq_m1');
   [~, gn] = min(fmdl.nodes(:,end));
   fmdl.gnd_node = gn; % make sure the ground node is away from surface electrodes

else % 3D fmdl
% NOTE Herve's interface won't let me build a box inside a box with different hmax_rec.. @#$%
%   body_geometry(1).ortho_brick = struct;
%   body_geometry(1).intersection.ortho_brick(1).opposite_corner_a = [ xllim  -10*elec_spacing      0 ];
%   body_geometry(1).intersection.ortho_brick(1).opposite_corner_b = [ xrlim  +10*elec_spacing zdepth ];
%%   body_geometry(1).intersection.ortho_brick(2).opposite_corner_a = [ xllim1+1 -5*elec_spacing     +1 ];
%%   body_geometry(1).intersection.ortho_brick(2).opposite_corner_b = [ xrlim1-1 +5*elec_spacing zdepth1];
%%   body_geometry(1).intersection.ortho_brick(2).complement_flag = 1;
%   body_geometry(1).max_edge_length = hmax_fwd;
%   body_geometry(2).ortho_brick.opposite_corner_a = [ xllim1 -1*elec_spacing      0 ];
%   body_geometry(2).ortho_brick.opposite_corner_b = [ xrlim1 +1*elec_spacing zdepth1];
%   body_geometry(2).max_edge_length = hmax_fwd/2.0;
%   for i=1:ne
%      x0 = xs+(ne-i)*elec_spacing;
%      electrode_geometry{i}.cylinder.radius        = elec_width/2;
%      electrode_geometry{i}.cylinder.top_center    = [x0, 0, +0.1];
%      electrode_geometry{i}.cylinder.bottom_center = [x0, 0, -0.1];
%      electrode_geometry{i}.max_edge_length = elec_width/(nodes_per_elec-1);
%   end
%   fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);

   xyz = [];
   for i=1:ne
      x0 = x0l(i); % calculated earlier
      xyz = [xyz; x0 0 0];
   end
   [R, X] = rot_line_to_xaxis(xyz);
   % rescale, centre electrodes so NetGen can be happy
   xyzc = (xyz - X)*R; % centre and scale electrodes: -1 to +1 y-axis

   % find the electrode locations
   % find the min/max and set limits for the domain that are
   %   20% bigger in every dimension
   xyz_box = [ min(xyzc); max(xyzc) ]; % a bounding box around the electrodes
   elec_mean_dist = norm(mean(diff(xyzc))); % average distance between electrodes
%   fprintf('  mean dist btw electrodes = %0.2f m\n', elec_mean_dist/norm(R));

   if 0
      rmi = elec_mean_dist*ne/4; % elec_mean_dist*4 - inner margin
      rmo = rmi.*[2 2 4]; % outer margin
      ri = xyz_box(1,:) - rmi; ri(2,:) = xyz_box(2,:) + rmi; ri(1,3)=-max(range(xyzc))/2;
      ro = xyz_box(1,:) - rmo; ro(2,:) = xyz_box(2,:) + rmo;
   else
      ro = ([ xllim   -extend_y*elec_spacing*ne zdepth;
              xrlim   +extend_y*elec_spacing*ne      1 ] - X(1:2,:))*R;
      ri = ([ xllim1  -extend_y*2/5*elec_spacing*ne zdepth1;
              xrlim1  +extend_y*2/5*elec_spacing*ne      2 ] - X(1:2,:))*R;
   end
   ps = [0 0 0; 0 0 1]; % surface plane
   ro_maxh   = hmax_fwd*norm(R);
   ri_maxh   = hmax_fwd*norm(R)/2.0;

   % build shape string for NetGen
   elec_pos   = [ xyzc(:,1:2), repmat([0, 0, 0, 1],ne,1) ]; % p(x,y,z=0), n(0,0,1)
   elec_shape = [elec_width*norm(R), 0, elec_width*norm(R)/(nodes_per_elec-1)]; % CEM, circular, maxh
   elec_obj     = 'ps';
   tlo = 'tlo ro'; % skip trailing ';\n'
   if FMDL_DIM == 3
      shape_str = [...
                   sprintf('solid ps = plane(%s);\n', a2s(ps)), ...
                   sprintf('solid bi = orthobrick(%s);\n', a2s(ri)), ...
                   sprintf('solid bo = orthobrick(%s);\n', a2s(ro)), ...
                   sprintf('solid ri = bi and ps -maxh=%f;\n', ri_maxh), ...
                   sprintf('solid ro = bo and ps and (not bi) -maxh=%f;\n', ro_maxh), ...
                   sprintf('solid mainobj = ri;\n')];
      % Note that ri must be the 'mainobj' so that it can intersect with the electrodes
      % additional top level objects for netgen
   else % netgen 2d model
      % to create the 2D slice we need to give NetGen something to work with
      %Need some width to let netgen work, but not too much so
      % that it meshes the entire region
      sw = range(ri(:,1)) / 5; % initial extimate
      sw = min(sw,2*ro_maxh); % coarse model slice width
      ri2d = ri; ri2d(:,2) = [-sw 0 ]; % [-sw/2; +sw/2 ];
      ro2d = ro; ro2d(:,2) = [-sw 0 ]; % [-sw/2; +sw/2 ];
      shape_str = [...
                   sprintf('solid ps = plane(%s);\n', a2s(ps)), ...
                   sprintf('solid bi = orthobrick(%s);\n', a2s(ri2d)), ...
                   sprintf('solid bo = orthobrick(%s);\n', a2s(ro2d)), ...
                   sprintf('solid ri = bi and ps -maxh=%f;\n', ri_maxh), ...
                   sprintf('solid ro = bo and ps and (not bi) -maxh=%f;\n', ro_maxh), ...
                   sprintf('solid mainobj = ri;\n')];
   end
   % fprintf('SHAPE_STR: %s', shape_str);
   [fmdl, mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj, tlo);
   if FMDL_DIM == 2 % 2D
      % now convert the roughly 2D slice into a true 2D plane
      [fmdl, mat_idx] = copy_mdl2d_from3d(fmdl, mat_idx, 'y');

      % reverse the centre and scaling
      nn = size(fmdl.nodes,1);
      Xn = repmat(X(1,[1 3]), nn, 1);
   else % 3D
      % c2f
      cmdl.mk_coarse_fine_mapping.f2c_offset  = [0 0 0];
      cmdl.mk_coarse_fine_mapping.f2c_project = [1 0 0; 0 0 1; 0 1 0];
      % duplicate parameters since mk_analytic/approx_c2f have different names...
      cmdl.mk_analytic_c2f.f2c_offset  = cmdl.mk_coarse_fine_mapping.f2c_offset;
      cmdl.mk_analytic_c2f.f2c_project = cmdl.mk_coarse_fine_mapping.f2c_project;

      % reverse the centre and scaling
      nn = size(fmdl.nodes,1);
      Xn = repmat(X(1,:), nn, 1);
   end
   fmdl.nodes = (fmdl.nodes / R) + Xn;

%   show_fem(fmdl); xlabel('x'); ylabel('y'); zlabel('z');

   [~, gn] = min(fmdl.nodes(:,end));
   fmdl.gnd_node = gn; % make sure the ground node is away from surface electrodes
   [~, gn] = min(cmdl.nodes(:,end));
   cmdl.gnd_node = gn; % make sure the ground node is away from surface electrodes
end

% stick electrode nodes into cmdl so that show_fem will plot them
for i=1:ne
   n=fmdl.electrode(i).nodes;
   nn=length(n);
   nx=fmdl.nodes(n,:);

   nnc = length(cmdl.nodes);
   cmdl.nodes = [cmdl.nodes; nx(:,[1 FMDL_DIM])];
   cmdl.electrode(i).nodes = (nnc+1):(nnc+nn);
   cmdl.electrode(i).z_contact = z_contact;
   fmdl.electrode(i).z_contact = z_contact;
end

% Note that the 2d fwd model mesh is a bit trashy around the electrode
% refinement...  there are bad quality elements just outside the per electrode
% refined areas

imdl= mk_common_model('a2d0c',ne); % 2d model
imdl.fwd_model = fmdl;
imdl.rec_model = cmdl;
% EIDORS "analytic_c2f" gets stuck, do an approximate one
%eidors_default('set','mk_coarse_fine_mapping','mk_analytic_c2f');
%eidors_default('set','mk_coarse_fine_mapping','mk_approx_c2f');
%[c2f,out] = mk_coarse_fine_mapping(fmdl,cmdl);
if ~skip_c2f
   [c2f,out] = mk_approx_c2f(fmdl,cmdl);
   imdl.fwd_model.coarse2fine = c2f;
   imdl.fwd_model.background = out;
end
imdl.name = ['EIDORS mk_geophysics_model ' str];
imdl.fwd_model.name = ['EIDORS mk_geophysics_model fwd model ' str];
imdl.rec_model.name = ['EIDORS mk_geophysics_model rec model ' str];

if save_model_to_disk
   save(filename,'imdl');
end

% convert 2x3 array to "x1,y1,z1;x2,y2,z2" string
% convert 1x3 array to "x1,y1,z1" string
% if 1x1 array, then use xy_ctr=x1,y1 and [0 0 1]=x2,y2,z2 (z+)
function s = a2s(a)
if length(a(:)) == 3
   s = sprintf('%f,%f,%f', ...
               a(1), a(2), a(3));
else
   s = sprintf('%f,%f,%f;%f,%f,%f', ...
                a(1,1), a(1,2), a(1,3), ...
               a(2,1), a(2,2), a(2,3));
end

% returns R rotation/scaling and X0 offset
% xyz1 = R * xyz + X; % rotate and scale to +/- 1
function [R,X] = rot_line_to_xaxis(xyz)
x = xyz(:,1); y=xyz(:,2); z=xyz(:,3);

% fit line to points
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/32502
p = mean(xyz);
[U,S,V] = svd([x-p(1), y-p(2), z-p(3)]);
if V(end,1) ~= 0
   N=1/V(end,1)*V(:,1);
else
   N=V(:,1);
end
A=p' + dot( xyz(1,  :) - p, N ) * N/norm(N)^2;
B=p' + dot( xyz(end,:) - p, N ) * N/norm(N)^2;

% rotate line to +y-axis
a = N/norm(N);
b = [ 1 0 0 ]; % +x
% http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
v = cross(a, b);
s = norm(v); c = dot(a, b); % sin, cos
xt = @(v) [   0  -v(3)  v(2); ... % skew symmetric cross-product of v
            v(3)    0  -v(1); ... % diagnoal is the scaling identity matrix
           -v(2)  v(1)    0];
if abs(s) < eps*1e3
   R = 1;
else
   R = (eye(3) + xt(v) + xt(v)^2 * (1-c)/s^2);
   R=R'; % for right multiply: xyz*R
end

X = repmat(p, size(xyz,1), 1);
R = R / max(range((xyz-X)*R)) * 2;

DEBUG=0;
if DEBUG
   clf;
   subplot(211);
   plot3(x,y,z,'bo');
   hold on;
   plot3(x(1),y(1),z(1),'go');
   plot3(p(1),p(2),p(3),'ro');
   plot3([A(1),B(1)],[A(2),B(2)],[A(3),B(3)]);
   grid; axis equal; hold off;
   xlabel('x'); ylabel('y'); zlabel('z');
   subplot(212);
   xx = (xyz-X)*R;
   plot3(xx(:,1), xx(:,2), xx(:,3),'bo');
   hold on;
   plot3(xx(1,1), xx(1,2), xx(1,3),'go');
   plot3(0,0,0,'ro');
   grid; axis equal; hold off;
   xlabel('x'); ylabel('y'); zlabel('z');
end

% the (max-min) range of a variable's values
function r = range(a)
r = max(a(:))-min(a(:));

function [mdl2,idx2] = copy_mdl2d_from3d(mdl3,idx3,xyz);
% AB: taken from EIDORS function ng_mk_gen_models() subfunction of the same name
% AB: NEW: xyz = 'x', 'y' or 'z' -- default was Z, we want X
   if xyz == 'x'
     % swap Z and X
     T = [ 0 0 1; 0 1 0; 1 0 0 ];
   elseif xyz == 'y'
     % swap Z and Y
     T = [ 1 0 0; 0 0 1; 0 1 0 ];
   elseif xyz == 'z'
     T = eye(3);
   else
     error('xyz must be "x", "y" or "z"');
   end
   mdl3.nodes = mdl3.nodes * T; % AB: SWAP axes

   % set name
   mdl2 = eidors_obj('fwd_model',sprintf('%s 2D',mdl3.name));

   % set nodes
   [bdy,idx] = find_boundary(mdl3.elems);
   vtx = mdl3.nodes;
   z_vtx = reshape(vtx(bdy,3), size(bdy) );
   z_vtx_thres = max(z_vtx(:))-10*eps*range(z_vtx(:));
   lay0  = find( all(z_vtx >= z_vtx_thres, 2) );
   bdy0  = bdy( lay0, :);

   vtx0  = unique(bdy0(:));
   mdl2.nodes = vtx(vtx0,1:2);

   % set elems
   nmap  = zeros(size(vtx,1),1); nmap(vtx0) = 1:length(vtx0);
   bdy0  = reshape(nmap(bdy0), size(bdy0) ); % renumber to new scheme
   mdl2.elems = bdy0;

   % set boundary
   mdl2.boundary = find_boundary( mdl2.elems);

   % set gnd_node
   mdl2.gnd_node = nmap(mdl3.gnd_node);

   % set material indices
   % TODO: vectorize code
   idx2 = {};
   idx0  = idx( lay0, :);
   for i=1:size(idx3,2)
     idx2{i} = [];
     ii = 1;
     for j=1:size(idx3{i},1)
         idx_tmp = find( idx0==idx3{i}(j) );
         if not(isempty(idx_tmp))
           idx2{i}(ii,1) = idx_tmp(1,1);
           ii = ii + 1;
         end
     end
   end

   % set electrode
   if isfield(mdl3,'electrode')
     mdl2.electrode = mdl3.electrode;
     for i=1:length(mdl2.electrode);
        enodes = nmap( mdl2.electrode(i).nodes );
        enodes(enodes==0) = []; % Remove 3D layers
        mdl2.electrode(i).nodes = enodes(:)';
     end
   end

   ignore = {'electrode', 'nodes', 'boundary', 'elems', 'gnd_node', 'boundary_numbers', 'mat_idx', 'mat_idx_reordered'};
   for n=fieldnames(mdl3)'
      if ~any(strcmp(n,ignore))
         mdl2.(n{:}) = mdl3.(n{:});
      end
   end

function do_unit_test
   ne = 16;
   imdl = mk_geophysics_model('h2p5a', ne);
   imdl.fwd_model.stimulation = stim_pattern_geophys(ne, 'Wenner');
   img = mk_image(imdl.fwd_model, 1);
   img.fwd_solve.get_all_meas = 1;
   vh = fwd_solve_halfspace(img);
   vd = fwd_solve(img);
   unit_test_cmp('h2a halfspace vs default TEST', norm(vh.meas - vd.meas), 0, 4e-3);
clf; h=plot([vh.meas vd.meas],'o--'); legend('analytic','FEM'); set(gca,'box','off'); set(h,'LineWidth',2);
