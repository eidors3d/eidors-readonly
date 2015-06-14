function fmdl = ng_mk_common_model(mdl_type, ...
               mdl_shape, elec_pos, elec_shape);
% NG_MK_COMMON_MODEL: utility to create common models
% fmdl = ng_mk_cyl_model(mdl_type,
%              mdl_shape, elec_pos, elec_shape);
%
% MDL_TYPE = 'CYL'
%    mdl_shape = {height, [radius, [maxsz]]}
%       if height = 0 -> calculate a 2D shape
%       radius (OPT)  -> (default = 1)
%       maxsz  (OPT)  -> max size of mesh elems (default = course mesh)
%   
%    ELECTRODE POSITIONS:
%     elec_pos = [n_elecs_per_plane,z_planes] 
%        OR
%     elec_pos = [degrees,z] centres of each electrode (N_elecs x 2)
%
% ELECTRODE SHAPES::
%    elec_shape = [width,height, maxsz]  % Rectangular elecs
%       OR
%    elec_shape = [radius, 0, maxsz ]    % Circular elecs
%       OR
%    elec_shape = [0, 0, maxsz ]         % Point elecs
%
% MDL_TYPE = 'CIRC'
%    Same as CYL interface, but uses netgen 2D interface. 
%    For this reason cyl height must be zero
%     elec_shape = [width, RFNUM]
%    where RFNUM is an optional refinement parameter (larger = more) 
%
%    For point electrodes, do not specify too small a maxh (ie < .1),
%      otherwise netgen has a problem and crashes
%
% USAGE EXAMPLES:
% Simple 3D cylinder. Radius = 1. No electrodes
%   fmdl= ng_mk_common_model('cyl',2,[8,.5],[0.05]); 

% (C) Andy Adler, 2015.  Licenced under GPL v2 or v3
% $Id$

if ischar(mdl_type) && strcmp(mdl_type,'UNIT_TEST'); do_unit_test; return; end

if isstruct(mdl_type)
   error('No code yet to process options structs');
else switch lower(mdl_type)
      case 'cyl'
         [body_geom, elec_geom, pp] = cyl_geom( mdl_shape, elec_pos, elec_shape);
         fmdl = ng_mk_geometric_models(body_geom, elec_geom);
         fmdl.body_geometry      = body_geom;
         fmdl.electrode_geometry = elec_geom;
      case 'circ'
         [shape, elec_pos, elec_shape] = circ_geom( mdl_shape, elec_pos, elec_shape);
         fmdl = ng_mk_2d_model(shape, elec_pos, elec_shape);
      otherwise
         error('mdl_type = (%s) not available (yet)', mdl_type);
end; end

   if mdl_dim(fmdl)==3 && pp.is2D
      fmdl = mdl2d_from3d(fmdl);
   end

% For 2D circle
function [shape, elec_pos, elec_shape] = circ_geom( mdl_shape, elec_pos, elec_shape);
   if mdl_shape(1) ~=0; error('specifying "circ" implies height of 0'); end
   if     length(mdl_shape)==1
      rad = 1;
      maxh = 2*pi*1 / 16;
   elseif length(mdl_shape)==2
      rad = mdl_shape(2);
      maxh = 2*pi*rad / 16;
   else
      rad = mdl_shape(2);
      maxh= mdl_shape(3);
   end


   if     size(elec_pos,1) == 0;
      theta= [];
   elseif size(elec_pos,1) == 1;
      theta = linspace( 0, 2*pi, elec_pos(1)+1)'; theta(end)=[];
   else
      theta = pi/180*elec_pos(:,1);
   end

   elec_pos= rad*[sin(theta), cos(theta)];

   % Point electrodes need refinement specified
   if length(elec_shape) == 1 && elec_shape(1) == 0
      elec_shape(2) = 10; % recommended default from 2d code
   end
    
  % If the boundary point is near the electrode, then unnecessary
  % refinement will occur. This could be avoided by placing the
  % boundary nodes carefully, but that introduces edge cases.
   
   theta_pts = ceil(2*pi*rad/maxh) - length(theta);
   if theta_pts > 0 % if we have extra points to distribute
      theta_  = [theta; 2*pi+theta(1)];
      per_interval = round( theta_pts * diff(theta_)/2/pi );
      for i=1:length(per_interval)
         new_pts = linspace( theta_(i), theta_(i+1), per_interval(i)+2);
         theta = [theta; new_pts(2:end-1)']; 
      end
      theta = sort(theta);
%  theta = linspace( 0, 2*pi, ceil(2*pi*rad/maxh))'; theta(end)=[];
   end

   % Must specify points starting like this: make then CCW
   shape = {rad*[sin(theta),-cos(theta)], maxh};


function [body_geom, elec_geom, pp] = cyl_geom( cyl_shape, elec_pos, elec_shape);
   [body_geom, pp] = parse_shape(cyl_shape);
   elec_geom = parse_elecs(elec_pos, elec_shape, pp);


function [body_geom,pp] = parse_shape(cyl_shape);
   tank_height = cyl_shape(1);
   tank_radius = 1;
   tank_maxh   = 0;
   is2D = 0;

   if length(cyl_shape)>1;
      tank_radius=cyl_shape(2);
   end
   if length(cyl_shape)>2; 
      tank_maxh  =cyl_shape(3);
   end
   if tank_height==0;
      is2D = 1;

      %Need some width to let netgen work, but not too much so
      % that it meshes the entire region
      tank_height = tank_radius/5; % initial extimate
      if tank_maxh>0
         tank_height = min(tank_height,2*tank_maxh);
      end
   end

   body_geom.cylinder.bottom_center = [0 0 0];
   body_geom.cylinder.top_center    = [0 0 tank_height];
   body_geom.cylinder.radius        = tank_radius;
   if tank_maxh > 0
      body_geom.max_edge_length     = tank_maxh;
   end

   pp.is2D = is2D; % put it here too.
   pp.tank_height= tank_height;
   pp.tank_radius= tank_radius;
   pp.body_geom  = body_geom;

% ELECTRODE POSITIONS:
%  elec_pos = [n_elecs_per_plane,z_planes] 
%     OR
%  elec_pos = [degrees,z] centres of each electrode (N_elecs x 2)
%
% ELECTRODE SHAPES::
%  elec_shape = [width,height, {maxsz}]  % Rectangular elecs
%     OR
%  elec_shape = [radius, {0, maxsz} ]  % Circular elecs
%     maxsz  (OPT)  -> max size of mesh elems (default = courase mesh)
% 
function elecs= parse_elecs(elec_pos, elec_shape, pp)
   is2D = pp.is2D;
   elec_geom = struct([]);
   hig = pp.tank_height;
   rad = pp.tank_radius;
   body_geom = pp.body_geom;
    
    n_elecs= size(elec_pos,1); 
    
    if n_elecs == 0
      elecs= {};
      return;
    end
   
   if is2D
      elec_pos(:,2) = hig/2;
   end

   % It never makes sense to specify only one elec
   % So elec_pos means the number of electrodes in this case
   if size(elec_pos,1) == 1
       % Parse elec_pos = [n_elecs_per_plane,z_planes] 
      n_elecs= elec_pos(1); % per plane
      th = linspace(0,2*pi, n_elecs+1)'; th(end)=[];

      on_elecs = ones(n_elecs, 1);
      el_th = []; 
      el_z  = []; 
      if length(elec_pos) == 1; % Height was forgotten, assume middle
         elec_pos(:,2) = hig/2;
      end
      for i=2:length(elec_pos)
        el_th = [el_th; th];
        el_z  = [el_z ; on_elecs*elec_pos(i)];
      end

   else
      el_th = elec_pos(:,1)*2*pi/360;
      el_z  = elec_pos(:,2);
   end
      
   n_elecs= size(el_z,1); 

   if size(elec_shape,1) == 1
      elec_shape = ones(n_elecs,1) * elec_shape;
   end

   elecs= {};
   for i= 1:n_elecs
     row = elec_shape(i,:); 
     elecs{i} = elec_spec( row, is2D, hig, rad, el_th(i), el_z(i) );
   end

% use clockwise coordinate system

function elec = elec_spec( row, is2D, hig, rad, el_th, el_z )
  xy_centre = rad*[sin(el_th),cos(el_th)];
  if     is2D
     if row(1) == 0;
        elec.point = [xy_centre, 0];
     else
        el_z = hig/2;
        row(2)=hig;
        elec.intersection.cylinder(1).top_center     = [0,0,2*row(1)];
        elec.intersection.cylinder(1).bottom_center  = [0,0,-2*row(1)];
        elec.intersection.cylinder(1).radius         = rad;
        elec.intersection.cylinder(1).complement_flag= 1;
        rad_dir = [sin(el_th);cos(el_th)];
        tan_dir = [0,1;-1,0]*rad_dir;
        elec.intersection.parallelepiped(1).vertex   = [0.97*xy_centre' - row(1)/2*tan_dir; el_z - row(2)/2];
        elec.intersection.parallelepiped(1).vector_a = [0; 0; 1];
        elec.intersection.parallelepiped(1).vector_b = [rad_dir; 0];
        elec.intersection.parallelepiped(1).vector_c = [tan_dir; 0];
        elec.intersection.parallelepiped(2).vertex   = [1.03*xy_centre' + row(1)/2*tan_dir; el_z + row(2)/2];
        elec.intersection.parallelepiped(2).vector_a =-[0; 0; 1];
        elec.intersection.parallelepiped(2).vector_b =-[rad_dir; 0];
        elec.intersection.parallelepiped(2).vector_c =-[tan_dir; 0];

        elec.enter_body_flag = 0;
     end
  else
     if row(1) == 0
        elec.point = [xy_centre, el_z];
     elseif length(row)<2 || row(2) == 0 % Circular electrodes 
        elec.cylinder.top_center     = [1.03*xy_centre, el_z];
        elec.cylinder.bottom_center  = [0.97*xy_centre, el_z];
        elec.cylinder.radius         = row(1);
        elec.enter_body_flag = 0;
     elseif row(2)>0      % Rectangular electrodes
% parallelepiped:     A parallelepiped is described by the following
%                     subfields: vertex ([0; 0; 0]), vector_a ([1; 0; 0]),
%                     vector_b ([0; 1; 0]), vector_c ([0; 0; 1]),
%                     complement_flag (false).
        rad_dir = [sin(el_th);cos(el_th)];
        tan_dir = [0,1;-1,0]*rad_dir;
        elec.intersection.parallelepiped(1).vertex   = [0.97*xy_centre' - row(1)/2*tan_dir; el_z - row(2)/2];
        elec.intersection.parallelepiped(1).vector_a = [0; 0; 1];
        elec.intersection.parallelepiped(1).vector_b = [rad_dir; 0];
        elec.intersection.parallelepiped(1).vector_c = [tan_dir; 0];
        elec.intersection.parallelepiped(2).vertex   = [1.03*xy_centre' + row(1)/2*tan_dir; el_z + row(2)/2];
        elec.intersection.parallelepiped(2).vector_a =-[0; 0; 1];
        elec.intersection.parallelepiped(2).vector_b =-[rad_dir; 0];
        elec.intersection.parallelepiped(2).vector_c =-[tan_dir; 0];
     else
        error('negative electrode width not supported');
     end
  end

  if length(row)>=3 && row(3) > 0
     elec.max_edge_length     = row(3);
  else
     % Do nothing
  end


function do_unit_test
  for tn = 1:do_test_number(0)
     fprintf('>>> ng_mk_cyl_models: TEST #%d\n',tn);
     fmdl= do_test_number(tn);
     subplot(3,3,rem(tn-1,9)+1)
     title(sprintf('test #%d',tn));
     show_fem(fmdl); drawnow
  end

function fmdl= do_test_number(tn)
   eidors_msg('ng_mk_cyl_models: UNIT_TEST #%d',tn,1);
   switch tn
   case 1;
% Simple 3D cylinder. Radius = 1. No electrodes
    fmdl= ng_mk_common_model('cyl',3,[0],[]); 

   case 2;
% Simple 2D cylinder. Radius = 2. Set minsize to refine
    fmdl= ng_mk_common_model('cyl',[0,2,.2],[0],[]); 

   case 3;
% 3D cylinder. Radius = 1. 2 planes of 8 elecs with radius 0.1
    fmdl= ng_mk_common_model('cyl',3,[8,1,2],[0.1]); 

   case 4;
% 3D cylinder. Radius = 1. 6 circ elecs with elec refinement
    fmdl= ng_mk_common_model('cyl',3,[7,1],[0.2,0,0.05]); 

   case 5;
% 3D cylinder. Radius = 1. 7 rect elecs with no refinement
    fmdl= ng_mk_common_model('cyl',3,[7,1],[0.2,0.3]); 

   case 6;
% 2D cylinder. Radius = 1. 11 rect elecs with refinement
    fmdl= ng_mk_common_model('cyl',0,[11],[0.2,0,0.05]); 

   case 7;
% 2D cylinder. Radius = 1.5. Refined(0.1). 11 elecs with refinement
    fmdl= ng_mk_common_model('cyl',[0,1,0.1],[11],[0.2,0,0.02]); 

   case 8;
% 2D cylinder. elecs at 0, 90 and 120 degrees
    fmdl= ng_mk_common_model('cyl',0,[0;90;120],[0.2,0,0.03]); 

   case 9;
% 2D cylinder. elecs at 0 (large,refined) and 120 (small) degrees
    fmdl= ng_mk_common_model('cyl',0,[0;120],[0.4,0,0.01;0.1,0,0.1]); 

   case 10;
% 3D cylinder. elecs at 0, 30, 60, 90 in planes
    fmdl= ng_mk_common_model('cyl',3,[0,0.5;30,1;60,1.5;90,2.0],[0.2,0,0.1]); 

   case 11;
% 3D cylinder. Various elecs at 0, 30, 60, 90 in planes
    el_pos = [0,0.5;30,1;60,1.5;90,2.0];
    el_sz  = [0.2,0,0.1;0.1,0,0.05;0.2,0.2,0.02;0.2,0.4,0.5];
    fmdl= ng_mk_common_model('cyl',3,el_pos,el_sz); 


   case 12;
    fmdl = ng_mk_common_model('circ', [0,1,.1], [4], [0.4,10]);
   case 13;
    fmdl = ng_mk_common_model('circ', 0, [9], [0.2,3]);
   case 14;
    fmdl = ng_mk_common_model('circ',[0,1,.10], [9], [0,1]);
   case 15;
    fmdl = ng_mk_common_model('circ',[0,1,.10], [9], [0]);
   case 16;
    fmdl = ng_mk_common_model('circ',[0,2,.05], [9], [0.1,100]);

   case 0; fmdl = 16; %%%% RETURN MAXIMUM
   otherwise;
     error('huh?')
   end
  
   if exist('img'); fmdl = img; end
