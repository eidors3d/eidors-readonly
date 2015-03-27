function [geom, fmdl] = ng_mk_cyl_model(mdl_type, ...
               mdl_shape, elec_pos, elec_shape);
% NG_MK_MODELS: utility to create common models
% [geom, fmdl] = ng_mk_cyl_model(mdl_type,
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
%
% USAGE EXAMPLES:
% Simple 3D cylinder. Radius = 1. No electrodes
%   fmdl= ng_mk_model('cyl',2,[8,.5],[0.05]); 

% (C) Andy Adler, 2015.  Licenced under GPL v2 or v3
% $Id$

if isstr(mdl_shape) && strcmp(mdl_shape,'UNIT_TEST'); do_unit_test; return; end

if isstruct(mdl_type)
   error('No code yet to process options structs');
else switch mdl_type
      case 'cyl'
         [body_geom, elec_geom] = cyl_geom( mdl_shape, elec_pos, elec_shape);
      otherwise
         error('mdl_type = (%s) not available', mdl_type);
end; end

fmdl = ng_mk_geometric_models(body_geom, elec_geom)
if geom.is2D
   fmdl = mdl2d_from3d(fmdl);
end

function fmdl = cyl_geom( cyl_shape, elec_pos, elec_shape);


   [tank_height, tank_radius, tank_maxh, is2D] = parse_shape(cyl_shape);
   [elecs, centres] = parse_elecs( elec_pos, elec_shape,  ...
                          tank_height, tank_radius, is2D );


function [tank_height, tank_radius, tank_maxh, is2D] = ...
              parse_shape(cyl_shape);
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
% OUTPUT:
%  elecs(i).pos   = [x,y,z]
%  elecs(i).shape = 'C' or 'R'
%  elecs(i).dims  = [radius] or [width,height]
%  elecs(i).maxh  = '-maxh=#' or '';
function [elecs, centres] = parse_elecs(elec_pos, elec_shape, hig, rad, is2D );
    
    n_elecs= size(elec_pos,1); 
    
    if n_elecs == 0
      elecs= struct([]); % empty
      centres= [];
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

   for i= 1:n_elecs
     row = elec_shape(i,:); 
     elecs(i) = elec_spec( row, is2D, hig, rad );
   end

%   %MC FIX 05.07.13 - COORDINATE SYSTEM - THETA = 0 AT 3PM AND THETA INCREASE ANTICLOCK
% NO - we will use clockwise coordinate system
   centres = [rad*sin(el_th),rad*cos(el_th),el_z];
   
   for i= 1:n_elecs; elecs(i).pos  = centres(i,:); end

   if n_elecs == 0
      elecs= struct([]); % empty
   end

function elec = elec_spec( row, is2D, hig, rad )
  if     is2D
     if row(1) == 0;
        elec.shape = 'P';
% To create a PEM, we make a square and take the corner. This isn't perfect, since
% the elec isn't quite where we asked for it, but that's as good is I can do. I tried
% asking for two rectangles to touch, but that freaks netgen out.
        elec.dims  =  [rad/20, hig]; 
     else
        elec.shape = 'R';
        elec.dims  = [row(1),hig];
     end
  else
     if row(1) == 0
        elec.shape = 'P' 
        elec.dims  = [rad/20, hig/10];
     elseif length(row)<2 || row(2) == 0 % Circular electrodes 
        elec.shape = 'C';
        elec.dims  = row(1);
     elseif row(2)>0      % Rectangular electrodes
        elec.shape = 'R';
        elec.dims  = row(1:2);
     else
        error('negative electrode width');
     end
  end

  if length(row)>=3 && row(3) > 0
     elec.maxh = sprintf('-maxh=%f', row(3));
  else
     elec.maxh = '';
  end


function do_unit_test
  for tn = 1:do_test_number(0)
     fprintf('>>> ng_mk_cyl_models: TEST #%d\n',tn);
     fmdl= do_test_number(tn);
     show_fem(fmdl);
  end

function fmdl= do_test_number(tn)
   eidors_msg('ng_mk_cyl_models: UNIT_TEST #%d',tn,1);
   switch tn
   case 1;
% Simple 3D cylinder. Radius = 1. No electrodes
    fmdl= ng_mk_models('cyl',3,[0],[]); 

   case 2;
% Simple 2D cylinder. Radius = 2. Set minsize to refine
    fmdl= ng_mk_models('cyl',[0,2,.2],[0],[]); 

   case 3;
% 3D cylinder. Radius = 1. 2 planes of 8 elecs with radius 0.1
    fmdl= ng_mk_models('cyl',3,[8,1,2],[0.1]); 

   case 4;
% 3D cylinder. Radius = 1. 6 circ elecs with elec refinement
    fmdl= ng_mk_models('cyl',3,[7,1],[0.2,0,0.05]); 

   case 5;
% 3D cylinder. Radius = 1. 7 rect elecs with no refinement
    fmdl= ng_mk_model('cyl',3,[7,1],[0.2,0.3]); 

   case 6;
% 2D cylinder. Radius = 1. 11 rect elecs with refinement
    fmdl= ng_mk_model('cyl',0,[11],[0.2,0,0.05]); 

   case 7;
% 2D cylinder. Radius = 1.5. Refined(0.1). 11 elecs with refinement
    fmdl= ng_mk_model('cyl',[0,1,0.1],[11],[0.2,0,0.02]); 

   case 8;
% 2D cylinder. elecs at 0, 90 and 120 degrees
    fmdl= ng_mk_model('cyl',0,[0;90;120],[0.2,0,0.03]); 

   case 9;
% 2D cylinder. elecs at 0 (large,refined) and 120 (small) degrees
    fmdl= ng_mk_model('cyl',0,[0;120],[0.4,0,0.01;0.1,0,0.1]); 

   case 10;
% 3D cylinder. elecs at 0, 30, 60, 90 in planes
    fmdl= ng_mk_model('cyl',3,[0,0.5;30,1;60,1.5;90,2.0],[0.2,0,0.1]); 

   case 11;
% 3D cylinder. Various elecs at 0, 30, 60, 90 in planes
    el_pos = [0,0.5;30,1;60,1.5;90,2.0];
    el_sz  = [0.2,0,0.1;0.1,0,0.05;0.2,0.2,0.02;0.2,0.4,0.5];
    fmdl= ng_mk_model('cyl',3,el_pos,el_sz); 

   case 0; fmdl = 11; %%%% RETURN MAXIMUM
   otherwise;
     error('huh?')
   end
  
   if exist('img'); fmdl = img; end
