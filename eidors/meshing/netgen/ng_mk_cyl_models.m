function fmdl = ng_make_cyl_models(cyl_shape, elec_pos, ...
                  elec_shape, extra_ng_code);
% NG_MAKE_CYL_MODELS: create cylindrical models using netgen
% fmdl = ng_make_cyl_models(cyl_shape, elec_pos, elec_shape, extra_ng_code);
% 
% cyl_shape = {height, [radius, [maxsz]]}
%    if height = 0 -> calculate a 2D shape
%    radius (OPT)  -> (default = 1)
%    maxsz  (OPT)  -> max size of mesh elems (default = courase mesh)
%
% ELECTRODE POSITIONS:
%  elec_pos = [n_elecs_per_plane,z_planes] 
%     OR
%  elec_pos = [x,y,z] centres of each electrode (N_elecs x 3)
%
% ELECTRODE SHAPES::
%  elec_shape = [width,height, {maxsz}]  % Rectangular elecs
%     OR
%  elec_shape = [radius, {0, maxsz} ]  % Circular elecs
%     maxsz  (OPT)  -> max size of mesh elems (default = courase mesh)
%     maxsz  (OPT)  -> max size of mesh elems (default = courase mesh)
%
% Specify either a common electrode shape or for each electrode
%
%
% USAGE EXAMPLES:
% Simple 3D cylinder. Radius = 1. No electrodes
%   fmdl= ng_make_cyl_models(3,[0],[]); 
% Simple 2D cylinder. Radius = 2. Set minsize to refine
%   fmdl= ng_make_cyl_models([0,2,.2],[0],[]); 
% 3D cylinder. Radius = 1. 2 planes of 8 elecs with radius 0.1
%   fmdl= ng_make_cyl_models(3,[8,1,2],[0.1]); 
% 3D cylinder. Radius = 1. 6 electrodes with elec refinement
%   fmdl= ng_make_cyl_models(3,[7,1],[0.2,0,0.05]); 

% (C) Andy Adler, 2009. Licenced under GPL v2 or v3
% $Id$

fnstem = tempname;
geofn= [fnstem,'.geo'];
meshfn= [fnstem,'.vol'];

[tank_height, tank_radius, tank_maxh, is2D] = parse_shape(cyl_shape);
[elecs, centres] = parse_elecs( elec_pos, elec_shape,  ...
                       tank_height, tank_radius );

write_geo_file(geofn, tank_height, tank_radius, tank_maxh, elecs);
call_netgen( geofn, meshfn);

fmdl = ng_mk_fwd_model( meshfn, centres, 'ng', []);
%delete(geofn); delete(meshfn); % remove temp files
if is2D
   fmdl = mdl2d_from3d(fmdl);
end

function write_geo_file(geofn, tank_height, tank_radius, tank_maxh, elecs)
   fid=fopen(geofn,'w');
   write_header(fid,tank_height,tank_radius,tank_maxh);

   n_elecs = length(elecs);
   %  elecs(i).pos   = [x,y,z]
   %  elecs(i).shape = 'C' or 'R'
   %  elecs(i).dims  = [radius] or [width,height]
   %  elecs(i).maxh  = '-maxh=#' or '';
   for i=1:n_elecs
      name = sprintf('elec%04d',i);
      pos = elecs(i).pos;
      if elecs(i).shape == 'C'
         write_circ_elec(fid,name, pos, pos,  ...
               elecs(i).dims, tank_radius, elecs(i).maxh);
      else
   write_rect_elec(fid,name,c, dirn,hw,d,maxh)
      end

      fprintf(fid,'solid cyl%04d = bigcyl    and %s; \n',i,name);
   end

   % SHOULD tank_maxh go here?
   fprintf(fid,'tlo bigcyl;\n');
   for i=1:n_elecs
      fprintf(fid,'tlo cyl%04d cyl -col=[1,0,0];\n ',i);
   end

   fclose(fid);

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
%  elec_pos = [x,y,z] centres of each electrode (N_elecs x 3)
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
function [elecs, centres] = parse_elecs( elec_pos, elec_shape, hig, rad );
   % It never makes sense to specify only one elec
   if size(elec_pos,1) == 1
       % Parse elec_pos = [n_elecs_per_plane,z_planes] 
      n_elecs= elec_pos(1); % per plane
      th = linspace(0,2*pi, n_elecs+1)'; th(end)=[];
      csth= rad*[sin(th),cos(th)]; % clockwise from TDC

      centres=[];
      on_elecs = ones(n_elecs, 1);
      for i=2:length(elec_pos)
        centres= [centres; [csth, on_elecs*elec_pos(i)]];
      end
   else
      centres = elec_pos;
   end
   n_elecs= size(centres,1);

   if size(elec_shape,1) == 1
      elec_shape = ones(n_elecs,1) * elec_shape;
   end
   elec_shape = [elec_shape, zeros(n_elecs, 2)]; % add default zeros

   elecs= struct; % empty
   for i= 1:n_elecs
     row = elec_shape(i,:);

     elecs(i).pos  = centres(i,:);
     if row(2) == 0
        elecs(i).shape = 'C';
        elecs(i).dims  = row(1);
     else
        elecs(i).shape = 'R';
        elecs(i).dims  = row(1:2);
     end

     if row(3) == 0
        elecs(i).maxh = '';
     else
        elecs(i).maxh = sprintf('-maxh=%f', row(3));
     end
   end
   


function write_header(fid,tank_height,tank_radius,maxsz);
   if maxsz==0; 
      maxsz = '';
   else
      maxsz = sprintf('-maxh=%f',maxsz);
   end

   fprintf(fid,'#Automatically generated by ng_make_cyl_models\n');
   fprintf(fid,'algebraic3d\n');
   fprintf(fid,'solid cyl=cylinder (0,0,0;0,0,%6.2f;%6.2f); \n', ...
           tank_height, tank_radius);
   fprintf(fid,['solid bigcyl= plane(0,0,0;0,0,-1)\n' ...
                'and  plane(0,0,%6.2f;0,0,1)\n' ...
                'and  cyl %s;\n'],tank_height,maxsz);  

function mdl2 = mdl2d_from3d(mdl3)
   mdl2 = eidors_obj('2D','fwd_model');
   bdy = find_boundary(mdl3.elems);
   vtx = mdl3.nodes;
   z_vtx = reshape(vtx(bdy,3), size(bdy) );
   lay0  = find( all(z_vtx==0,2) );
   bdy0  = bdy( lay0, :);
   
   vtx0  = unique(bdy0(:));
   mdl2.nodes = vtx(vtx0,1:2);

   nmap  = zeros(size(vtx,1),1); nmap(vtx0) = 1:length(vtx0);
   bdy0  = reshape(nmap(bdy0), size(bdy0) ); % renumber to new scheme
   mdl2.elems = bdy0;

   mdl2.gnd_node = nmap(mdl3.gnd_node);

% Manage Electrodes
   if ~isfield(mdl3,'electrode'); return; end

   mdl2.electrode = mdl3.electrode;
   for i=1:length(mdl2.electrode);
      enodes = nmap( mdl2.electrode(i).nodes );
      enodes(enodes==0) = []; % Remove 3D layers
      mdl2.electrode(i).nodes = enodes;
   end


function write_rect_elec(fid,name,c, dirn,hw,d,maxh)
% writes the specification for a netgen cuboid on fid, named name, centerd on c,
% in the direction given by vector dirn,
% hw = [height, width]  and depth d
% direction is in the xy plane
   h = hw(1); w= hw(2);
   dirnp = [-dirn(2),dirn(1),0];
   dirnp = dirnp/norm(dirnp);

   bl = c - (d/2)* dirn + (w/2)*dirnp - [0,0,h/2];
   tr = c + (d/2)* dirn - (w/2)*dirnp + [0,0,h/2];
   fprintf(fid,'solid %s  =plane (%6.3f,%6.3f,%6.3f;0, 0, -1 )\n ', ...
           name,bl(1),bl(2),bl(3));
   fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )\n ', ...
           bl(1),bl(2),bl(3),-dirn(1),-dirn(2),0);
   fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )\n ', ...
           bl(1),bl(2),bl(3),dirnp(1),dirnp(2),0);
   fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;0, 0, 1  )\n ', ...
           tr(1),tr(2),tr(3));
   fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )\n ', ...
           tr(1),tr(2),tr(3),dirn(1),dirn(2),0);
   fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )%s;\n ', ...
           tr(1),tr(2),tr(3),-dirnp(1),-dirnp(2),0,maxh);

function write_circ_elec(fid,name,c, dirn,rd,ln,maxh)
% writes the specification for a netgen cylindrical rod on fid,
%  named name, centerd on c,
% in the direction given by vector d, radius rd  lenght ln
% direction is in the xy plane
% the direction vector
   dirn(3) = 0; dirn = dirn/norm(dirn);

   inpt = c - dirn.*(ln/2);
   outpt =c + dirn.*(ln/2);

   fprintf(fid,'solid %s  = ', name);
   fprintf(fid,'  plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f) and\n', ...
         inpt(1),inpt(2),inpt(3),-dirn(1),-dirn(2),-dirn(3));
   fprintf(fid,'  plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f) and\n', ...
         outpt(1),outpt(2),outpt(3),dirn(1),dirn(2),dirn(3));
   fprintf(fid,' cylinder(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f;%6.3f) %s;\n', ...
         inpt(1),inpt(2),inpt(3),outpt(1),outpt(2),outpt(3), rd,maxh);
