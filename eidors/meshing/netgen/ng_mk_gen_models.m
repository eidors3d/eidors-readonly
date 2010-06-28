function [fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos,  elec_shape);
% NG_MAKE_ELLIP_MODELS: create elliptical models using netgen
%[fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape);
%
% INPUT:
% shape_str = {height, [x_radius, y_radius, [maxsz]]}
%    if height = 0 -> calculate a 2D shape
%    x_radius, y_radius (OPT)  -> elliptical eccentricity in x,y directions(default = 1)
%    maxsz  (OPT)  -> max size of mesh elems (default = course mesh)
%
% ELECTRODE POSITIONS:
%  elec_pos = [x,y,z,nx,ny,nz] centres of each electrode (N_elecs x 6)
%   [x,y,z] is the position, and [nx,ny,nz] is the surface normal
%
% ELECTRODE SHAPES::
%  elec_shape = [width,height, maxsz]  % Rectangular elecs
%     OR
%  elec_shape = [radius, 0, maxsz ]    % Circular elecs
%     OR 
%  elec_shape = [0, sz, maxsz ]         % Point electrodes
%    (point elecs does some tricks with netgen, using sz square, so the elecs aren't exactly where you ask)
%
% Specify either a common electrode shape or for each electrode
%
% OUTPUT:
%  fmdl    - fwd_model object
%  mat_idx - indices of materials
%
% USAGE EXAMPLES:

% (C) Andy Adler, 2010. Licenced under GPL v2 or v3
% $Id$

if isstr(shape_str) && strcmp(shape_str,'UNIT_TEST'); do_unit_test; return; end

if nargin < 4; extra_ng_code = {'',''}; end
cache_obj = { shape_str, elec_pos, elec_shape };

fmdl = eidors_obj('get-cache', cache_obj, 'ng_mk_gen_models' );
if isempty(fmdl);
   fmdl = mk_ellip_model( shape_str, elec_pos, elec_shape);
   eidors_cache('boost_priority', -2); % netgen objs are low priority
   eidors_obj('set-cache', cache_obj, 'ng_mk_gen_models', fmdl);
   eidors_cache('boost_priority', +2); % return values
end

mat_idx = fmdl{2};
fmdl = fmdl{1};

function [fmdl_mat_idx] = mk_ellip_model( shape_str, elec_pos, elec_shape);

   fnstem = tempname;
   geofn= [fnstem,'.geo'];
   ptsfn= [fnstem,'.msz'];
   meshfn= [fnstem,'.vol'];

   is2D = 0;
   [elecs, centres] = parse_elecs( elec_pos, elec_shape, is2D );

   n_pts = write_geo_file(geofn, ptsfn, shape_str, elecs);
   if n_pts == 0 
      call_netgen( geofn, meshfn);
   else
      call_netgen( geofn, meshfn, ptsfn);
   end

   [fmdl,mat_idx] = ng_mk_fwd_model( meshfn, centres, 'ng', []);

   delete(geofn); delete(meshfn); delete(ptsfn); % remove temp files
   if is2D
      [fmdl,mat_idx] = mdl2d_from3d(fmdl,mat_idx);
   end

   % convert CEM to PEM if so configured
   % TODO shunt model is unsupported
   if isfield(fmdl,'electrode');
   fmdl.electrode = pem_from_cem(elecs, fmdl.electrode, fmdl.nodes);
   end

   fmdl_mat_idx = {fmdl,mat_idx};

% for the newest netgen, we can't call msz file unless there are actually points in  it
function n_pts_elecs = write_geo_file(geofn, ptsfn, shape_str, elecs);
   fid=fopen(geofn,'w');
   write_header(fid, shape_str);

   n_elecs = length(elecs);
   %  elecs(i).pos   = [x,y,z]
   %  elecs(i).shape = 'C' or 'R' or 'P'
   %  elecs(i).dims  = [radius] or [width,height]
   %  elecs(i).maxh  = '-maxh=#' or '';
   pts_elecs_idx = []; 

   tank_radius = norm( std( vertcat( elecs(:).pos ), 1), 2);
   for i=1:n_elecs
      name = sprintf('elec%04d',i);
      pos = elecs(i).pos;
      dirn= elecs(i).dirn;
      switch elecs(i).shape
       case 'C'
         write_circ_elec(fid,name, pos, dirn,  ...
               elecs(i).dims, tank_radius, elecs(i).maxh);
       case 'R'
         write_rect_elec(fid,name, pos, dirn,  ...
               elecs(i).dims, tank_radius, elecs(i).maxh);
       case 'P'
         if 0 % Netgen doesn't put elecs where you ask
            pts_elecs_idx = [ pts_elecs_idx, i]; 
            continue; % DON'T print solid cyl
         end
         write_rect_elec(fid,name, pos, dirn,  ...
               elecs(i).dims, tank_radius, elecs(i).maxh);

       otherwise; error('huh? shouldnt get here');
      end
      fprintf(fid,'solid cyl%04d = mainobj    and %s; \n',i,name);
   end

   % SHOULD tank_maxh go here?
   fprintf(fid,'tlo mainobj;\n');
   for i=1:n_elecs
      if any(i == pts_elecs_idx); continue; end
      fprintf(fid,'tlo cyl%04d cyl -col=[1,0,0];\n ',i);
   end

%  if ~isempty(extra_ng_code{1})
%     fprintf(fid,'tlo %s  -col=[0,1,0];\n',extra_ng_code{1});
%  end

   fclose(fid); % geofn
% From Documentation: Syntax is
% np
% x1 y1 z1 h1
% x2 y2 z2 h2
   n_pts_elecs= length(pts_elecs_idx);
   fid=fopen(ptsfn,'w');
   fprintf(fid,'%d\n',n_pts_elecs);
   for i = pts_elecs_idx;
      posxy = elecs(i).pos(1:2);
      fprintf(fid,'%10f %10f 0 %10f\n', posxy, elecs(i).dims(1) );
   end
   fclose(fid); % ptsfn


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
function [elecs, centres] = parse_elecs(elec_pos, elec_shape, is2D );

   n_elecs= size(elec_pos,1); 
   if n_elecs == 0
      elecs= struct([]); % empty
      centres= [];
      return;
   end


   if size(elec_shape,1) == 1
      elec_shape = ones(n_elecs,1) * elec_shape;
   end

   for i= 1:n_elecs
     elecs(i) = elec_spec( elec_shape(i,:), elec_pos(i,:) );
   end
   
   centres = elec_pos(:,1:3);

function elec = elec_spec( row, posrow );
  elec.pos = posrow(1:3);
  elec.dirn= posrow(4:6);

  if row(1) == 0
     elec.shape = 'P' 
     elec.dims  = row(2)*[1,1];
  elseif length(row)<2 || row(2) == 0 % Circular electrodes 
     elec.shape = 'C';
     elec.dims  = row(1);
  elseif row(2)>0      % Rectangular electrodes
     elec.shape = 'R';
     elec.dims  = row(1:2);
  else
     error('negative electrode width');
  end

  if length(row)>=3 && row(3) > 0
     elec.maxh = sprintf('-maxh=%f', row(3));
  else
     elec.maxh = '';
  end


function write_header(fid, shape_str);
   fprintf(fid,'#Automatically generated by ng_mk_gen_models\n');
   fprintf(fid,'algebraic3d\n');
   fprintf(fid,shape_str);

function [mdl2,idx2] = mdl2d_from3d(mdl3,idx3);
   % set name
   mdl2 = eidors_obj('fwd_model',sprintf('%s 2D',mdl3.name));

   % set nodes
   [bdy,idx] = find_boundary(mdl3.elems);
   vtx = mdl3.nodes;
   z_vtx = reshape(vtx(bdy,3), size(bdy) );
   lay0  = find( all(z_vtx==0,2) );
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

   % copy other fields
   if isfield(mdl3,'stimulation'); mdl2.stimulation= mdl3.stimulation; end
   %if isfield(mdl3,'solve');       mdl2.solve = mdl3.solve;            end
   mdl2.solve = 'aa_fwd_solve'; % FIXME? can't use default np_fwd_solve
   if isfield(mdl3,'jacobian');    mdl2.jacobian = mdl3.jacobian;      end
   %if isfield(mdl3,'system_mat');  mdl2.system_mat = mdl3.system_mat;  end
   mdl2.system_mat = 'aa_calc_system_mat'; % FIXME? can't use default np_calc_system_mat

   % update cache
   mdl2 = eidors_obj('fwd_model',mdl2);

function write_rect_elec(fid,name,c, dirn,wh,d,maxh)
% writes the specification for a netgen cuboid on fid, named name, centerd on c,
% in the direction given by vector dirn,
% hw = [height, width]  and depth d
% direction is in the xy plane
   d= min(d);
   w = wh(1); h= wh(2);
   dirn(3) = 0; dirn = dirn/norm(dirn);
   dirnp = [-dirn(2),dirn(1),0];
   dirnp = dirnp/norm(dirnp);

   bl = c - (d/2)* dirn + (w/2)*dirnp - [0,0,h/2];
   tr = c + (d/2)* dirn - (w/2)*dirnp + [0,0,h/2];
   fprintf(fid,'solid %s  = ', name);
   fprintf(fid,' plane (%6.3f,%6.3f,%6.3f;0, 0, -1) and\n', ...
           bl(1),bl(2),bl(3));
   fprintf(fid,' plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f) and\n', ...
           bl(1),bl(2),bl(3),-dirn(1),-dirn(2),0);
   fprintf(fid,' plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f) and\n', ...
           bl(1),bl(2),bl(3),dirnp(1),dirnp(2),0);
   fprintf(fid,' plane(%6.3f,%6.3f,%6.3f;0, 0, 1) and\n', ...
           tr(1),tr(2),tr(3));
   fprintf(fid,' plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f) and\n', ...
           tr(1),tr(2),tr(3),dirn(1),dirn(2),0);
   fprintf(fid,' plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )%s;\n', ...
           tr(1),tr(2),tr(3),-dirnp(1),-dirnp(2),0,maxh);

function write_circ_elec(fid,name,c, dirn,rd,ln,maxh)
% writes the specification for a netgen cylindrical rod on fid,
%  named name, centerd on c,
% in the direction given by vector d, radius rd  lenght ln
% direction is in the xy plane
% the direction vector
   dirn(3) = 0; dirn = dirn/norm(dirn);

   ln = min(ln);
 % I would divide by 2 here (shorted tube in cyl), but ng doesn't like
 % That - it fails for 16 (but no 15 or 17) electrodes
   inpt = c - dirn.*(ln/1);
   outpt =c + dirn.*(ln/1);

   fprintf(fid,'solid %s  = ', name);
   fprintf(fid,'  plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f) and\n', ...
         inpt(1),inpt(2),inpt(3),-dirn(1),-dirn(2),-dirn(3));
   fprintf(fid,'  plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f) and\n', ...
         outpt(1),outpt(2),outpt(3),dirn(1),dirn(2),dirn(3));
   fprintf(fid,'  cylinder(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f;%6.3f) %s;\n', ...
         inpt(1),inpt(2),inpt(3),outpt(1),outpt(2),outpt(3), rd,maxh);


function electrode = pem_from_cem(elecs, electrode, nodes)
% elecs = electrode structure of model, from the parse_elecs function
% electrode = the forward electrode model
% nodes = the coordinates for the nodes
% Can only have one node per electrode so we get a Point Electrode Model.
% Choose the node with the greatest angle, so we atlest pick a consistent
% side of the electrode: NetGen seems to give a random order to the nodes
% in the electrode listing so we can't just pick the first one.
% The nodes aside from those on the edges are not garanteed to be at any
% particular location, so won't be consistent between meshes.
% TODO should probably also adjust contact impedance too: its found later
% by taking the average of the edges around the PEM's node, and those
% will vary for each mesh -- should adjust so all electrodes get a
% consistent effective impedance later.
  Ne = length(electrode);
  for i = 1:Ne
    if elecs(i).shape == 'P'
      % find the angles of the nodes for this electrode relative to (0,0)
      xy = nodes(electrode(i).nodes,:);
      ang = atan2(xy(:,2),xy(:,1));
      % if the angles cover more than 180 degrees, must be an angle
      % roll-over from -pi to +pi, so take all the negative angles
      % and move them up
      if (max(ang) - min(ang)) > pi
        ang = ang + (ang <0)*2*pi;
      end
      % choose the counter-clockwise most node only
      if size(xy,2) == 3 ; ang = ang - xy(:,3); end
      [jnk, ind] = max(ang);
      electrode(i).nodes = electrode(i).nodes(ind);
    end
  end


function do_unit_test
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
              'solid mainobj= plane(0,0,0;0,0,-1)\n' ...
                        'and  plane(0,0,2;0,0,1)\n' ...
                        'and  cyl -maxh=0.2;\n'];
 elec_pos = [  1,  0,  1,   1,  0,  0;
               0,  1,1.2,   0,  1,  0]; 
 elec_shape=[0.1];
 ng_mk_gen_models(shape_str, elec_pos, elec_shape)

