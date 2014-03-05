function [fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos,  elec_shape, elec_obj, extra_ng_code);
% NG_MK_GEN_MODELS: create generic models using netgen
%[fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
%[fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj, extra_ng_code);
%
% INPUT:
% shape_str = string of netgen *.geo file code, with mainobj as the main object
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
% ELECTRODE DEFITIONS:
%  elec_obj = 'obj_name' or {'name','for','each','elec'} where
%    the name is the primitive netgen object (cylinder, plane, etc) which the electrode intersects
%
% OUTPUT:
%  fmdl    - fwd_model object
%  mat_idx - indices of materials
%
% USAGE EXAMPLES:
%shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
%             'solid bottom = plane(0,0,0;0,0,-1);\n' ...
%             'solid top    = plane(0,0,2;0,0,1);\n' ...
%             'solid mainobj= top and bottom and cyl -maxh=0.3;\n'];
%elec_pos = [  1,  0,  1,   1,  0,  0;
%              0,  1,1.2,   0,  1,  0;
%              0.8,  0,  0, 0,  0, -1]; 
%elec_shape=[0.1];
%elec_obj = {'cyl','cyl','bottom'};
%fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

% (C) Andy Adler, 2010. (C) Alistair Boyle, 2013. Licenced under GPL v2 or v3
% $Id$

if isstr(shape_str) && strcmp(shape_str,'UNIT_TEST'); do_unit_test; return; end

if nargin <= 4; extra_ng_code = ''; end
cache_obj = { shape_str, elec_pos, elec_shape, elec_obj, extra_ng_code };

fmdl = eidors_obj('get-cache', cache_obj, 'ng_mk_gen_models' );
if isempty(fmdl);
   fmdl = mk_gen_model( shape_str, elec_pos, elec_shape, elec_obj, extra_ng_code);
   eidors_cache('boost_priority', -2); % netgen objs are low priority
   eidors_obj('set-cache', cache_obj, 'ng_mk_gen_models', fmdl);
   eidors_cache('boost_priority', +2); % return values
end

mat_idx = fmdl.mat_idx_reordered;

function fmdl = mk_gen_model( shape_str, elec_pos, elec_shape, ...
                         elec_obj, extra_ng_code);

   fnstem = tempname;
   geofn= [fnstem,'.geo'];
   ptsfn= [fnstem,'.msz'];
   meshfn= [fnstem,'.vol'];

   is2D = 0;
   [elecs, centres] = parse_elecs( elec_pos, elec_shape, is2D, elec_obj );

   n_pts = write_geo_file(geofn, ptsfn, shape_str, elecs, extra_ng_code);
   if n_pts == 0 
      call_netgen( geofn, meshfn);
   else
      call_netgen( geofn, meshfn, ptsfn);
   end

   fmdl = ng_mk_fwd_model( meshfn, centres, 'ng', []);

   delete(geofn); delete(meshfn); delete(ptsfn); % remove temp files
   if is2D
      fmdl = mdl2d_from3d(fmdl);
   end

   % convert CEM to PEM if so configured
   % TODO shunt model is unsupported
   if isfield(fmdl,'electrode');
     fmdl.electrode = pem_from_cem(elecs, fmdl.electrode, fmdl.nodes);
   end

% for the newest netgen, we can't call msz file unless there are actually points in  it
function n_pts_elecs = write_geo_file(geofn, ptsfn, shape_str, ...
                          elecs, extra_ng_code);
   fid=fopen(geofn,'w');
   write_header(fid, shape_str);

   n_elecs = length(elecs);
   %  elecs(i).pos   = [x,y,z]
   %  elecs(i).shape = 'C' or 'R' or 'P'
   %  elecs(i).dims  = [radius] or [width,height]
   %  elecs(i).maxh  = '-maxh=#' or '';
   pts_elecs_idx = []; 

   if n_elecs > 1
      elec_depth = min(nonzeros(distmat(vertcat( elecs(:).pos ))))/2;
      % tank_radius = norm( std( vertcat( elecs(:).pos ), 1), 2);
      % NOTE: all functions but the point electrode used to use
      % tank_radius/4. Point electrode used just tank_radius.
   end
   for i=1:n_elecs
      name = sprintf('elec%04d',i);
      pos = elecs(i).pos;
      dirn= elecs(i).dirn;
      switch elecs(i).shape
       case 'C'
         write_circ_elec(fid,name, pos, dirn,  ...
               elecs(i).dims, elec_depth, elecs(i).maxh);
       case 'R'
         write_rect_elec(fid,name, pos, dirn,  ...
               elecs(i).dims, elec_depth, elecs(i).maxh);
       case 'P'
         if 0 % Netgen doesn't put elecs where you ask
            pts_elecs_idx = [ pts_elecs_idx, i]; 
            continue; % DON'T print solid cyl
         end
         write_rect_elec(fid,name, pos, dirn,  ...
               elecs(i).dims, elec_depth, elecs(i).maxh);

       case 'U'
         eidors_msg('user defined electrode %d',i, 4);
         continue;
 
       otherwise; error('huh? shouldnt get here');
      end
      fprintf(fid,'solid cyl%04d = mainobj    and %s; \n',i,name);
   end

   % In some cases, ng can't handle an extra ';'
   if ~isempty(extra_ng_code)
      fprintf(fid,'%s;\n', extra_ng_code);
   end
   % SHOULD tank_maxh go here?
   fprintf(fid,'tlo mainobj;\n');
   for i=1:n_elecs
      if any(i == pts_elecs_idx); continue; end
      if elecs(i).shape == 'U';   continue; end % USER WILL DEFINE
      fprintf(fid,'tlo cyl%04d %s -col=[1,0,0];\n',i, elecs(i).ngobj);
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
function [elecs, centres] = parse_elecs( elec_pos, elec_shape, is2D, elec_obj );

   n_elecs= size(elec_pos,1); 
   if n_elecs == 0
      elecs= struct([]); % empty
      centres= [];
      return;
   end


   if size(elec_shape,1) == 1
      elec_shape = ones(n_elecs,1) * elec_shape;
   end
   
   centres = elec_pos(:,1:3);

   for i= 1:n_elecs
     if isstr(elec_obj); 
        elecs(i) = elec_spec( elec_shape(i,:), elec_pos(i,:), elec_obj );
     else
        elecs(i) = elec_spec( elec_shape(i,:), elec_pos(i,:), elec_obj{i}  );
     end
   end

function elec = elec_spec( row, posrow, elec_obj );
  elec.pos =  posrow(1:3);
  elec.dirn=  posrow(4:6);

  if all(isnan(elec.dirn))
     elec.shape = 'U'; %user defined
     elec.dims = NaN;
     elec.ngobj = '';
     elec.maxh = '';
     elec = orderfields(elec); % UNBELIEVABLY STUPID MATLAB
     return
  end

  elec.ngobj= elec_obj;

  if row(1) == 0
     elec.shape = 'P';
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

     elec = orderfields(elec); % UNBELIEVABLY STUPID MATLAB

function write_header(fid, shape_str);
   fprintf(fid,'#Automatically generated by ng_mk_gen_models\n');
   fprintf(fid,'algebraic3d\n');
   fprintf(fid,shape_str);

function write_rect_elec(fid,name,c, dirn,wh,d,maxh)
% writes the specification for a netgen cuboid on fid, named name, centerd on c,
% in the direction given by vector dirn,
% hw = [height, width]  and depth d
% direction is in the xy plane
   d= min(d);
   w = wh(1); h= wh(2);
   dirn = dirn/norm(dirn);
   dirnp = [-dirn(2),dirn(1),0];
if any(isnan(dirnp));
   error('ng_mk_gen_models: how to define width and height for vertical electrodes?')
end
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
   dirn = dirn/norm(dirn);

   ln = min(ln);
 % I would divide by 2 here (shorted tube in cyl), but ng doesn't like
 % That - it fails for 16 (but no 15 or 17) electrodes
   inpt = c - dirn.*(ln/1);
   outpt =c + dirn.*(ln/1);

   fprintf(fid,'solid %s  = ', name);
   fprintf(fid,'  plane(%f,%f,%f;%f,%f,%f) and\n',       inpt, -dirn);
   fprintf(fid,'  plane(%f,%f,%f;%f,%f,%f) and\n',       outpt, dirn);
   fprintf(fid,'  cylinder(%f,%f,%f;%f,%f,%f;%f) %s;\n', inpt, outpt, rd,maxh);


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
  for tn = 14 %1:14
     eidors_msg('ng_mk_gen_models: unit_test %02d',tn,1);
     fmdl= do_test_number(tn);
     show_fem(fmdl); drawnow
  end

function fmdl= do_test_number(tn)
switch tn
   case 1;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
              'solid mainobj= orthobrick(-2,-2,0;2,2,2) and cyl -maxh=0.3;\n'];
 elec_pos = []; elec_shape = []; elec_obj = {};
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 2;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
              'solid mainobj= plane(0,0,0;0,0,-1)\n' ...
                        'and  plane(0,0,2;0,0,1)\n' ...
                        'and  cyl -maxh=0.3;\n'];
 elec_pos = [  1,  0,  1,   1,  0,  0;
               0,  1,1.2,   0,  1,  0]; 
 elec_shape=[0.1];
 elec_obj = 'cyl';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 3;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
              'solid mainobj= orthobrick(-2,-2,0;2,2,2) and cyl -maxh=0.3;\n'];
 th = linspace(0,2*pi,15)'; th(end) = [];
 cs = [cos(th), sin(th)];
 elec_pos = [  cs, th/2/pi + 0.5, cs, 0*th];
 elec_shape=[0.1];
 elec_obj = 'cyl';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 4;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
              'solid mainobj= orthobrick(-2,-2,0;2,2,2) and cyl -maxh=0.3;\n'];
 th = linspace(0,2*pi,15)'; th(end) = [];
 cs = [cos(th), sin(th)];
 elec_pos = [  cs, th/2/pi + 0.5, cs, 0*th];
 elec_shape=[0.1*th/2/pi + 0.05];
 elec_obj = 'cyl';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 5;
 shape_str = ['solid cyl    = cylinder (0,0,0; 1,0,0; 1); \n', ...
              'solid mainobj= plane(0,0,0;-1,0,0)\n' ...
                        'and  plane(2,0,0;1,0,0)\n' ...
                        'and  cyl -maxh=0.3;\n'];
 elec_pos = [  1,  0,  1,   0,  0,  1;
             1.2,  1,  0,   0,  1,  0]; 
 elec_shape=[0.1];
 elec_obj = 'cyl';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 6;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
              'solid bottom = plane(0,0,0;0,0,-1);\n' ...
              'solid top    = plane(0,0,2;0,0,1);\n' ...
              'solid mainobj= top and bottom and cyl -maxh=0.3;\n'];
 elec_pos = [  1,  0,  1,   1,  0,  0;
               0,  1,1.2,   0,  1,  0;
               0.8,  0,  0, 0,  0, -1]; 
 elec_shape=[0.1];
 elec_obj = {'cyl','cyl','bottom'};
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 7;
 shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
              'solid mainobj= top and orthobrick(-2,-2,-2;2,2,0);\n'];
 elec_pos = [  1,  0,  0,   0,  0,  1;
               0,  0,  0,   0,  0,  1;
              -1,  0,  0,   0,  0,  1];
 elec_shape=[0.1];
 elec_obj = 'top';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 8;
 shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
              'solid cyl    = ellipticcylinder(0,0,0;2.5,0,0;0,1,0);\n' ...
              'solid mainobj= top and cyl and orthobrick(-2,-2,-2;2,2,0);\n'];
 elec_pos = [  1,  0,  0,   0,  0,  1;
               0,  0,  0,   0,  0,  1;
              -1,  0,  0,   0,  0,  1;
               1, -1,-1.2,  0, -1,  0;
               0, -1,-1.0,  0, -1,  0;
              -1, -1,-0.8,  0, -1,  0];
 elec_shape=[0.1];
 elec_obj = {'top','top','top','cyl','cyl','cyl'};
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 9;
 shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
              'solid ball   = sphere(-1.25,0,-1;0.5); tlo ball;\n' ...
              'solid mainobj= top and orthobrick(-2,-1,-2;2,1,0) and not ball -maxh=0.5;\n'];
 elec_pos = linspace( -1.5,1.5,5)';
 elec_pos = [  elec_pos, elec_pos*[0,0,0,0], elec_pos*0+1];
 elec_shape=[0.3];
 elec_obj = 'top';
 [fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
 img = mk_image( fmdl, 1);
 img.elem_data(mat_idx{2}) = 1.1; 
 
 fmdl = img; % so that the code shows the image

  case 10;
 shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
              'solid mainobj= top and orthobrick(-3,-3,-2;3,3,0) -maxh=0.5;\n'];
 [elec_pos_x,elec_pos_y] = meshgrid(linspace( -1.5,1.5,5),linspace(-2,2,7));
 elec_pos = [  elec_pos_x(:), elec_pos_y(:), ones(size(elec_pos_x(:)))*[0,0,0,1] ];
 elec_shape=[0.2];
 elec_obj = 'top';
 [fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 11;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1.0); \n', ...
              'solid tank   = orthobrick(-2,-2,0;2,2,0.4) and cyl; \n', ...
              'solid fish   = ellipsoid(0.2,0.2,0.2;0.2,0,0;0,0.1,0;0,0,0.1); tlo fish;\n', ...
              'solid mainobj= tank and not fish -maxh=0.3;\n'];
 n_elec = 7;
 th = linspace(0,2*pi,n_elec+1)'; th(end) = [];
 cs = [cos(th), sin(th)];
 elec_pos = [  cs, 0.2+0*th, cs, 0*th];
 elec_shape=[0.05];
 for i=1:n_elec; elec_obj{i} = 'cyl'; end
 i=i+1;elec_pos(i,:) = [ 0  ,0.2,0.2,-1,0,0]; elec_obj{i} = 'fish';
 i=i+1;elec_pos(i,:) = [ 0.4,0.2,0.2, 1,0,0]; elec_obj{i} = 'fish';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 12;
shape_str = ['solid top     = ellipsoid(0,0,0; 0,0,1; 1,0,0; 0,1,0); \n' ...
    'solid mainobj= top and orthobrick(-2,-2,0;2,2,2) -maxh=0.1;\n'];
deg2rad = pi/180;
th = (-70:20:70)'*deg2rad;
 elec_pos = [0*th,sin(th),cos(th),0*th,sin(th),cos(th); ...
             sin(th),0*th,cos(th),sin(th),0*th,cos(th)];
 elec_shape=[0.05];
 elec_obj = 'top';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 13;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,1,0; 1); \n', ...
              'solid bottom = plane(0, 0,0;0,-1,0);\n' ...
              'solid top    = plane(0,10,0;0, 1,0);\n' ...
              'solid cut1   = plane(0, 4,0;0,-1,0);\n' ...
              'solid cut2   = plane(0, 6,0;0, 1,0);\n' ...
              'solid ball   = cyl and cut1 and cut2;  tlo ball;\n' ...
              'solid mainobj= ( top and (not cut2) and cyl ) or ' ...
                      '(bottom      and (not cut1) and cyl ) -maxh=0.8;\n'];
 elec_pos = [ 0, 10,  0, 0,  1,  0; 
              0,  0,  0, 0, -1,  0]; 
 elec_shape=[1.0];
 elec_obj = {'top','bottom'};
 [fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
 fmdl = mk_image(fmdl,1); 
 fmdl.elem_data(mat_idx{2}) = 1.1;

   case 14;
 shape_str = [ ...
  'solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
  'solid mainobj= plane(0,0,0;0,0,-1)\n' ...
        'and  plane(0,0,2;0,0,1)\n' ...
        'and  cyl -maxh=0.3;\n' ...
  'solid in_elec = sphere(0,-1,1;0.2)' ...
        'and not    sphere(0,-1,1;0.15) -maxh=0.05;\n' ...
        'solid in_elec0= in_elec  and mainobj;\n' ...
        'tlo in_elec0 cyl;\n' ...
  'solid out_elec = sphere(0,-1,1;0.4)' ...
        'and not    sphere(0,-1,1;0.35) -maxh=0.05;\n' ...
        'solid out_elec0= out_elec  and mainobj;\n' ...
        'tlo out_elec0 cyl;\n'];
 elec_pos = [  0, -1,   0, NaN,NaN,NaN; % get rest of tank first, to remove
               1,  0,   1,   1,  0,  0;
               0,  1, 1.2,   0,  1,  0;
               0, -1, 1.2, NaN,NaN,NaN;
               0, -1, 1.4, NaN,NaN,NaN];
 elec_shape=[0.1];
 elec_obj = 'cyl';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
 fmdl.electrode = fmdl.electrode(2:end); % Throw away the first
  otherwise;
     error('huh?')
end
