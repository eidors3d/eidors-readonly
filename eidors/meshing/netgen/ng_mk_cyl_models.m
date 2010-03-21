function [fmdl,mat_idx] = ng_mk_cyl_models(cyl_shape, elec_pos, ...
                  elec_shape, extra_ng_code);
% NG_MAKE_CYL_MODELS: create cylindrical models using netgen
%[fmdl,mat_idx] = ng_mk_cyl_models(cyl_shape, elec_pos, ...
%                 elec_shape, extra_ng_code);
% INPUT:
% cyl_shape = {height, [radius, [maxsz]]}
%    if height = 0 -> calculate a 2D shape
%    radius (OPT)  -> (default = 1)
%    maxsz  (OPT)  -> max size of mesh elems (default = course mesh)
%
% ELECTRODE POSITIONS:
%  elec_pos = [n_elecs_per_plane,z_planes] 
%     OR
%  elec_pos = [degrees,z] centres of each electrode (N_elecs x 2)
%
% ELECTRODE SHAPES::
%  elec_shape = [width,height, {maxsz, pem}]  % Rectangular elecs
%     OR
%  elec_shape = [radius, {0, maxsz, pem} ]  % Circular elecs
%     maxsz  (OPT)  -> max size of mesh elems (default = course mesh)
%     pem  (OPT)  -> 1: Point Electrode Model, 0: Complete Electrode Model (default)
%
% Specify either a common electrode shape or for each electrode
%
% EXTRA_NG_CODE
%   string of extra code to put into netgen geo file. Normally this
%   would be to insert extra materials into the space
%
% OUTPUT:
%  fmdl    - fwd_model object
%  mat_idx - indices of materials (if extra_ng_code is used)
%    Note mat_idx does not work in 2D. Netgen does not provide it.
%
%
% USAGE EXAMPLES:
% Simple 3D cylinder. Radius = 1. No electrodes
%   fmdl= ng_mk_cyl_models(3,[0],[]); 
% Simple 2D cylinder. Radius = 2. Set minsize to refine
%   fmdl= ng_mk_cyl_models([0,2,.2],[0],[]); 
% 3D cylinder. Radius = 1. 2 planes of 8 elecs with radius 0.1
%   fmdl= ng_mk_cyl_models(3,[8,1,2],[0.1]); 
% 3D cylinder. Radius = 1. 6 circ elecs with elec refinement
%   fmdl= ng_mk_cyl_models(3,[7,1],[0.2,0,0.05]); 
% 3D cylinder. Radius = 1. 7 rect elecs with no refinement
%   fmdl= ng_mk_cyl_models(3,[7,1],[0.2,0.3]); 
% 2D cylinder. Radius = 1. 11 rect elecs with refinement
%   fmdl= ng_mk_cyl_models(0,[11],[0.2,0,0.05]); 
% 2D cylinder. Radius = 1.5. Refined(0.1). 11 elecs with refinement
%   fmdl= ng_mk_cyl_models([0,1,0.1],[11],[0.2,0,0.02]); 
% 2D cylinder. elecs at 0, 90 and 120 degrees
%   fmdl= ng_mk_cyl_models(0,[0;90;120],[0.2,0,0.03]); 
% 2D cylinder. elecs at 0 (large,refined) and 120 (small) degrees
%   fmdl= ng_mk_cyl_models(0,[0;120],[0.4,0,0.01;0.1,0,0.1]); 
% 3D cylinder. elecs at 0, 30, 60, 90 in planes
%   fmdl= ng_mk_cyl_models(3,[0,0.5;30,1;60,1.5;90,2.0],[0.2,0,0.1]); 
% 3D cylinder. Various elecs at 0, 30, 60, 90 in planes
%   el_pos = [0,0.5;30,1;60,1.5;90,2.0];
%   el_sz  = [0.2,0,0.1;0.1,0,0.05;0.2,0.2,0.02;0.2,0.4,0.5];
%   fmdl= ng_mk_cyl_models(3,el_pos,el_sz); 
% Simple 3D cylinder with a ball
%   extra={'ball','solid ball = sphere(0.5,0.5,2;0.4);'}
%   [fmdl,mat_idx]= ng_mk_cyl_models(3,[0],[],extra); 
%   img= eidors_obj('image','ball'); img.fwd_model= fmdl;
%   img.elem_data(mat_idx{1}) = 1; img.elem_data(mat_idx{2}) = 2;
% 3D cylinder with 8 electrodes and cube
%   extra={'cube','solid cube = orthobrick(0.5,0.5,0.5;0,0,1.5);'}
%   [fmdl,mat_idx]= ng_mk_cyl_models(2,[8,0.5,1.5],[0.1],extra); 
% 3D cylinder with inner cylinder
%   extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,1;1,1,2) -maxh=0.05;'}
%   [fmdl,mat_idx]= ng_mk_cyl_models(3,[0],[],extra); 
% 2D cylinder with 8 electrodes and hole
%   extra={'ball','solid ball = sphere(0.2,0.2,0;0.2) -maxh=0.05;'}
%   fmdl= ng_mk_cyl_models(0,[8],[0.1,0,0.05],extra); 
% 2D cylinder with 9 electrodes and inner cylinder
%   extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.03;'}
%   fmdl= ng_mk_cyl_models(0,[9],[0.2,0,0.05],extra); 
%   img= eidors_obj('image','ball'); img.fwd_model= fmdl;
%   ctr = interp_mesh(fmdl); ctr=(ctr(:,1)-0.2).^2 + (ctr(:,2)-0.2).^2;
%   img.elem_data = 1 + 0.1*(ctr<0.2^2);

% (C) Andy Adler, 2009. Licenced under GPL v2 or v3
% $Id$

if nargin < 4; extra_ng_code = {'',''}; end
cache_obj = { cyl_shape, elec_pos, elec_shape, extra_ng_code};

fmdl = eidors_obj('get-cache', cache_obj, 'ng_mk_cyl_models' );
if isempty(fmdl);
   fmdl = mk_cyl_model( cyl_shape, elec_pos, elec_shape, extra_ng_code );
   eidors_cache('boost_priority', -2); % netgen objs are low priority
   eidors_obj('set-cache', cache_obj, 'ng_mk_cyl_models', fmdl);
   eidors_cache('boost_priority', +2); % return values
end

function fmdl = mk_cyl_model( cyl_shape, elec_pos, elec_shape, extra_ng_code );

   fnstem = tempname;
   geofn= [fnstem,'.geo'];
   ptsfn= [fnstem,'.msz'];
   meshfn= [fnstem,'.vol'];

   [tank_height, tank_radius, tank_maxh, is2D] = parse_shape(cyl_shape);
   [elecs, centres] = parse_elecs( elec_pos, elec_shape,  ...
                          tank_height, tank_radius, is2D );

   n_pts = write_geo_file(geofn, ptsfn, tank_height, tank_radius, ...
                  tank_maxh, elecs, extra_ng_code);
   if n_pts == 0 
      call_netgen( geofn, meshfn);
   else
      call_netgen( geofn, meshfn, ptsfn);
   end

   [fmdl,mat_idx] = ng_mk_fwd_model( meshfn, centres, 'ng', []);

   delete(geofn); delete(meshfn); delete(ptsfn); % remove temp files
   if is2D
      [fmdl,max_idx] = mdl2d_from3d(fmdl,mat_idx);
   end

   % convert CEM to PEM if so configured
   % TODO shunt model is unsupported
   if isfield(fmdl,'electrode');
   fmdl.electrode = pem_from_cem(elecs, fmdl.electrode, fmdl.nodes);
   end


% for the newest netgen, we can't call msz file unless there are actually points in  it
function n_pts_elecs = write_geo_file(geofn, ptsfn, tank_height, tank_radius, ...
                        tank_maxh, elecs, extra_ng_code);
   fid=fopen(geofn,'w');
   write_header(fid,tank_height,tank_radius,tank_maxh,extra_ng_code);

   n_elecs = length(elecs);
   %  elecs(i).pos   = [x,y,z]
   %  elecs(i).shape = 'C' or 'R'
   %  elecs(i).dims  = [radius] or [width,height]
   %  elecs(i).maxh  = '-maxh=#' or '';
   pts_elecs_idx = []; 
%^keyboard
   for i=1:n_elecs
      name = sprintf('elec%04d',i);
      pos = elecs(i).pos;
      switch elecs(i).shape
       case 'C'
         write_circ_elec(fid,name, pos, pos,  ...
               elecs(i).dims, tank_radius, elecs(i).maxh);
       case 'R'
         write_rect_elec(fid,name, pos, pos,  ...
               elecs(i).dims, tank_radius, elecs(i).maxh);
       case 'P'
         pts_elecs_idx = [ pts_elecs_idx, i]; 
         continue; % DON'T print solid cyl

       otherwise; error('huh? shouldnt get here');
      end
      fprintf(fid,'solid cyl%04d = bigcyl    and %s; \n',i,name);
   end

   % SHOULD tank_maxh go here?
   fprintf(fid,'tlo bigcyl;\n');
   for i=1:n_elecs
      if any(i == pts_elecs_idx); continue; end
      fprintf(fid,'tlo cyl%04d cyl -col=[1,0,0];\n ',i);
   end

   if ~isempty(extra_ng_code{1})
      fprintf(fid,'tlo %s  -col=[0,1,0];\n',extra_ng_code{1});
   end

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
      fprintf(fid,'%10f %10f 0 %10f\n', posxy, elecs(i).dims );
   end
   fclose(fid); % ptsfn

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

   if is2D
      elec_pos(:,2) = hig/2;
   end

   % It never makes sense to specify only one elec
   % So elec_pos means the number of electrodes in this case
   if size(elec_pos,1) == 1
       % Parse elec_pos = [n_elecs_per_plane,z_planes] 
      n_elecs= elec_pos(1); % per plane
      th = linspace(0,2*pi, n_elecs+1)'; th(end)=[];
      csth= rad*[sin(th),cos(th)]; % clockwise from TDC

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
   
   centres = [rad*sin(el_th),rad*cos(el_th),el_z];
   for i= 1:n_elecs; elecs(i).pos  = centres(i,:); end

   if n_elecs == 0
      elecs= struct([]); % empty
   end

function elec = elec_spec( row, is2D, hig, rad )
  if     is2D
     if length(row)>=2 && row(2) == -1 % Point electrodes
        % Create rectangular electrodes with bottom, cw point where we want
        elec.shape = 'P' 
        if length(row)>=3 && row(3) > 0
           elec.dims  =  row(3);
        else
           elec.dims  =  rad; % Make big if unspecified
        end
     else
        elec.shape = 'R';
        elec.dims  = [row(1),hig];
     end
  else
     if length(row)<2 || row(2) == 0 % Circular electrodes 
        elec.shape = 'C';
        elec.dims  = row(1);
     elseif row(2) == -1 % Point electrodes
        % Create rectangular electrodes with bottom, cw point where we want
        elec.shape = 'P' 
        if length(row)>=3 && row(3) > 0
           elec.dims  =  row(3);
        else
           elec.dims  =  rad; % Make big if unspecified
        end
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

  if length(row)<4 || row(4) == 0
     elec.model = 'cem'; % Complete Electrode Model (CEM)
  else
     elec.model = 'pem'; % Point Electrode Model (PEM)
  end
  %TODO support Shunt Electrode Model (SEM)

function write_header(fid,tank_height,tank_radius,maxsz,extra);
   if maxsz==0; 
      maxsz = '';
   else
      maxsz = sprintf('-maxh=%f',maxsz);
   end

   if ~isempty( extra{1} )
      extra{1} = [' and not ',extra{1}];
   end

   fprintf(fid,'#Automatically generated by ng_mk_cyl_models\n');
   fprintf(fid,'algebraic3d\n');
   fprintf(fid,'%s\n',extra{2}); % Define extra stuff here
   fprintf(fid,'solid cyl=cylinder (0,0,0;0,0,%6.2f;%6.2f); \n', ...
           tank_height, tank_radius);
   fprintf(fid,['solid bigcyl= plane(0,0,0;0,0,-1)\n' ...
                'and  plane(0,0,%6.2f;0,0,1)\n' ...
                'and  cyl %s %s;\n'],tank_height,extra{1},maxsz);  

function [mdl2,idx2] = mdl2d_from3d(mdl3,idx3);
   % set name
   mdl2 = eidors_obj('fwd_model',sprintf('%s 2D',mdl3.name));

   % set nodes
   bdy = find_boundary(mdl3.elems);
   vtx = mdl3.nodes;
   z_vtx = reshape(vtx(bdy,3), size(bdy) );
   lay0  = find( all(z_vtx==0,2) );
   bdy0  = bdy( lay0, :);
   
   vtx0  = unique(bdy0(:));
   mdl2.nodes = vtx(vtx0,1:2);

   % TODO set boundary

   % set elems
   nmap  = zeros(size(vtx,1),1); nmap(vtx0) = 1:length(vtx0);
   bdy0  = reshape(nmap(bdy0), size(bdy0) ); % renumber to new scheme
   mdl2.elems = bdy0;

   % set boundary
   mdl2.boundary = find_boundary( mdl2.elems);

   % set gnd_node
   mdl2.gnd_node = nmap(mdl3.gnd_node);

   idx2 = []; % Don't know how to manage edges accurately

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
    if elecs(i).model == 'pem'
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
      [jnk, ind] = max(ang);
      electrode(i).nodes = electrode(i).nodes(ind);
    end
  end

