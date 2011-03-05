function mapping = mk_c2f_circ_mapping( mdl, xyzr );
% MK_C2F_CIRC_MAPPING: create a mapping matrix from circles/spheres to FEM
% mapping= mk_c2f_circ_mapping( mdl, xyzr );
%
% Mapping approximates elem_data_fine from elem_data_coase
%   elem_data_model = Mapping * elem_data_circles
%
% mdl is coarse fwd_model
% xyzr is the 3xN matrix (2D) or 4xN matrix (3D) of
%      circle centres and radii
%
% this function approximates using points interpolated into elements
%   use mdl.interp_mesh.n_points to control interpolation density  
%
% if a 3xN matrix is specified for a 3D model, then cylindrical
%  shapes (circle extruded in z) are selected
%

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end

c_obj = cache_obj(mdl, xyzr);

mapping = eidors_obj('get-cache', c_obj, 'circle_mapping');
if ~isempty(mapping)
    eidors_msg('mk_c2f_circ_mapping: using cached value', 3);
else

    mdl = fix_model(mdl);
    switch size(xyzr,1)
      case 3; mapping = contained_elems_2d( mdl, xyzr );
      case 4; mapping = contained_elems_3d( mdl, xyzr );
      case 5: error('size of xyzr incorrect');
    end

    eidors_obj('set-cache', c_obj, 'circle_mapping', mapping);
    eidors_msg('mk_coarse_fine_mapping: setting cached value', 3);
end

% Mapping depends only on nodes and elems - remove the other stuff
function c_obj = cache_obj(mdl, xyzr)
   c_obj = {mdl.nodes, mdl.elems, xyzr};

function mapping = contained_elems_2d_new( mdl, xyr );
   % We fill sparse by columns, (ie adding in CCS storage)
   Nc = size(xyr,      2); % Num circs
   too_far = elems_too_far( mdl, xyr );

   mapping = sparse( num_elems(mdl) , Nc );
   for i=1:Nc
     mapping(:,i) = circ_in_elem_2d(mdl, find( ~too_far(:,i)), ...
                xyr(1,i), xyr(2,i), xyr(3,i));
   end

% look only at elements 'look'
function mapping = circ_in_elem_2d( mdl, look, xc, yc, rc);
   Nt = elem_dim(mdl) + 1; % nodes per simplex
   pirc2 = pi*rc^2;
% start assuming no content
   mapping = sparse(num_elems(mdl),1);
% For each element, find nodes inside
   els = mdl.elems(look,:);
   ndx = reshape(mdl.nodes(els,1) - xc, size(els));
   ndy = reshape(mdl.nodes(els,2) - yc, size(els));
   n_in = (ndx.^2 + ndy.^2) < rc^2; % node inside
% triangles with 3 nodes inside are all in, stop looking at them
   all_n_in = sum(n_in,2) == Nt;
   mapping(look(all_n_in)) = 1;
   look(all_n_in)= []; n_in(all_n_in,:)= [];
% find distance inside each face
   f_in = zeros( length(look), Nt); 
   k=1;  for i= look(:)';
      faces = mdl.elem2face(i,:);
      out   =~mdl.inner_normal(i,:);
      f_norm= mdl.normals( faces, :);
      f_norm(out,:) = -f_norm(out,:);

      f_ctr = mdl.face_centre( faces,:);
      v_ctr = repmat([xc,yc],Nt,1) - f_ctr;
      v_ctr = sum(v_ctr .* f_norm,2)/rc; % >1 in, <-1 out, 
      f_in(k,:) = v_ctr';
   k=k+1;end

% triangles with any sides outside are out
   any_s_out= any(f_in<-1,2);
   look(any_s_out)= [];
   n_in(any_s_out,:) = [];
   f_in(any_s_out,:) = [];

% triangles with 3 sides inside are all in.
   all_s_in = sum(f_in>1,2) == Nt;
   mapping(look(all_s_in)) = pirc2 / ...
            mdl.elem_volume(look(all_s_in));
   look(all_s_in)= [];
   n_in(all_s_in) = [];
   f_in(all_s_in) = [];

% Now, all triangles should be partially in
% Calculate area chopped out 
   fin1 = f_in<1;
   a_out = zeros(size(fin1));
   a_out(fin1) = acos(f_in(fin1));
   a_out(fin1) = (a_out(fin1) - cos(a_out(fin1)).*sin(a_out(fin1)))/pi;
   
   % start with the default. This is accurate if there are
   % no contained nodes, otherwise we need to add back
   % those fractions
   mapping(look) = pirc2 / mdl.elem_volume(look);

   k=1; for i= look(:)';
      vol = pi*rc^2 / mdl.elem_volume(i);
      switch sum(n_in(k,:))
         case 0; % already do
   
         case 1; 

         case 2; 

         otherwise; error('cant get here'); 
      end
   k=k+1; end
   
   keyboard
     

%function mapping = contained_elems_2d_old( mdl, xyr );
 function mapping = contained_elems_2d( mdl, xyr );
   Ne = size(mdl.elems,1); % Num elems
   Nc = size(xyr,      2); % Num circs
   % We fill sparse by columns, due to CCS storage, this is fairly efficient
   mapping = sparse( Ne, Nc );

   % Interpolate
   n_interp =  7-size(mdl.nodes,2);
   m_pts = interp_mesh( mdl, n_interp); 
   for i=1:Nc
     xc = m_pts(:,1,:) - xyr(1,i);
     yc = m_pts(:,2,:) - xyr(2,i);
     inr= xc.^2 + yc.^2 < xyr(3,i)^2;
     frac= mean(inr,3);
     mapping(:,i) = frac;
   end

   % 1. Get furthest node in each element
   % 2. Get the max edge length of each elem
   % 3. Find elems that are too far
function too_far = elems_too_far( mdl, xyr );
   Ne = num_elems(mdl);
   Nc = size(xyr, 2); % Num circs
   Nt = elem_dim(mdl) + 1; % Elements per simplex
   nodes = repmat(mdl.nodes,[1,1,Nc]);
   targets = repmat(xyr(1:mdl_dim(mdl),:),[1,1,num_nodes(mdl)]);
   targets = shiftdim(targets,2);
   dist = nodes - targets;
   dist = sqrt(sum(dist.^2,2));
   node_target_dist = squeeze(dist);
   furthest_elem_node_dist = node_target_dist(mdl.elems,:);
   furthest_elem_node_dist = reshape(furthest_elem_node_dist,Ne,Nt,Nc);
   [furthest_elem_node_dist, furthest_elem_node]= max(furthest_elem_node_dist,[],2);
   furthest_elem_node_dist = squeeze(furthest_elem_node_dist);
   furthest_elem_node = squeeze(furthest_elem_node);
   
   max_edge_len =  repmat(mdl.max_edge_len,1,Nc);
   radius = ones(Ne,1)*xyr(Nt,:);
   too_far = (furthest_elem_node_dist - max_edge_len) > radius;

function mapping = contained_elems_3d( mdl, xyr );
   Ne = size(mdl.elems,1); % Num elems
   Nc = size(xyr,      2); % Num circs
   % We fill sparse by columns, due to CCS storage, this is fairly efficient
   mapping = sparse( Ne, Nc );

%    % INterpolate
%    n_interp = 4; % 7-df
%    m_pts = interp_mesh( mdl, n_interp); 
%    for i=1:Nc
%      mapping(:,i) = contained_elem_pts(m_pts, xyr(:,i));
%    end

   % 4. Make a tmp model with only the remaining elems
   % 5. Interpolate
   % 6. Merge

   too_far = elems_too_far( mdl, xyr );
   
   tmp = eidors_obj('fwd_model','tmp','nodes',mdl.nodes,'elems',mdl.elems);
   mapping = sparse( Ne, Nc );
   n_interp = 6;
   for i=1:Nc
       good = ~too_far(:,i);
       tmp.elems = mdl.elems(good,:);
       m_pts = interp_mesh( tmp, n_interp);
       mapping(good,i) = contained_elem_pts(m_pts, xyr(:,i));
   end
   
   
   
function frac= contained_elem_pts(m_pts, xyr);
% This is more clear
%    xc = m_pts(:,1,:) - xyr(1);
%    yc = m_pts(:,2,:) - xyr(2);
%    zc = m_pts(:,3,:) - xyr(3);
%    inr= xc.^2 + yc.^2 + zc.^2 < xyr(4)^2;

% But this is how to stop matlab from wasting memory
     inr = (m_pts(:,1,:) - xyr(1)).^2 + ...% xc =
           (m_pts(:,2,:) - xyr(2)).^2 + ...% yc =
           (m_pts(:,3,:) - xyr(3)).^2;     % zc =
     inpts = inr < xyr(4)^2;
     frac= mean( inpts ,3);
     if sum(inpts(:))==0
         eidors_msg(['mk_c2f_circ_mapping: Interpolation failed: increase ', ...
                         'fwd_model.interp_mesh.n_interp']);
     end

function do_unit_test
   %2D example
   imdl = mk_common_model('a2c2',16); fmdl=imdl.fwd_model;
   xyc = [0,0.27,0.18;0,-0.1,0.03;0,0.1,0.2;0.1,0.37,0.1]';
   th=linspace(0,2*pi,20)';
   xx=[0*th+1]*xyc(1,:)+sin(th)*xyc(3,:);
   yy=[0*th+1]*xyc(2,:)+cos(th)*xyc(3,:);
   show_fem(fmdl,[0,0,1]); set(line(xx,yy),'Color',[1,0,0],'LineWidth',2);

   c2f= mk_c2f_circ_mapping( fmdl, xyc );
   t1= all( abs(c2f(1:4)-0.2857)<.001 ) & all( c2f(5:end)==0 );
   fprintf('2D ex 1: [%d]\n',full(t1));

   c2f= mk_c2f_circ_mapping( fmdl, [.0;.05;0.03]); 
   t2= abs( c2f(1) - 0.1429) < .001 & all( c2f(2:end)==0 );
   fprintf('2D ex 2: [%d]\n',full(t2));
      
   %3D example - cylinder
   imdl = mk_common_model('a3cr',16); fmdl=imdl.fwd_model;
   c2f= mk_c2f_circ_mapping( fmdl, [0;0;0.1]); 
   t3= all( abs(c2f(1:4)-0.1714)<.001 ) & all( c2f(5:64)==0 );
   fprintf('3D ex 1: [%d]\n',full(t3));

   %3D example - cylinder
   imdl = mk_common_model('a3cr',16); fmdl=imdl.fwd_model;
   c2f= mk_c2f_circ_mapping( fmdl, [0;0;0;0.1]); 
   t4= all( abs(c2f(193:196)-0.0595)<.001 ) & all( c2f(1:64)==0 );
   fprintf('3D ex 2: [%d]\n',full(t4));
   
   %3D example - cylinder - 2 pts
   imdl = mk_common_model('a3cr',16); fmdl=imdl.fwd_model;
   c2f= mk_c2f_circ_mapping( fmdl, [0 0;0 0;0 0;0.1 0.2]); 
   t4= all( abs(c2f(193:196,1)-0.0595)<.001 ) & all( c2f(1:64)==0 );
   fprintf('3D ex 3: [%d]\n',full(t4));
